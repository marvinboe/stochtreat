/*
 *  kernel.cpp
 *  stochmut
 *
 *  Created by Tom Lenaerts on 2/5/11.
 *  Copyright 2011 Universite Libre de Bruxelles. All rights reserved.
 *
 */

#include "kernel.h"



void Kernel::printAll(){
    std::cout << std::endl << "#pool " << std::endl;
    std::cout << _pool<< std::endl;
    std::cout << "#reactions " << std::endl;
    for(unsigned r = 0; r < _allr.size() ; ++r){
        std::cout << r << "\t" << *_allr[r] << std::endl;
    }
    std::cout << "#Number of reactions " << _allr.size() << std::endl;
    std::cout << "#Sum propensity " << _allr.propSum() << std::endl;
    std::cout << std::endl << _depend << std::endl;
}

bool Kernel::directMethod(RanGen& ran){
    return true;
}

bool Kernel::nextMethod(RanGen& ran){
    // x. all a_{i} (done when dependency graph was determined)
    // y. all putative times tau_{i} are created when the indexed queue was constructed.
    // 1. the top element is the next reaction in the stochastic process.
    QueueElement* next = _queue.top();
    double prev_t = next->tau();
    Reaction* r = _allr[next->index()];

    //2. apply the reaction to the compartment/cells in the pool
    bool lsc_moved=false;
    try{
        lsc_moved = r->apply(_pool, prev_t);
    }
    catch (std::exception & e) {
        std::cout <<"reactions exception: '"<<e.what()<<"'"<<std::endl;
        std::cout<<"debug output <tau> <index> <reactantType> <reactantComp> <reactionType> <queueindex> <propensity>"<<std::endl;
        std::cout <<next->tau()<<" "<<next->index()<<" " <<r->reactant()<<" "<<r->reactantComp()<<" "
            <<r->inType()<<" "<<r->inGraph()<<" "<<r->propensity()<<std::endl;
        exit(1);
    }
    if((_lsctime < 0.)&&lsc_moved &&_pool.getC(0)==0 ) {
        _lsctime = prev_t;
    }
    //4. update every reaction dependend upon the reaction that was just used.
    //first set a_{i} of the current reaction;
    double next_estimate =std::numeric_limits<double>::infinity();
    if(r->sufficientReactants(_pool)){
        r->setPropensity(r->reactantFactor(_pool)); 
        next_estimate = r->calcPutativeTime(ran.randouble(), prev_t);
    }
    else { // store some information relevant for later calculations.
        if(r->getTZero() == -1.0){
            r->setLasta();
            r->setTZero(prev_t);
            r->setLastPtime();
            r->setPropensity(0.0);  // only now  allowed
        }
    }
    next->setTau(next_estimate); 
    _queue.update(next->queueLoc()); //repair the heap tree 

    //then set the a_{j} and tau_{j} of the dependent reactions
    unsigned type = r->inType();
    unsigned graph = r->inGraph();
    DependencyNode* node = _depend.get(type, graph);
    std::vector<DependencyNode*>::iterator start = node->begin();
    std::vector<DependencyNode*>::iterator stop = node->end();
    while (start != stop){
        if((*start)->reaction() != next->index()){
            next_estimate =std::numeric_limits<double>::infinity();
            Reaction *rd = _allr[(*start)->reaction()];
            QueueElement *qed = _queue[(*start)->reaction()];
            if(!rd->sufficientReactants(_pool) && qed->tau() < std::numeric_limits<double>::infinity()){
                if(rd->getTZero() == -1.0){ // if it is the first time that the reactants are not sufficient
                    rd->setLasta();
                    rd->setTZero(prev_t);
                    rd->setLastPtime();
                    rd->setPropensity(0.0);  // only now is allowed
                }
            }
            else if(rd->sufficientReactants(_pool)) {
                double prevtau = qed->tau(); // can also be infinity
                double t1 = prev_t;
                double t2 = prev_t;
                double a_old = rd->propensity();		
                if(rd->getTZero() >= 0.0){ //if the reactants were zero before and qed->tau = infinity
                    a_old = rd->last_a();
                    t1 = rd->getTZero();
                    prevtau = rd->last_ptime();
                    rd->setTZero(-1.0);
                }	
                rd->setPropensity(rd->reactantFactor(_pool));
                //queuelement should be  in the same position in the indexedqueue as in the reaction pool
                if(rd->getTZero() < 0.0 && qed->tau() == std::numeric_limits<double>::infinity()){ 
                    //rection has never been used before
                    next_estimate = rd->calcPutativeTime(ran.randouble(), prev_t);
                }
                else next_estimate = (a_old / rd->propensity()) * (prevtau - t1) + t2;	
                rd->setPutativeTime(next_estimate);

            }

            qed->setTau(next_estimate); 
            _queue.update(qed->queueLoc());
        }
        ++start;
    }
    return lsc_moved;
}

void Kernel::reinitialize(RanGen& ran,double prev_t){

    for (unsigned int r_id =0 ; r_id < _allr.size(); ++r_id){
        Reaction* rd=_allr[r_id];
        QueueElement* qed=_queue[r_id];
        double next_estimate =std::numeric_limits<double>::infinity();

        if(!rd->sufficientReactants(_pool) && qed->tau() < std::numeric_limits<double>::infinity()){
            if(rd->getTZero() == -1.0){ // if this is the first time that the reactants are not sufficient
                rd->setLasta();
                rd->setTZero(prev_t);
                rd->setLastPtime();
                rd->setPropensity(0.0); 
            }
        }
        else if(rd->sufficientReactants(_pool)) {
            double prevtau = qed->tau(); // can also be infinity
            double t1 = prev_t;
            double t2 = prev_t;
            double a_old = rd->propensity();		
            rd->setPropensity(rd->reactantFactor(_pool));

            if(qed->tau() == std::numeric_limits<double>::infinity()){ 
                if(rd->propensity() > 0.0){
                    next_estimate = rd->calcPutativeTime(ran.randouble(), prev_t);
                }
            }
            else 
                next_estimate = (a_old / rd->propensity()) * (prevtau - t1) + t2;	

            rd->setPutativeTime(next_estimate);
        }

        qed->setTau(next_estimate); 
        _queue.update(qed->queueLoc());
    }
}

void Kernel::detUpdate(){
    //update the deterministic compartments using the standard formula's
    double treatamount=(_data.dt() * _pool.getTreatRate());

    for(unsigned k = _pool.numStoch(); k < _pool.numComp(); ++k){
        _pool.updateDet(k,_data);
        _pool.treatDeterministically(k,treatamount);
    }
}

double Kernel::execute(RanGen& ran, double t, int sim_type){

    _time = (t*365.0); // _time is in the function expressed in days
    double t_max=_data.getTmax_in_years()*365.;
    double  _time_step = _data.dt();
    if(sim_type==TREATMENTRUN) {
        t_max=_time+_data.treatment_dur()*365.;
    }


    //######turn treatment on or off
    if (sim_type==TREATMENTRUN) {
        _pool.setTreatRate(_data.treatment_rate());
        _doctor.calc_initial_reference(_time,_pool);
    }
    else
        _pool.setTreatRate(0.);

    for (unsigned int r_id=0; r_id< _allr.size(); ++r_id){
        if (_allr[r_id]->inType()==3){
            _allr[r_id]->setRate(_pool.getTreatRate());
        }
    }
    reinitialize(ran,_time);
    //#####done switching treatment
    
    int iters =0.;// (int)ceil(_time / _data.dt());

    double next_stoch = (_queue.top())->tau(); //when occurs the next stochastic reaction

    while(_time<t_max && !stopsim(_time,sim_type) ){

        _doctor.consult(_time,_pool);
        //start new update
        _pool.memorize(); //every time we update the state is stored (calculations are performed on these states)
        _time += _time_step; 
        int reactions_count=0;
        while(_time >= next_stoch)	{
            nextMethod(ran);
            next_stoch = (_queue.top())->tau();
            ++reactions_count;
        }
        detUpdate();
        iters++;
    }
    _doctor.consult(_time,_pool);
    if(sim_type==DIAGNOSISRUN) {
        _pool.setDiagRes((_time/365.0));
    }

    return ( _time / 365.0);

}

float Kernel::get_nolsctime() const{
    return _lsctime/365.;
}


bool Kernel::hasLSC() const{
    return _pool.containsLSC();
}

void Kernel::addStochCompSizes(std::vector<double>& data) const{
    for (unsigned i=0; i < data.size(); ++i) {
        data[i] += _pool.getN(i);
    }
}

std::ostream& Kernel::writeModel(std::ostream& output){
    output << _pool;
    return output;
}

std::istream& Kernel::readModel(std::istream& input){
    input >> _pool;
    return input;
}

void Kernel::reset_treatment(RanGen& ran,double t){

    _pool.reset_treatment();
    reinitialize(ran,t*365.);

}


void Kernel::introduce_immunity_inlowest(){

    unsigned int k =0;
    bool cancer_cell_replaced=false;
    while (k< _data.ncompartments() &&  !cancer_cell_replaced){
        cancer_cell_replaced=_pool.manual_mutation(k);
        ++k;
    }


}

void Kernel::introduce_resistance(unsigned k){
    if (k==0){
        if (_pool.getC(0)>0)
            _pool.decr(0,C,1.);
        else 
            _pool.decr(0,H,1.);
        _pool.incr(0,I,1.);
    }
    else{
        _pool.incr(k,I,1.);
    }

}

bool Kernel::stopsim(double t,int sim_type){
    if (sim_type==TREATMENTRUN){
        if (_doctor.reduction_reached()){
            if (_stoptimer<=0.){
                _stoptimer=t;
            }
            if (t >= _stoptimer+_data.required_reduction_time()){
                _stoptimer=-1.;
                return true;
            }
        }
        else{
            // _stoptimer=-1.;
            return false;
        }
    }
    else if (sim_type==DIAGNOSISRUN){
        if (_doctor.diagnosis_reached()) return true;
        else return false;
    }
    else if (sim_type==RELAPSERUN){
        if (_doctor.relapse_reached()) return true;
        else return false;
    }
    return false;
}

bool Kernel::read_model(std::string path, int runid, int i){
    // writeModel(std::cout);

    //#################read compartment data from file##########
    std::stringstream ssin;
    ssin << path<< "patient-"<< runid << "-"<< i << ".txt";
    std::ifstream input(ssin.str().c_str());
    if(!input.is_open()){
        return false;
    }
    else{//input is open
        readModel(input);
        // ker.writeModel(std::cout);
        input.close();
    }
    return true;
    //#############end reading model data ######################
}

bool Kernel::write_model(std::string path, int runid, int i){

    //##########write compartment data to file
    std::stringstream ss;
    ss << path<< "patient-"<< runid << "-"<< i << ".txt";
    std::ofstream output(ss.str().c_str());
    if(!output.is_open()){
        std::cout << " error: unable to open output file " << ss.str() << std::endl;
        throw std::runtime_error("unable to open output");
    }
    writeModel(output);
    output.close();
    //######end writing patient data to file
    return true;
}
