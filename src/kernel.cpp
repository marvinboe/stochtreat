/*
 *  kernel.cpp
 *  stochmut
 *
 *  Created by Tom Lenaerts on 2/5/11.
 *  Copyright 2011 Universite Libre de Bruxelles. All rights reserved.
 *
 */

#include "kernel.h"
#include <algorithm>
#include <limits>
#include <iomanip>



void Kernel::printAll(){
    cout << endl << "#pool " << endl;
    cout << _pool<< endl;
    cout << "#reactions " << endl;
    for(unsigned r = 0; r < _allr.size() ; ++r){
        cout << r << "\t" << *_allr[r] << endl;
    }
    cout << "#Number of reactions " << _allr.size() << endl;
    cout << "#Sum propensity " << _allr.propSum() << endl;
    cout << endl << _depend << endl;
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
    catch (...) {
        std::cout <<"error reaction without reactants: debug output <tau> <index> <reactantType> <reactantComp> <reactionType> <queueindex> <propensity>"<<std::endl;
        std::cout <<next->tau()<<" "<<next->index()<<" " <<r->reactant()<<" "<<r->reactantComp()<<" "
            <<r->inType()<<" "<<r->inGraph()<<" "<<r->propensity()<<std::endl;
        std::cout <<"further diag time: " <<getDiagnosis()*365.<<" "<<_pool.get(r->reactantComp(),r->reactant())<< " <cs in 0comp>="<<_pool.getC(0)<<" <Cs in 1comp>="<<_pool.getC(1)<<std::endl;
        exit(1);
    }
    if(lsc_moved & (_lsctime == -1)) {
        _lsctime = prev_t;
    }
    //4. update every reaction dependend upon the reaction that was just used.
    //first set a_{i} of the current reaction;
    double next_estimate =numeric_limits<double>::infinity();
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
    vector<DependencyNode*>::iterator start = node->begin();
    vector<DependencyNode*>::iterator stop = node->end();
    while (start != stop){
        if((*start)->reaction() != next->index()){
            next_estimate =numeric_limits<double>::infinity();
            Reaction *rd = _allr[(*start)->reaction()];
            QueueElement *qed = _queue[(*start)->reaction()];
            if(!rd->sufficientReactants(_pool) && qed->tau() < numeric_limits<double>::infinity()){
                if(rd->getTZero() == -1.0){ // if it is the first time that the reactants are not sufficient
                    rd->setLasta();
                    rd->setTZero(prev_t);
                    rd->setLastPtime();
                    rd->setPropensity(0.0);  // only now is allowed
                }
            }
            else if(rd->sufficientReactants(_pool)) {//TODO unnesseccary "if"?
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
                if(rd->getTZero() < 0.0 && qed->tau() == numeric_limits<double>::infinity()){ 
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

void Kernel::reinitialize(Model& pool,RanGen& ran,double prev_t){

    for (unsigned int r_id =0 ; r_id < _allr.size(); ++r_id){
        Reaction* rd=_allr[r_id];
        QueueElement* qed=_queue[r_id];
        // std::cout <<"reinitialize "<<prev_t<<" "<<r_id<<" "<<_allr[r_id]->propensity()<<" " <<_allr[r_id]->rate()<<" "<<qed->tau()<<" "<<_allr[r_id]->reactantFactor(pool)<<" "; 
        double next_estimate =numeric_limits<double>::infinity();

        if(!rd->sufficientReactants(_pool) && qed->tau() < numeric_limits<double>::infinity()){
            if(rd->getTZero() == -1.0){ // if it is the first time that the reactants are not sufficient
                rd->setLasta();
                rd->setTZero(prev_t);
                rd->setLastPtime();
                rd->setPropensity(0.0);  // only now is allowed TODO?
            }
        }
        else if(rd->sufficientReactants(_pool)) {//TODO unnesseccary "if"?
            double prevtau = qed->tau(); // can also be infinity
            double t1 = prev_t;
            double t2 = prev_t;
            double a_old = rd->propensity();		
            // if(rd->getTZero() >= 0.0){ //if the reactants were zero before and qed->tau = infinity
            //     a_old = rd->last_a();
            //     t1 = rd->getTZero();
            //     prevtau = rd->last_ptime();
            //     rd->setTZero(-1.0);
            // }	
            rd->setPropensity(rd->reactantFactor(_pool));
            //queuelement should be  in the same position in the indexedqueue as in the reaction pool
                //rd->getTZero() < 0.0 && //rection has never been used before
            if(qed->tau() == numeric_limits<double>::infinity()){ 
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

        // std::cout <<"after "<< setprecision(16) <<rd->propensity()<<" " <<rd->rate() <<" "<<qed->tau()<<std::endl; 
    }
    // std::cout <<"############## first time: "<<_queue.top()->tau()<<" "<<prev_t<<std::endl;

}

void Kernel::adjustReactions(RanGen& ran, unsigned compartment, unsigned type){
    unsigned which[]={C,B}; //C = Cancer cells, B=Treated cancer cells
    // cout << "RESET REACTIONS " << endl;
    for(int i=0; i < 2; ++i){
        DependencyNode* node = _depend.get(type, ((compartment-1)*4)+which[i]);
        Reaction* r = _allr[node->reaction()];
        QueueElement* elm = _queue[node->reaction()];
        double prev_t = _time; 
        double next_estimate =numeric_limits<double>::infinity();
        // cout << setprecision(8) << *r << "\t" <<  _time << " \t" << elm->tau() << "\t" << elm->queueLoc() << endl;

        if(!r->sufficientReactants(_pool) && elm->tau() < numeric_limits<double>::infinity()){
            if(r->getTZero() == -1.0){ // if it is the first time that the reactants are not sufficient
                r->setLasta();
                r->setTZero(prev_t);
                r->setLastPtime();
                r->setPropensity(0.0);  // only now is allowed
            }
        }
        else if(r->sufficientReactants(_pool)) {  
            double prevtau = elm->tau(); // can also be infinity
            double t1 = prev_t;
            double t2 = prev_t;
            double a_old = r->propensity();		
            if(r->getTZero() >= 0.0){ //if the reactants were zero before and qed->tau = infinity
                a_old = r->last_a();
                t1 = r->getTZero();
                prevtau = r->last_ptime();
                r->setTZero(-1.0);
            }	
            r->setPropensity(r->reactantFactor(_pool));

            //queuelement should be  in the same position in the indexedqueue as in the reaction pool
            if(r->getTZero() < 0.0 && elm->tau() == numeric_limits<double>::infinity()){ 
                //rection has never been used before
                next_estimate = r->calcPutativeTime(ran.randouble(), prev_t);
            }
            else next_estimate = (a_old / r->propensity()) * (prevtau - t1) + t2;	
            r->setPutativeTime(next_estimate);

        }


        //adjust the queue, find the location of the reactin in the queue and change and update tau
        elm->setTau(next_estimate); 
        _queue.update(elm->queueLoc()); 	
        // cout << setprecision(8)<< *r << "\t" << prev_t << "\t" << elm->tau() << "\t" << elm->queueLoc() << endl;
    }
}

void Kernel::treatCells(RanGen& ran){
    //difference between stochastic and deterministic treatment !

    //calculate probability of treatment, which is rate x time
    float tobetreated=(_data.dt() * _data.treatment_rate());
    for(unsigned k = 1; k < _pool.numStoch(); ++k){
        bool changed = _pool.treatStochastically(k,tobetreated, ran);
        //adjust the selfrenewal and differentiation reactions accordingly !
        if(changed){
            adjustReactions(ran, k, SELF_RENEWAL);
            adjustReactions(ran, k, DIFFERENTATION);
        }
    }	
    for(unsigned k = _pool.numStoch(); k < _pool.numComp(); ++k){
        _pool.treatDeterministically(k,tobetreated);
    }	
}

void Kernel::detUpdate(){
    //update the deterministic compartments using the standard formula's
    double treatamount=(_data.dt() * _pool.getTreatRate());

    for(unsigned k = _pool.numStoch(); k < _pool.numComp(); ++k){
        _pool.updateDet(k,_data);
        _pool.treatDeterministically(k,treatamount);
    }

    // for(unsigned k = _pool.numStoch(); k < _pool.numComp(); ++k){
    // }	
}

double Kernel::execute(RanGen& ran, double t, bool treat){

    _time = (t*365.0); // _time is in the function expressed in days
    double t_max=_data.getTmax_in_years()*365.;
    double  _time_step = _data.dt();
    if(treat) {
        t_max=_time+_data.treatment_dur()*365.;
    }
    // std::cout <<"#timing debug: "<<_time<<" "<<t_max<<" "<<_data.getTmax_in_years()<<" "<<_data.treatment_dur()<<std::endl;

    if (treat) 
        _pool.setTreatRate(_data.treatment_rate());
    else
        _pool.setTreatRate(0.);

    //turn treatment on or off
    for (unsigned int r_id=0; r_id< _allr.size(); ++r_id){
        if (_allr[r_id]->inType()==3){
            _allr[r_id]->setRate(_pool.getTreatRate());
        }
    }
    reinitialize(_pool,ran,_time);

    //	cout << "##execute starts " << _time << endl;
    int iters =0.;// (int)ceil(_time / _data.dt());

    double next_stoch = (_queue.top())->tau(); //when occurs the next stochastic reaction
    // while(_time<t_max && ( (!treat && !_pool.diagnosis(_data)) || (treat && !_pool.reduction(_data)) )){
    while(_time<t_max &&  (treat || (!treat && !_pool.diagnosis(_data)))){

        //start new update
        _pool.memorize(); //every time we update the state is stored (calculations are performed on these states)
        _pool.check_LSCvanished(_time);
        // if (_time/365. >= yearscounter ){
        //     _pool.print_cells(std::cout,_time);
        //     ++yearscounter;
        // }
        _time += _time_step; 
        int reactions_count=0;
        double starttime_reacts=next_stoch;
        while(_time >= next_stoch)	{
            nextMethod(ran);
            // std::cout <<std::setprecision(9)<<next_stoch<<" ";
            next_stoch = (_queue.top())->tau();
            ++reactions_count;
        }
        // std::cout <<std::endl;
        // std::cout <<"reaction debug: "<<reactions_count<<" "<<next_stoch-starttime_reacts<<" "<<_time<<std::endl;
        detUpdate();

        // std::cout <<"########debug output "<<_time<<std::endl;
        // for (int k=0; k<_data.ncompartments();++k){
        //     std::cout << k << " : "<< _pool.getH(k) << " " << _pool.getC(k) << " " << _pool.retrieveH(k) << " " << _pool.retrieveC(k) << endl;
        // }
        // std::cin.ignore();

        iters++;
    }
    if(!treat) {
        _pool.calcAlpha(); // required to recalculate disease burden
        _pool.setDiagRes((_time/365.0));
    }
    else {
        _pool.setWhenReduction((_time/365.0));
    }
    //	cout << "laste iter  " << iters << "\t  " << (_time/365.0) << endl;

    // std::cout <<"debug after: "<<_pool.diagnosis(_data)<<" "<<_pool.lastN()<<" "<<_pool.containsLSC()<<" "<<_pool.diseaseBurden()<<std::endl;
    return ( _time / 365.0);

}

float Kernel::burden(){
    return _pool.diseaseBurden();
}

bool Kernel::reachedDiagnosis(){
    return _pool.diagnosis(_data);
}

float Kernel::getDiagnosis(){
    return _pool.diagRes();
}


bool Kernel::reachedReduction(){
    return _pool.reduction(_data);
}

float Kernel::getReduction(){
    return _pool.getReduction();
}

float Kernel::whenReduction(){
    return _pool.when();
}


float Kernel::get_nolsctime(){
    return _pool.get_nolsctime();
}


bool Kernel::hasLSC(){
    return _pool.containsLSC();
}

void Kernel::addStochCompSizes(double* data){
    for (unsigned i=0; i < (_pool.numStoch()+1); ++i) {
        data[i] += _pool.getN(i);
    }
}

ostream& Kernel::writeModel(ostream& output){
    output << _pool;
    return output;
}

std::istream& Kernel::readModel(std::istream& input){
    input >> _pool;
    return input;
}

void Kernel::reset_treatment(RanGen& ran,double t){

    _pool.reset_treatment();
    reinitialize(_pool,ran,t);

}



