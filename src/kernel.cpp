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
    std::cout << endl << "#pool " << endl;
    std::cout << _pool<< endl;
    std::cout << "#reactions " << endl;
    for(unsigned r = 0; r < _allr.size() ; ++r){
        std::cout << r << "\t" << *_allr[r] << endl;
    }
    std::cout << "#Number of reactions " << _allr.size() << endl;
    std::cout << "#Sum propensity " << _allr.propSum() << endl;
    std::cout << endl << _depend << endl;
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
        std::cout <<"further diag time: " <<getDiagnosisTime()*365.<<" "<<_pool.get(r->reactantComp(),r->reactant())<< " <cs in 0comp>="<<_pool.getC(0)<<" <Cs in 1comp>="<<_pool.getC(1)<<std::endl;
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
        double next_estimate =numeric_limits<double>::infinity();

        if(!rd->sufficientReactants(_pool) && qed->tau() < numeric_limits<double>::infinity()){
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

double Kernel::execute(RanGen& ran, double t, bool treat){

    _time = (t*365.0); // _time is in the function expressed in days
    double t_max=_data.getTmax_in_years()*365.;
    double  _time_step = _data.dt();
    if(treat) {
        t_max=_time+_data.treatment_dur()*365.;
    }

    //######turn treatment on or off
    if (treat) 
        _pool.setTreatRate(_data.treatment_rate());
    else
        _pool.setTreatRate(0.);

    for (unsigned int r_id=0; r_id< _allr.size(); ++r_id){
        if (_allr[r_id]->inType()==3){
            _allr[r_id]->setRate(_pool.getTreatRate());
        }
    }
    reinitialize(_pool,ran,_time);
    //#####done switching treatment

    int iters =0.;// (int)ceil(_time / _data.dt());

    double next_stoch = (_queue.top())->tau(); //when occurs the next stochastic reaction

    // while(_time<t_max && ( (!treat && !_pool.diagnosis(_data)) || (treat && !_pool.reduction(_data)) )){
    while(_time<t_max &&  (treat || (!treat && !_pool.diagnosis(_data)))){

        //start new update
        _pool.memorize(); //every time we update the state is stored (calculations are performed on these states)
        _pool.check_LSCvanished(_time);
        _time += _time_step; 
        int reactions_count=0;
        double starttime_reacts=next_stoch;
        while(_time >= next_stoch)	{
            nextMethod(ran);
            next_stoch = (_queue.top())->tau();
            ++reactions_count;
        }
        detUpdate();
        iters++;
    }
    if(!treat) {
        _pool.calcAlpha(); // required to recalculate disease burden
        _pool.setDiagRes((_time/365.0));
    }
    else {
        _pool.setWhenReduction((_time/365.0));
    }

    return ( _time / 365.0);

}

float Kernel::burden(){
    return _pool.diseaseBurden();
}

bool Kernel::reachedDiagnosis(){
    return _pool.diagnosis(_data);
}

float Kernel::getDiagnosisTime(){
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
    reinitialize(_pool,ran,t*365.);

}



