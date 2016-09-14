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

class description {
public:
	description(unsigned int& count):_count(count){};
	void operator()(Reaction* r){
		cout << _count << "\t" << *r << endl;_count++;
	}
	unsigned count() {return _count;}
private:
	unsigned int& _count;
};


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
	bool lsc_moved = r->apply(_pool, prev_t);
	if(lsc_moved & (_lsctime == -1)) {
		_lsctime = prev_t;
	}
	//4. update every reaction dependend upon the reaction that was just used.
	//first set a_{i} the current reaction;
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
	_queue.update(next->queueLoc()); //verify whether heapify is correctly executed

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

void Kernel::adjustReactions(RanGen& ran, unsigned compartment, unsigned type){
	unsigned which[]={C,B}; //C = Cancer cells, B=Treated cancer cells
	cout << "RESET REACTIONS " << endl;
	for(int i=0; i < 2; ++i){
		DependencyNode* node = _depend.get(type, ((compartment-1)*4)+which[i]);
		Reaction* r = _allr[node->reaction()];
		QueueElement* elm = _queue[node->reaction()];
		double prev_t = _time; 
		double next_estimate =numeric_limits<double>::infinity();
		cout << setprecision(8) << *r << "\t" <<  _time << " \t" << elm->tau() << "\t" << elm->queueLoc() << endl;
	
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
		cout << setprecision(8)<< *r << "\t" << prev_t << "\t" << elm->tau() << "\t" << elm->queueLoc() << endl;
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
	for(unsigned k = _pool.numStoch(); k < _pool.numComp(); ++k){
			_pool.updateDet(k,_data);
	}
}

double Kernel::execute(RanGen& ran, double t, bool treat){
	_time = (t*365.0); // _time is in the function expressed in days
	double  _time_step = _data.dt();
//	cout << "##execute starts " << _time << endl;
	int iters = (int)ceil(_time / _data.dt());
	int endsim = _data.ntimes();

	if(treat) {
		endsim = iters + (int) (_data.treatment()*365.0/_data.dt());	
//		cout << "treatment ends at " <<  treatmentend << " iterations ("<< (t + _data.treatment()) << " years)" << endl;
	}
	double next_stoch = (_queue.top())->tau(); //when occurs the next stochastic reaction
	
        std::cout <<"debug before: "<<_pool.diagnosis(_data)<<" "<<_pool.lastN()<<" "<<_pool.containsLSC()<<" "<<_pool.diseaseBurden()<<std::endl;
	while(iters < endsim && ( (!treat && !_pool.diagnosis(_data)) || (treat && !_pool.reduction(_data)) )){
		
		//treat cells -> affects reactions in priorityqueue !!
		if(treat) {
			treatCells(ran); // treat fraction of cells
			next_stoch = (_queue.top())->tau();
		}

		//start new update
		_pool.memorize(); //every time we update the state is stored (calculations are performed on these states)
                _pool.check_LSCvanished(_time);
                // if (_time/365. >= yearscounter ){
                //     _pool.print_cells(std::cout,_time);
                //     ++yearscounter;
                // }
		while(_time >= next_stoch)	{
			nextMethod(ran);
			next_stoch = (_queue.top())->tau();
		}
		detUpdate();
		
		iters++;
		_time += _time_step; 
	}
	if(!treat) {
		_pool.calcAlpha(); // required to recalculate disease burden
		_pool.setDiagRes((_time/365.0));
	}
	else {
		_pool.setWhenReduction((_time/365.0));
	}
//	cout << "laste iter  " << iters << "\t  " << (_time/365.0) << endl;

        std::cout <<"debug after: "<<_pool.diagnosis(_data)<<" "<<_pool.lastN()<<" "<<_pool.containsLSC()<<" "<<_pool.diseaseBurden()<<std::endl;
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




