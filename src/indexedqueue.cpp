/*
 *  indexedqueue.cpp
 *  stochmut
 *
 *  Created by Tom Lenaerts on 2/2/11.
 *  Copyright 2011 Universite Libre de Bruxelles. All rights reserved.
 *
 */

#include "indexedqueue.h"
#include <limits>

QueueElement::QueueElement(const QueueElement& other){
	_idx=other.index();
	_tau=other.tau();
	_queloc=other.queueLoc();
}

QueueElement& QueueElement::operator=(const QueueElement& other){
	_idx=other.index();
	_tau=other.tau();
	_queloc=other.queueLoc();
	return *this;
}

std::ostream& QueueElement::display(std::ostream& os){
	os << "["<< _idx << " " <<_tau << " " << _queloc << "]";
	return os;
}



IndexedQueue::IndexedQueue(RanGen& ran, Model& pool,AllReactions& all, double t, double dt):_pqueue(&_compare, &_swap),_sentinel(NULL){
    init(ran,all);
}

void IndexedQueue::init(RanGen& ran, AllReactions& all){

	_pqueue.clear();
	for(unsigned i=0; i < _indices.size() ; ++i){
		QueueElement* elm = _indices[i];
		delete elm;
		_indices[i] = NULL;
	}
        _indices.clear();


	//create empty PQueue using information in all.
	if(!_sentinel){
		_sentinel = new QueueElement();
		_pqueue.setSentinel(_sentinel);  // sentinel is not added to the indexvector
	}
        else if (_pqueue.size()==0){//if init() is called again, sentinel has to be set again!
            _pqueue.setSentinel(_sentinel);
        }
	
	for (unsigned r=0; r < all.size() ; ++r){
		// cout << r << "\t" << all[r]->propensity() << "\t";
		double time_i = std::numeric_limits<double>::infinity();
		if(all[r]->propensity() > 0.0){
			double rval = ran.randouble();
			time_i = all[r]->calcPutativeTime(rval); 
			// cout << "Rval = " << rval << " for reaction " << *(all[r]) << " ";//" produces pututative time = "<<  time_i << endl;
		}
		else all[r]->setPutativeTime(time_i);
		// cout << all[r]->putativeTime()<<"\t" << endl;
		
		QueueElement *elm = new QueueElement(r, time_i);
		_indices.push_back(elm);
		elm->setQueueLoc(_pqueue.size()); //starts at 1!!!
		try {
			_pqueue.push(elm);
		}
		catch (std::exception & e){
                    std::cout <<" HEAP exception: " << e.what() << std::endl;
			exit(0);
		}
	}
        // all.print(std::cout);
}

IndexedQueue::~IndexedQueue(){
	_pqueue.clear();
	for(unsigned i=0; i < _indices.size() ; ++i){
		QueueElement* elm = _indices[i];
		delete elm;
		_indices[i] = NULL;
	}
	delete _sentinel;
}

