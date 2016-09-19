/*
 *  indexedqueue.h
 *  stochmut
 *
 *  Created by Tom Lenaerts on 2/2/11.
 *  Copyright 2011 Universite Libre de Bruxelles. All rights reserved.
 *
 */

#ifndef __INDEXEDQUEUE_H
#define __INDEXEDQUEUE_H

#include <string>
#include <istream>
#include <iostream>
#include <iomanip>
#include <vector>
#include "PriorityQueue.h"
#include "reactions.h"
#include "model.h"
#include "rangen.h"

using  namespace std;



class QueueElement {
public:
	QueueElement():_idx(999),_tau(-1.0){};
	QueueElement(unsigned idx, double tau):_idx(idx), _tau(tau){};
	QueueElement(const QueueElement& other);
	
	QueueElement& operator=(const QueueElement& other);
	
	unsigned index() const {return _idx;}
	void setIndex(unsigned idx) {_idx=idx;}
	double tau() const {return _tau;}
	void setTau(double tau) {_tau=tau;}
	unsigned queueLoc() const {return _queloc;}
        /** set location of queue element to l*/
	void setQueueLoc(unsigned l) { _queloc = l;}
	
	
	void print() {cout <<  "["<< _idx << " " <<_tau << " " << _queloc << "]"<< endl;}

	friend ostream & operator<<(ostream &o, QueueElement& qe){return qe.display(o);}
	
protected:
	ostream& display(ostream& os);
	unsigned _idx;
	double _tau;
	unsigned _queloc;
};

class QL_swap {
public:
	void exchange(QueueElement*& a, QueueElement*& b){
		QueueElement* temp = a;
		a = b;
		b = temp;
		unsigned ql = a->queueLoc();
		a->setQueueLoc(b->queueLoc()); 
		b->setQueueLoc(ql);
		
	}
};


class QL_Compare {
public:
	bool smaller(QueueElement* first, QueueElement* second){
		return first->tau() < second->tau();
	}
};

class IndexedQueue {
public:
	IndexedQueue():_pqueue(&_compare,&_swap),_sentinel(NULL){_indices.clear();}
	IndexedQueue(RanGen& ran, Model& pool, AllReactions& all, double t, double dt);
	~IndexedQueue();
	
	QueueElement* top() {return _pqueue.top();}
	QueueElement* pop() {return _pqueue.pop();}
	
	unsigned size() const {return _pqueue.size();}
	
	QueueElement* operator[](int location) {return _indices[location];}

	void update(int index) {_pqueue.update(index);}

        
        void init(RanGen& ran, AllReactions & all);
	
	void printQ() {_pqueue.print();} 
	
protected:
	QL_Compare _compare;
	QL_swap _swap;
        /** this is the queue which internally uses a heap.*/
	PriorityQueue<QueueElement*,QL_Compare, QL_swap> _pqueue;
        /** TODO why is this called indices?
         * it stores the queue elements... */
	vector<QueueElement*> _indices;
	QueueElement* _sentinel;
};

#endif

