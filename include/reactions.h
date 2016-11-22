/*
 *  reactions.h
 *  stochmut
 *
 *  Created by Tom Lenaerts on 2/1/11.
 *  Copyright 2011 Universite Libre de Bruxelles. All rights reserved.
 *
 */

#ifndef __REACTIONS_H
#define __REACTIONS_H

#include <string>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <cmath>
#include "model.h"
#include "rangen.h"


class Reaction { //X0 ->X0 + X1 or X0 -> X0 + X0
public:
	Reaction():_comp_reactant(0), _reactant(0),_comp_product1(0), _product1(0), _comp_product2(0), 
		_product2(0), _r(-1.0), _a(-1.0), _ptime(-1.0), _lasta(-1.0), _lastptime(-1.0), _time_zero(-1.0), _intype(-1), _ingraph(-1), _used(0){};

        /** Constructor for reaction. parameters:
         * cr   -> compartment of reactant
         * rct  -> reactant cell type 
         * cp1  -> compartment of product 1
         * p1   -> product 1 cell type
         * cp2  -> compartment of product2
         * p2   -> product 2 cell type
         * r    -> reaction rate
         * cell types: 
         * 0: healthy, 1: cancer, 2: immune, 3: treated.*/
	Reaction(int cr, int rct, int cp1, int p1, int cp2, int p2, double r):_comp_reactant(cr), 
		_reactant(rct),_comp_product1(cp1), _product1(p1), _comp_product2(cp2), _product2(p2), 
		_r(r), _a(-1.0), _ptime(-1.0), _lasta(-1.0), _lastptime(-1.0), _time_zero(-1.0), _intype(-1), _ingraph(-1), _used(0){};
	Reaction(const Reaction& other);
	virtual ~Reaction() {};
	
	unsigned int reactantComp() const {return _comp_reactant;}
	unsigned int reactant() const {return _reactant;}
	unsigned int product1Comp() const {return _comp_product1;}
	unsigned int product1() const {return _product1;}
	unsigned int product2Comp() const {return _comp_product2;}
	unsigned int product2() const {return _product2;}
	double rate() const {return _r;}
	
        /** Returns the index of the reaction in the dependency graph.*/
	unsigned int inGraph() const {return _ingraph;}

        /** Returns type of reaction:
         * MORAN=0,SELF_RENEWAL=1,DIFFERENTATION=2,TREATMENT=3. */
	unsigned int inType() const {return _intype;}

        /** Storing node type and the index of the node (in the graph)*/
	void setDG(unsigned type, unsigned loc) { _intype=type; _ingraph = loc;}
	
	double propensity() const {return _a;}
	double probability(double dt) const {return _a * dt;}

        /** Sets the propensity for this this reaction.
         * x -> number of cells for this reaction
         * needs to be set: _r -> reaction rate per cell */
	void setPropensity(double x) {_a = x *_r;}
	
        void setRate(double v) {_r=v;};

	double putativeTime() const {return _ptime;}

        /** calculates the putative time for this reaction,
         * adds to time "prev" and stores result in queue. */
	double calcPutativeTime(double ranval, double prev=0.0) { 
		_ptime = prev + (1.0 / _a ) * std::log(1.0/ranval);
		return _ptime;
	}
	double setPutativeTime(double v) { 
		_ptime = v;
		return _ptime;
	}
	
	void setLasta() { _lasta = _a;}
	double last_a() const {return _lasta;}

	void setLastPtime() { _lastptime= _ptime;}
	double last_ptime() const {return _lastptime;}

	void setTZero(double t) {_time_zero = t;} 
	double getTZero() const {return _time_zero;}
	
	unsigned used() const {return _used;}
	void incrUsed() {_used += 1;}
	void resetUsed() {_used = 0;}

	virtual Reaction& operator=(const Reaction& other);
	friend std::ostream & operator<<(std::ostream &o, Reaction& r){return r.display(o);}
	
        /** applies this reaction to the stem cell pool.
         * returns true if lsc vanished. */
	virtual bool apply(Model& pool,double time);
	virtual bool sufficientReactants(Model& pool);
        
        /** returns number of cells for given reaction type */
	virtual double reactantFactor(Model& pool);
	
	
protected:
	virtual std::ostream& display(std::ostream& os);
	unsigned int _comp_reactant;
	unsigned int _reactant;
	unsigned int _comp_product1;
	unsigned int _product1;
	unsigned int _comp_product2;
	unsigned int _product2;
	double _r; // rate constants of reaction (= compartment_rate * eps or compartmentr_rate * (1-eps))
	double _a; // _r * number of cells of type i
	double _ptime; // exponential function _a * exp (- (time*sum_a)) * dt
	double _lasta;
	double _lastptime;
	double _time_zero; // < timepoint when propensity reached zero
	unsigned _intype;
	unsigned _ingraph;
	unsigned _used;
};


class SelfRenewal : public Reaction {
public:
	SelfRenewal(int compartment, int ct, double rate):Reaction(compartment,ct,compartment,ct,compartment,ct,rate){};
	SelfRenewal(const SelfRenewal& other):Reaction(other){};
	virtual ~SelfRenewal() {};
	virtual SelfRenewal& operator=(const SelfRenewal& other);
	
};

class Differentation : public Reaction {
public:
	Differentation(int compartment, int ct, double rate):Reaction(compartment,ct,compartment+1,ct,compartment+1,ct,rate){};
	Differentation(const Differentation& other):Reaction(other){};
	virtual ~Differentation() {};
	virtual Differentation& operator=(const Differentation& other);
	
};


class Treatment : public Reaction {
public:
	Treatment(int compartment, double rate):Reaction(compartment,1,compartment,3,-1,-1,rate){};
	Treatment(const Treatment& other):Reaction(other){};
	virtual ~Treatment() {};
	virtual Treatment& operator=(const Treatment& other);
        virtual bool apply(Model & pool, double time);
};


class StemCellRenewal : public Reaction {
public:
	StemCellRenewal(RanGen* ran, int ct, double rate):Reaction(0,ct,0,ct,0,ct,rate),_ran(ran){};
	StemCellRenewal(const StemCellRenewal& other): Reaction(other),_ran(other.ran()){};
	virtual ~StemCellRenewal() {_ran = NULL;}
	virtual StemCellRenewal& operator=(const StemCellRenewal& other);
	RanGen* ran() const {return _ran;}
	virtual bool apply(Model& pool, double time);
	
private:
	RanGen* _ran;
};


class MoranReaction : public Reaction { //X0 + Y0 -> 2X0 + 2Y1
public:
	MoranReaction():Reaction(),_comp_reactant2(0), _reactant2(0) {};
	MoranReaction(int cr1, int rct1, int cr2, int rct2,int cp1, int p1, int cp2, int p2, double r):Reaction(cr1,rct1,cp1,p1,cp2,p2,r),_comp_reactant2(cr2), _reactant2(rct2){};
	MoranReaction(const MoranReaction& other);
	virtual ~MoranReaction () {};
	
	unsigned int reactant2Comp() const {return _comp_reactant2;}
	unsigned int reactant2() const {return _reactant2;}
	
	
	virtual MoranReaction& operator=(const MoranReaction& other);
	virtual bool apply(Model& pool, double time);
	virtual bool sufficientReactants(Model& pool);
//	virtual double reactantFactor(Model& pool);
	
protected:
	virtual std::ostream& display(std::ostream& os);
	unsigned int _comp_reactant2;
	unsigned int _reactant2;
};


class MoranRenewal : public MoranReaction {
public:
	MoranRenewal(int compartment, int ct, double rate):MoranReaction(compartment,ct,compartment,ct,compartment,ct,compartment+1,ct,rate){};
	MoranRenewal(const MoranRenewal& other):MoranReaction(other){};
	virtual ~MoranRenewal () {};
	virtual MoranRenewal& operator=(const MoranRenewal& other);
	
};

class MoranDifferentiation : public MoranReaction {
public:
	MoranDifferentiation(int compartment, int ct1, int ct2, double rate):MoranReaction(compartment,ct1,compartment,ct2,compartment,ct1,compartment+1,ct2,rate){};
	MoranDifferentiation(const MoranDifferentiation& other):MoranReaction(other){};
	virtual ~MoranDifferentiation () {};
	virtual MoranDifferentiation& operator=(const MoranDifferentiation& other);
	
};


class AllReactions  {
public:
	AllReactions():_sumprop(0.0){_all.clear();}
	~AllReactions(){
		while(_all.size() > 0){
			Reaction* r = _all[_all.size()-1];
			_all.pop_back();
			delete r;
		}
	}
	
	unsigned int size() const {return (unsigned)_all.size();}

        /** Adds reaction to the end of _all.
         * returns index of this element (=_all.size()-1) */
	unsigned int add(Reaction*);

        /** Returns pointer to reaction saved in _all[pos] */
	Reaction* operator[](unsigned pos);
	
	double propSum() const {return _sumprop;}
	void setPropSum(double v) { _sumprop = v;}

        void print(std::ostream &);
	
        std::vector<Reaction*>::iterator begin() {return _all.begin();}
        std::vector<Reaction*>::iterator end() {return _all.end();}
        std::vector<Reaction*>::const_iterator begin() const{return _all.begin();}
        std::vector<Reaction*>::const_iterator end() const {return _all.end();}
	
protected:
        std::vector<Reaction*> _all;
	double _sumprop;
};

#endif

