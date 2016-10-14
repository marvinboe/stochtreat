/*
 *  Model.h
 *  stochmut
 *
 *  Created by Tom Lenaerts on 1/31/11.
 *  Copyright 2011 Universite Libre de Bruxelles. All rights reserved.
 *
 */

#ifndef __Model_H
#define __Model_H

#include <string>
#include <istream>
#include <cmath>
#include <iostream>
#include "data.h"
#include "rangen.h"


enum celltypes {H=0,C,I,B}; //H=healthy, C=cancer, I=resitant to TKI and B=Bound to TKI.


class Model {
public:
	Model(Data data, unsigned int numstoch);
	Model(Data data, std::istream & is);
	Model(const Model& other);
	~Model(){
		delete[] _compartments;
		delete[] _rates;
		delete[] _previous;
	}
	
	unsigned int numStoch() const {return _numstoch;}
	unsigned int numComp() const {return _numcomp;}
	
	void setRate(unsigned int k, double v);
	double getRate(unsigned k) const;

	void setTreatRate( double v){_treatrate=v;};
	double getTreatRate() const {return _treatrate;};
	
        /** returns numbers of cells in compartment k.*/
	double getN(unsigned int k) const;
	void setH(unsigned int k, double v); 
	double getH(unsigned int k) const; 
	void setC(unsigned int k, double v); 
	double getC(unsigned int k) const; 
	void setI(unsigned int k, double v); 
	double getI(unsigned int k) const; 
	void setB(unsigned int k, double v); 
	double getB(unsigned int k) const; 

	double retrieveN(unsigned int k) const;
	void storeH(unsigned int k, double v); 
	double retrieveH(unsigned int k) const; 
	void storeC(unsigned int k, double v); 
	double retrieveC(unsigned int k) const; 
	void storeI(unsigned int k, double v); 
	double retrieveI(unsigned int k) const; 
	void storeB(unsigned int k, double v); 
	double retrieveB(unsigned int k) const; 
	
	double lastN() const {return getN(_numcomp-1);}
	double lastH() const {return getH(_numcomp-1);}
	double lastC() const {return getC(_numcomp-1);}
	double lastI() const {return getI(_numcomp-1);}
	double lastB() const {return getB(_numcomp-1);}
	
        /** Returns the time of diagnosis.*/
	double diagRes() const {return _diagnosis;}

        /** Sets the time of diagnosis.*/
	void setDiagRes(double v) {_diagnosis = v;}

        /** Returns time of CML reduction. */
	double when() const {return _when;}
        /** Sets the time when reduction is reached. */
	void setWhenReduction(double v) {_when = v;}

	friend std::ostream & operator<<(std::ostream &o, Model& p){return p.display(o);}
	friend std::istream & operator>>(std::istream &i, Model& p){return p.read(i);}
	Model& operator=(const Model&);

        /** Returns the number of cells of type t in compartment k after last deterministic
         * time step. Returns values from _previous container.*/
	double retrieve(unsigned k, unsigned t);
        /** Sets the number of cells of type t in compartment k in _previous container.*/
	void store(unsigned k, unsigned t, double v);


        /** Makes all bound cells to cancer cells again*/
        void reset_treatment();

        /** Returns actual number of cells in pool of type t and 
         * in compartment k. */
	double get(unsigned k, unsigned t);

        /** Sets number of cells of type t in compartment k to v.*/
	void set(unsigned k, unsigned t, double v);
        /** Increments number of cells of type t in compartment k by v.*/
	void incr(unsigned k, unsigned t, double v);
        /** Decrements number of cells of type t in compartment k by v.*/
	void decr(unsigned k, unsigned t, double v);
	
        /** Stores all compartments cell counts.
         * Stores the cell counts of all cell types in each compartment
         * in the corresponding vector. !!!also sets the cell count 
         * of each cell type in the deterministic compartments to zero!!!.*/
	void memorize();

	bool updateDet(unsigned k, Data& data);

        /** returns True if diagnosis is reached.*/
	bool diagnosis(const Data& data) const;
	bool reduction(const Data& data) const;
	float getReduction() const;

        /** returns true if LSC is in stem cell pool */
	bool containsLSC() const;

	bool treatDeterministically(unsigned k, double amount);
	bool treatStochastically(unsigned k, double rate, RanGen& ran);

        /** Calculates alpha from cell numbers. Required to 
         * calculate reduction of disease burden later.*/
	void calcAlpha();
	double getAlpha() const {return _alpha;}
	
	double diseaseBurden() const;

        /** print cell numbers.
         * <HSC> <LSC> */
        void print_cells(std::ostream &,double _time);

        /** makes one cancer cell in compartment k immune. Returns true 
         * if successfull (at least one cell existed), otherwise false.*/
        bool manual_mutation(unsigned int k,unsigned int celltype_from=C,unsigned int celltype_to=I);
	
private:
	double myround(double val);
	double mylog(double p1, double base) const;

        std::ostream& display(std::ostream&) const;
        std::istream& read(std::istream& is);

	unsigned int _endstoch;
	double* _compartments;
	double* _previous;
	double* _rates;
	unsigned _numstoch;
	unsigned _numcomp;
	double _alpha;
        double _treatrate;

        /** time of diagnosis (in years?!)
         */
	double _diagnosis;

        /** time of reduction (in years?!)
         */
	double _when;
};

#endif
