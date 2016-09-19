/*
 *  kernel.h
 *  stochmut
 *
 *  Created by Tom Lenaerts on 2/5/11.
 *  Copyright 2011 Universite Libre de Bruxelles. All rights reserved.
 *
 */

#ifndef __KERNEL_H
#define __KERNEL_H

#include <fstream>
#include "model.h"
#include "reactions.h"
#include "dependency.h"
#include "indexedqueue.h"
#include "rangen.h"

class Kernel {
public:
    /** Kernel contructor:
     * 1. Calls constructor of Model to initialize system in _pool.
     * 2. Calls constructor of DependencyGraph to initialize list of reactions (_allr) and the dependency graph _depend.
     * 3. Calls constructor of IndexedQueue with _pool and _allr as argument to initialize reaction queue.
     */
	Kernel(RanGen& ran, Data& data, unsigned numstochs, double time=0.0):_time(time),_data(data), _dt(data.dt()), _pool(data,numstochs), 
	_depend(_pool, data, _allr, ran), _queue(ran, _pool, _allr, time, data.dt()),_endtime(data.ntimes()*(data.dt()/365.0)),_lsctime(-1){};


	Kernel(RanGen& ran, Data& data,std::istream& is, double time=0.0):_time(time),_data(data), _dt(data.dt()),
            _pool(data,is),_depend(_pool, data, _allr, ran), _queue(ran, _pool, _allr, time, data.dt()),_endtime(data.ntimes()*(data.dt()/365.0)),_lsctime(-1){};

        /** Reinitializes the kernel if cell counts have changed. */
        void reinitialize(Model& pool,RanGen& ran);

	void printAll();
	double execute(RanGen& ran, double t, bool treat);
	bool reachedDiagnosis();
	float getDiagnosis();
	bool reachedReduction();
	float getReduction();
	float whenReduction();

        /**  get first timepoint without LSC in population.  */
        float get_nolsctime();

        /** Returns "true" if LSC is present in stem cell pool */
	bool hasLSC();

	float burden();
	void addStochCompSizes(double* data);

        /** write all model data to the std::ostream */
	ostream& writeModel(ostream&);

        /** Reads all model data from the std::istream.
         * Requires same format as writeModel(ostream&).  */
        std::istream& readModel(std::istream& input);
	
private:
	void adjustReactions(RanGen& ran,unsigned compartment, unsigned type);
	bool directMethod(RanGen& ran);
	bool nextMethod(RanGen& ran);
	void detUpdate();
	void treatCells(RanGen& ran);
	
	
	double _time;
	Data _data;
	double _dt;
	Model _pool;
	AllReactions _allr;
	DependencyGraph _depend;
	IndexedQueue _queue;
	double _endtime;
	double _lsctime;
	
};

#endif
