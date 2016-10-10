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
#include "doctor.h"
#include "reactions.h"
#include "dependency.h"
#include "indexedqueue.h"
#include "rangen.h"

class Kernel {
    public:
        /** Kernel contructor:
         * 1. Calls constructor of Model to initialize system in _pool.
         * 2. Calls constructor of DependencyGraph to initialize list of reactions (_allr) and the dependency graph _depend.
         * 3. Calls constructor of IndexedQueue with _pool and _allr as argument to initialize reaction queue.  */
        Kernel(RanGen& ran, Data& data, unsigned numstochs, double time=0.0):_time(time),_data(data), _dt(data.dt()), _pool(data,numstochs), 
        _depend(_pool, data, _allr, ran), _queue(ran, _pool, _allr, time, data.dt()),_endtime(data.getTmax_in_years()),_lsctime(-1),_doctor(){};

        /** Kernel contructor that reads model data from std::istream& is instead of initialising new. Steps:
         * 1. Calls constructor of Model to initialize system in _pool.
         * 2. Calls constructor of DependencyGraph to initialize list of reactions (_allr) and the dependency graph _depend.
         * 3. Calls constructor of IndexedQueue with _pool and _allr as argument to initialize reaction queue.  */
        Kernel(RanGen& ran, Data& data,std::istream& is, double time=0.0):_time(time),_data(data), _dt(data.dt()),
        _pool(data,is),_depend(_pool, data, _allr, ran), _queue(ran, _pool, _allr, time, data.dt()),_endtime(data.getTmax_in_years()),_lsctime(-1){};


        /** Prints cell pool and reactions data to std::cout.*/
        void printAll();

        /** Executes the simulation starting with the current state of the
         * cell pool. Returns the time after simulation stopped in years.
         * parameters:
         * RanGen& ran -> random number generator for stochastic updates
         * double t -> start time of the simulation (This has to correspond 
         *             to the cell pool time, otherwise stochastic updates
         *             won't work!)
         * bool treat -> wether or not to apply treatment in this run.  */
        double execute(RanGen& ran, double t, bool treat);

        /** Returns true if cell count in cell pool reached diagnosis level. */
        bool reachedDiagnosis();
        /** Returns the time of (first) diagnosis for the cell pool. */
        float getDiagnosisTime();
        /** Returns true if reduction is reached in cell pool. Depends on the number of cells.*/
        bool reachedReduction();
        /** Returns the reduction after (or while treatment).
         * reduction = 2.0 - log_10(burden)*/
        float getReduction();

        /** Returns time when required reduction level is reached. */
        float whenReduction();

        /**  get first timepoint without LSC in population.  */
        float get_nolsctime();

        /** Returns "true" if LSC is present in stem cell pool */
        bool hasLSC();

        /** sets the time (in years) until simulations  stop*/
        void set_ntime(double t){ _data.setTmax(t);}

        /** Makes all bound cells to cancer cells again*/
        void reset_treatment(RanGen& ran,double t);

        /** Calculates and returns the disease burden for this particular
         * patient. Based on "alpha" which needs to be calculated before,
         * typically when diagnosis is reached. */
        float burden();

        void addStochCompSizes(double* data);

        /** write all model data to the std::ostream */
        ostream& writeModel(ostream&);

        /** Reads all model data from the std::istream.
         * Requires same format as writeModel(ostream&).  */
        std::istream& readModel(std::istream& input);

        /** Returns the treatment response at start of treatment. */
        double initial_treatment_response(){ return _doctor.calc_response(); }

        /** Prints full doctors report to ostream. */
        void print_full_doctors_report(std::ostream& os) {_doctor.print_patient_record(os);}


        const Doctor & doctor() { return _doctor; }
    private:
        bool directMethod(RanGen& ran);
        bool nextMethod(RanGen& ran);
        void detUpdate();

        /** Reinitializes the kernel if cell counts or rates have changed. */
        void reinitialize(Model& pool,RanGen& ran,double simtime);

        double _time;
        Data _data;
        double _dt;
        Model _pool;
        AllReactions _allr;
        DependencyGraph _depend;
        IndexedQueue _queue;
        double _endtime;
        double _lsctime;

        Doctor _doctor;

};

#endif
