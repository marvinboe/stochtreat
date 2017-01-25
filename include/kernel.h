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

#include "model.h"
#include "doctor.h"
#include "reactions.h"
#include "dependency.h"
#include "indexedqueue.h"
#include "rangen.h"
#include <fstream>
#include <sstream>

enum sim_type {DIAGNOSISRUN=0,TREATMENTRUN,RELAPSERUN};
class Kernel {
    public:
        /** Kernel contructor:
         * 1. Calls constructor of Model to initialize system in _pool.
         * 2. Calls constructor of DependencyGraph to initialize list of reactions (_allr) and the dependency graph _depend.
         * 3. Calls constructor of IndexedQueue with _pool and _allr as argument to initialize reaction queue.  */
        Kernel(RanGen& ran, Data& data, double time=0.0):_time(time),_data(data), _dt(data.dt()), _pool(data), 
        _depend(_pool, data, _allr, ran), _queue(ran, _pool, _allr, time, data.dt()),
        _endtime(data.getTmax_in_years()),_lsctime(-1.),_doctor(data.diagnosis_level(),data.reduction(),data.relapse_reduction()){};

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
         * int sim_type -> DIAGNOSIS, TREATMENT or RELAPSE.  */
        double execute(RanGen& ran, double t, int sim_type);

        /** Returns the reduction after (or while treatment).
         * reduction = 2.0 - log_10(burden)*/
        float getReduction() const;

        /**  get first timepoint without LSC in population.  */
        float get_nolsctime() const;

        /** Returns "true" if LSC is present in stem cell pool */
        bool hasLSC() const;

        /** sets the time (in years) until simulations  stop*/
        void set_ntime(double t){ _data.setTmax(t);}

        /** Makes all bound cells to cancer cells again*/
        void reset_treatment(RanGen& ran,double t);

        /** adds the sizes of stochastic compartments to the given data vector. */
        void addStochCompSizes(std::vector<double>& data) const;

        /** write all model data to the std::ostream */
        std::ostream& writeModel(std::ostream&);

        /** Reads all model data from the std::istream.
         * Requires same format as writeModel(ostream&).  */
        std::istream& readModel(std::istream& input);

        /* Reads model data from file given by path, runid and i. */
        bool read_model(std::string path, int runid, int i);

        /* Writes model data to file with filename given by path, runid and i. */
        bool write_model(std::string path, int runid, int i);

        /** Prints full doctors report to ostream. */
        void print_full_doctors_report(std::ostream& os) const{_doctor.print_patient_record(os);}

        /** Returns reference to doctor to give treatment information.*/
        const Doctor & doctor() const{ return _doctor; }

        /** Replaces a cancer cell with an immune cell in the lowest possible compartment */
        void introduce_immunity_inlowest();

        /** Introduces resistance in compartment k.*/
        void introduce_resistance(unsigned k);

    private:
        bool directMethod(RanGen& ran);
        bool nextMethod(RanGen& ran);
        void detUpdate();

        /** Reinitializes the kernel if cell counts or rates have changed. */
        void reinitialize(RanGen& ran,double simtime);

        bool stopsim(double time,int sim_type);

        double _time;
        double _stoptimer;
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
