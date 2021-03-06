/* Output relevent data
Copyright 2016 Marvin A. Böttcher

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#ifndef __OUTPUT_H
#define __OUTPUT_H

#include "kernel.h"
#include "data.h"
#include "doctor.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>


struct Print_specifiers{

    Print_specifiers(std::string output_choice);
    bool per_patient=false;
    bool nolsctime=false;
    bool initialresponse=false;
    bool timetodiagnosis=false;
    bool timetoreduction=false;
    bool fullburden=false;
    bool yearlyburden=false;
    bool overview_at_end=true;
    bool relapsetime=false;
    bool timepointburden_median=false;
    bool timepointburden_full=false;
    bool treat_dynamics=false;
    bool treat_dynamics_all=false;
    bool relapse_dynamics=false;
    bool relapse_dynamics_all=false;
    bool last_burden_output=false;
    bool print_simparams=false;

    /** Returns "true" if single patient data will be printed,
     * i.e. decides if a newline is printed after each patient.
     * !Update this if you add more output cases! */
    operator bool() const;

};

struct Timepoint_measurements{
    const int len=6;
    Timepoint_measurements():v(len,std::vector<double>(0)){}
    typedef std::vector<double> rundata;
    std::vector<double> t {0.25,0.5,0.75,1.0,1.25,1.5};
    std::vector<rundata> v;
    std::vector<double> return_av() const; 
    std::vector<double> return_std() const; 
};

class Stats_Output{
    public:

        Stats_Output(std::string output_choice,Simulation_Parameters simparams);

        /** Initialises per patient variables .*/
        void initialize_per_patient(int patient_id);

        /** Prints out information at the end of simulation.
         * This can include statistics over patients etc.*/
        void print_at_end() const;

        /**  Prints out information for each patient. Stores 
         * values relevant for later output.*/
        void print_save_patient(const Kernel & ker);

        /**  Prints out information for each patient. */
        void print_patient(const Kernel& ker) const;

        /**  Stores values relevant for later output.*/
        void save_patient(const Kernel & ker);

        /** stores relevant data after diagnosis is reached.*/
        void save_data_after_diagnosisrun(const Kernel & ker, double time);

        /** stores relevant data after treatment.*/
        void save_data_after_treatment(const Kernel & ker, double time);

        /** stores relevant data after relapse run.*/
        void save_data_after_relapse(const Kernel & ker, double time);

    private:
        double calc_median(std::vector<double> x) const;
        double lin_interpol(double v0, double v1, double t) const;
        double calc_quantile(const std::vector<double>& invec, double quant) const;
        bool print_dynamics( const std::vector<std::vector<double>>& burden2dvector) const;
        bool print_dynamics_all( const std::vector<std::vector<double>>& burden2dvector) const;
        bool print_relapse_dynamics_all() const;

        clock_t _timer;
        int _patients;

        double _nolsc;
        double _diagnosed_nolsc;
        double _diagnosis_time;
        double _total_diagnosis_time;
        double _diagnosed;
        double _reachedreduction;
        double _total_timetoreduction;
        double _timetoreduction;
        double _timetorelapse;
        double _timebeforerelapserun;
        double _burden_after_treatment;
        double _c_instoch_after_treatment;
        double _c_inneutral_after_treatment;
        double _resshare_treat;
        double _initialburden_alpha;
        double _resshare_relapse;
        std::string _yearlyburden;
        std::vector<double> _avgsize;
        std::vector<std::vector<double>> burden_record;
        std::vector<std::vector<double>> _relapse_burden_record;
        std::vector<std::vector<double>> _norelapse_burden_record;

        double _treat_dynamics_interval;

        std::vector<double> _redresult;
        Timepoint_measurements _timepoints_measure;
        int _no_recurrence_patients;
        unsigned _recurrence_count;
        unsigned _nolsc_recurrence_count;

        bool _lsc_at_diagnosis;
        bool _nolsc_treattest;
        bool _diagnosis_reached;


        Simulation_Parameters _simparams;

        Print_specifiers _print;

};

#endif
