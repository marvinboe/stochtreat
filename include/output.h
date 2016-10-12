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
#include "doctor.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


struct Print_specifiers{

    Print_specifiers(std::string output_choice);
    bool per_patient=false;
    bool nolsctime=false;
    bool initialresponse=false;
    bool timetodiagnosis=false;
    bool fullburden=false;
    bool yearlyburden=false;
    bool overview_at_end=true;

    operator bool() const {return (per_patient||nolsctime||initialresponse||timetodiagnosis||yearlyburden);}

};

class Stats_Output{
    public:

        Stats_Output(std::string output_choice,unsigned no_stochcomps,bool treattest);

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

    clock_t _timer;
    int _output_specifier;
    int _patients;

    double _nolsc;
    double _diagnosed_nolsc;
    double _diagnosis_time;
    double _total_diagnosis_time;
    double _diagnosed;
    double _reachedreduction;
    double _total_timetoreduction;
    double _timetoreduction;
    double _burden_after_treatment;
    std::string _yearlyburden;
    std::vector<double> _avgsize;

    std::vector<double> _redresult;

    int _no_recurrence_patients;
    unsigned _recurrence_count;
    unsigned _nolsc_recurrence_count;

    bool _lsc_at_diagnosis;
    bool _nolsc_treattest;
    bool _diagnosis_reached;

    bool _treattest;

    Print_specifiers _print;

};

#endif
