/* process patient data for stochtreat.
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

#ifndef __DOCTOR_H
#define __DOCTOR_H

#include "model.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <vector>
#include <iostream>

class Doctor{
    public:
        typedef double recorddata;
        Doctor();
        Doctor(double diagnosis_level, double full_reduction, double relapse_reduction);
        
        /** Reads the patient data at a single time point and adds to the
         * patients record in _data. */
        void take_bloodsample(double t, const Model& patient);

        /** Caclulates and returns the slope of burden decline based on the given patient record. */
        double calc_response(double from_time=-1., double end_time=-1., double timespan=-1.) const;

        /** Sets the initial reference value per patient to calculate the tumor
         * burden later. 1/initial_burden=_alpha = (NC + NB + (2.0 * NH)) / (NC + NB). */
        void calc_initial_reference(double time,const Model& patient);

        /** Caculates and returns the tumor burden based on patient data.*/
        double calc_tumor_burden(const Model& patient) const;

        /** Caculates and returns the tumor burden at a specific already recorded timepoint (defaults to the last).*/
        double get_tumor_burden(double t=-1.) const;

        /** Returns number of cancer cells in the stochastic compartments
         * at a specific already recorded timepoint (defaults to the last).*/
        double get_c_instoch(double t=-1.) const;

        /** Returns number of cancer cells in the neutral compartments
         * at a specific already recorded timepoint (defaults to the last).*/
        double get_c_inneutral(double t=-1.) const;

        /** Returns the full tumor burden at specified
         * interval times.*/
        std::vector<double> get_burden_at_interval(double intdays) const;

        /** Takes a bloodsample every _sampling_timestep timesteps .*/
        void consult(double t, const Model & patient);

        /** Returns patient record to specified outstream. */
        void print_patient_record(std::ostream &os) const;

        /** Returns burden in yearly steps in a std::vector. */
        std::vector<double> get_yearly_burden() const;

        /** Calculates and returns the share of resistant cells in patient. */
        double calc_resistant_share(const Model& patient) const;

        /** Caculates and returns share of resistant cells at timepoint (defaults to the last).*/
        double get_resistant_share(double t=-1.) const;

        /** Calculates the log reduction of tumor burden. Returns -1 if burden=0. */
        double get_logreduction(double t=-1.) const;

        /** Returns true if required reduction level is reached.*/
        bool reduction_reached(double l=-1.,double t=-1.) const;

        /** Returns the timepoint when reduction was reached (for the first time).*/
        double reduction_time(double l=4.) const;

        /** Returns if cell count reaches diagnosis level.*/
        bool diagnosis_reached(double level=-1., double t=-1.) const;

        /** Returns if burden reaches relapse level.*/
        bool relapse_reached(double level=-1., double t=-1.) const;

        /** Returns the starting burden of cancer cells. */
        double return_initial_cratio() const {return _alpha;}
    private:
        std::vector<double> _burden_data;
        std::vector<double> _res_share_data;
        std::vector<double> _lastn_data;
        std::vector<double> _timepoints;
        std::vector<double> _bcrabl_instoch;
        std::vector<double> _bcrabl_inneutral;

        double _next_timepoint;
        double _sampling_timestep;
        double _starttime;
        double _starttime_treatment;
        double _slope_timeintervall;

        double _diagnosis_level;
        double _full_reduction;
        double _relapse_reduction;

        double _alpha;

        bool _first_time_consulted;

        /** returns index corresponding to 
         * timepoint t. TODO handling of missing values. */
        int find_timepoint(double t) const;

        double slope(const std::vector<recorddata>::const_iterator x_begin,const std::vector<recorddata>::const_iterator y_begin,unsigned int elements) const;

};
#endif
