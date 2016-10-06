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
#include <vector>

class Doctor{
    public:
        typedef double recorddata;
        Doctor();
        
        /** Reads the patient data at a single time point and adds to the
         * patients record in _data. */
        void take_bloodsample(double t, const Model& patient);

        /** Caclulates and returns the slope of burden decline based on the given patient record. */
        double calc_response(double from_time=-1., double end_time=-1., double timespan=-1.);

        /** Sets the initial reference value per patient to calculate the tumor
         * burden later. 1/initial_burden=_alpha = (NC + NB + (2.0 * NH)) / (NC + NB). */
        void calc_initial_reference(const Model& patient);

        /** Caculates and returns the tumor burden based on patient data.*/
        double calc_tumor_burden(const Model& patient);

        /** Caculates and returns the tumor burden at a specific already recorded timepoint (defaults to the last).*/
        double get_tumor_burden(double t=-1.);

        /** Takes a bloodsample every _sampling_timestep timesteps .*/
        void consult(double t, const Model & patient);

        /** Returns patient record to specified outstream. */
        void print_patient_record(std::ostream &os);


    private:
        std::vector<recorddata> _data;
        std::vector<double> _timepoints;

        double _next_timepoint;
        double _sampling_timestep;
        double _starttime;
        double _slope_timeintervall;

        double _alpha;

        bool _first_time_consulted;

        int find_timepoint(double t);

        double slope(const std::vector<recorddata>::const_iterator x_begin,const std::vector<recorddata>::const_iterator y_begin,unsigned int elements);

};
#endif
