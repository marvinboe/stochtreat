
#include "doctor.h"

Doctor::Doctor(){
    //do nothing
    _next_timepoint=0.;
    _sampling_timestep=1.;
    _slope_timeintervall=62.;
    _starttime=0.;
    _first_time_consulted=true;
}


double Doctor::get_tumor_burden(double t) const{
    if (t<0.) t=(_timepoints.size()>0?_timepoints.back():0.);

    unsigned int i=find_timepoint(t);
    return _data[i];
    // do nothing
}


void Doctor::calc_initial_reference(const Model& patient){
	double NC=patient.lastC();
	double NB=patient.lastB();
	double NH=patient.lastH();
	_alpha = (NC + NB + (2.0 * NH)) / (NC + NB);
}

double Doctor::calc_tumor_burden(const Model& patient) const{
	double NC=patient.lastC();
	double NB=patient.lastB();
	double NH=patient.lastH();

	double burden = _alpha* ((NC + NB) / (NC + NB + (2.0 * NH)));//(_alpha*(100.0 * ((NC + NB) / (NC + NB + (2.0 * NH)))));
//	cout << burden << std::endl;
	return burden;
}

void Doctor::take_bloodsample(double t, const Model & patient){
    if (_data.size()==0){
        calc_initial_reference(patient);
    }
    double burden = calc_tumor_burden(patient);
    _data.push_back(burden);
    _timepoints.push_back(t);
}


int Doctor::find_timepoint(double t) const{
    int i =_timepoints.size()-1;
    while (i > 0 && _timepoints[i] > t) --i; 
    return i;
}

double Doctor::calc_response(double from_time, double end_time, double timespan) const{

    if (from_time<0. && end_time <0. && timespan <0.){


        const auto x_begin=_timepoints.begin();
        const auto y_begin=_data.begin();
        
        int i = find_timepoint(_starttime+_slope_timeintervall);

        if (i<0){
            return 100.;
            //this is bad TODO
        }

        return slope(x_begin,y_begin,i);

    }
    else {
        //to be implemented TODO
        return from_time;
    }
}

double Doctor::slope(const std::vector<recorddata>::const_iterator x_begin,const std::vector<recorddata>::const_iterator y_begin,unsigned int elements) const {
    const auto n    = elements;
    const auto s_x  = std::accumulate(x_begin,x_begin+elements, 0.0);
    const auto s_y  = std::accumulate(y_begin,y_begin+elements, 0.0);
    const auto s_xx = std::inner_product(x_begin,x_begin+elements, x_begin, 0.0);
    const auto s_xy = std::inner_product(x_begin,x_begin+elements, y_begin, 0.0);
    // std::cout <<"#slope debug: sx="<<s_x<<" sy="<<s_y<<" s_xx="<<s_xx<<" s_xy"<<s_xy<<std::endl;
    const auto a    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
    return a;
}

void Doctor::consult(double t, const Model& patient){
    if (_first_time_consulted){
        _starttime=t;
        _next_timepoint=t;
        _first_time_consulted=false;
    }

    while (_next_timepoint <= t){
        take_bloodsample(_next_timepoint,patient);
        _next_timepoint+=_sampling_timestep;
    }
}


void Doctor::print_patient_record(std::ostream &os){
    os <<"# patient data: <time> <burden>"<<std::endl;
    for (unsigned int i=0; i< _timepoints.size(); ++i){
        os <<_timepoints[i]<<" "<<_data[i]<<std::endl;
    }
}

