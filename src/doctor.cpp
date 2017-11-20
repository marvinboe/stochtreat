
#include "doctor.h"

Doctor::Doctor(){
    //do nothing
    _next_timepoint=0.;
    _sampling_timestep=1.;
    _slope_timeintervall=62.;
    _diagnosis_level=1.e12;
    _starttime=0.;
    _starttime_treatment=-1.;
    _first_time_consulted=true;
    _alpha=0.;
}

Doctor::Doctor(double diagnosis_level, double full_reduction, double relapse_reduction):
    _diagnosis_level(diagnosis_level),_full_reduction(full_reduction),_relapse_reduction(relapse_reduction){
    //do nothing
    _next_timepoint=0.;
    _sampling_timestep=1.;
    _slope_timeintervall=62.;
    _starttime=0.;
    _starttime_treatment=-1.;
    _first_time_consulted=true;
    _alpha=0.;
}


double Doctor::get_tumor_burden(double t) const{
    if (t<0. && _timepoints.size()>0)
        return _burden_data.back();

    int i=find_timepoint(t);
    if (i>=0 && _timepoints.size()>0)
        return _burden_data[i];
    return -1;
    // do nothing 
}

double Doctor::get_c_instoch(double t) const{
    if (t<0. && _timepoints.size()>0)
        return _bcrabl_instoch.back();

    int i=find_timepoint(t);
    if (i>=0 && _timepoints.size()>0)
        return _bcrabl_instoch[i];
    return -1;
    // do nothing 
}


void Doctor::calc_initial_reference(double time,const Model& patient){
	double NC=patient.lastC()+patient.lastI();
	double NB=patient.lastB();
	double NH=patient.lastH();
	_alpha = (NC + NB + (2.0 * NH)) / (NC + NB);
        _starttime_treatment=time;
        int i = find_timepoint(time);
        _burden_data[i]=100.;
        // std::cout <<"init alpha debug: "<<_alpha<<" "<<NC<<" "<<NB<<" "<<NH<<std::endl;
}

double Doctor::calc_tumor_burden(const Model& patient) const{
	double NC=patient.lastC()+patient.lastI();
	double NB=patient.lastB();
	double NH=patient.lastH();

	double burden = 0.;
        if (_alpha > 0.)
	    burden = _alpha*100.0 * ((NC + NB) / (NC + NB + (2.0 * NH)));//_alpha* ((NC + NB) / (NC + NB + (2.0 * NH)));//
        // std::cout << burden << std::endl;
	return burden;
}

double Doctor::calc_resistant_share(const Model& patient) const{
	double NC=patient.lastC();
        double NI=patient.lastI();
	double NB=patient.lastB();
	double NH=patient.lastH();

	return NI/(NC+NB+NH+NI);
}

void Doctor::take_bloodsample(double t, const Model & patient){

    double burden = calc_tumor_burden(patient);
    double share = calc_resistant_share(patient);
    _burden_data.push_back(burden);
    _res_share_data.push_back(share);
    _lastn_data.push_back(patient.lastN());
    _timepoints.push_back(t);
    _bcrabl_instoch.push_back(patient.get_c_instoch());
}


int Doctor::find_timepoint(double t) const{
    if (t<0.) return -1;
    int i =_timepoints.size()-1;
    if (i>0 && t > (2*_timepoints.back()-_timepoints.rbegin()[1])) return -1; //if t > t_max+t_step 
    while (i > 0 && _timepoints[i] > t) --i; 
    return i;
}

double Doctor::calc_response(double from_time, double end_time, double timespan) const{

    if (from_time<0. && end_time <0. && timespan <0.){


        auto x_begin=_timepoints.begin();
        auto y_begin=_burden_data.begin();
        
        int i_start = find_timepoint(_starttime_treatment);
        int reglength = find_timepoint(_starttime_treatment+_slope_timeintervall)-i_start;
        std::advance(x_begin,i_start);
        std::advance(y_begin,i_start);

        if (reglength<0||i_start<0){
            return std::numeric_limits<double>::infinity();
            //this is bad TODO
        }

        return slope(x_begin,y_begin,reglength);

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


void Doctor::print_patient_record(std::ostream &os) const{
    os <<"# patient data: <time> <burden>"<<std::endl;
    for (unsigned int i=0; i< _timepoints.size(); ++i){
        os <<_timepoints[i]/365.<<" "<<_burden_data[i]<<" "<<_bcrabl_instoch[i]<<std::endl;
    }
}

std::vector<double> Doctor::get_yearly_burden() const{
    std::vector<double> returnvec;
    double t=_starttime_treatment+365.; //TODO something is one day off here???
    while (t<=_timepoints.back()){
        returnvec.push_back(get_tumor_burden(t));
        t+=365.;
    }
    // returnvec.push_back(get_tumor_burden(-1.));
    // std::cout <<t<<" "<<_timepoints.back()<<std::endl;
    return returnvec;
}

std::vector<double> Doctor::get_burden_at_interval(double intdays) const{
    std::vector<double> returnvec;
    double t=_starttime_treatment+0.; //TODO something is one day off here???
    while (t<=_timepoints.back()){
        returnvec.push_back(get_tumor_burden(t));
        t+=intdays;
    }
    // std::cout <<t<<" "<<_timepoints.back()<<std::endl;
    return returnvec;
}

double Doctor::get_resistant_share(double t) const{
    if (t<0.) t=(_timepoints.size()>0?_timepoints.back():0.);

    int i=find_timepoint(t);
    if (i>=0)
        return _res_share_data[i];
    else
        return -1.;
    // do nothing
}

double Doctor::get_logreduction(double t) const{ 
    double burden=get_tumor_burden(t);
    if (burden >0.)
        return -std::log10(burden/100.); //!burden starts at 100
    else 
        return -1.;
}

double Doctor::reduction_time(double l) const{

    for (unsigned int ti=0 ; ti < _timepoints.size(); ++ti){
        if (get_logreduction(_timepoints[ti]) >= l) return _timepoints[ti]/365.;
    }
    return -1.;
}

bool Doctor::reduction_reached(double l, double t) const {
    if (l<0.) l=_full_reduction;
    if (t<0.) t=(_timepoints.size()>0?_timepoints.back():-1.);
    if (t<0.) return false; 

    if (l<0.) return false; //for negative l: reduction never reached
    // std::cout <<"reduction_reached debug: "<<get_logreduction(t)<<" "<<l<<std::endl;
    return (get_logreduction(t)>=l);
}


bool Doctor::diagnosis_reached( double level, double t) const{
    if (level < 0.) level=_diagnosis_level;

    if (t<0.) t=(_timepoints.size()>0?_timepoints.back():-1.);
    if (t<0.) return false; 
    int i=find_timepoint(t);

    // std::cout <<"debug diagnosis: "<<_lastn_data[i]<<" "<<level<<" "<<i<<" "<<t<<std::endl;
    if ((i>=0) && (_lastn_data[i] >= std::pow(10.,level))){
        return true;
    }
    else return false;
}

bool Doctor::relapse_reached(double l, double t) const {
    if (l<0.) l=_relapse_reduction;
    if (t<0.) t=(_timepoints.size()>0?_timepoints.back():-1.);
    if (t<0.) return false; 

    // std::cout <<"relapse_reached debug: "<<l<<" "<<t<<" "<<get_logreduction(t)<<" "<<get_tumor_burden(t)<<std::endl;
    return (get_logreduction(t)<l && get_logreduction(t) >=0.);
}
