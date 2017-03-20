
#include "output.h"

Print_specifiers::Print_specifiers(std::string output_choice){

    if (output_choice.find("patient")!=std::string::npos){
        per_patient=true;
    }
    if (output_choice.find("nolsctime")!=std::string::npos){
        nolsctime=true;
    }
    if (output_choice.find("diagtime")!=std::string::npos){
        timetodiagnosis=true;
    }
    if (output_choice.find("reductiontime")!=std::string::npos){
        timetoreduction=true;
    }
    if (output_choice.find("initresponse")!=std::string::npos){
        initialresponse=true;
    }
    if (output_choice.find("fullburden")!=std::string::npos){
        fullburden=true;
    }
    if (output_choice.find("yearlyburden")!=std::string::npos){
        yearlyburden=true;
    }
    if (output_choice.find("relapsetime")!=std::string::npos){
        relapsetime=true;
    }
    if (output_choice.find("nooverview")!=std::string::npos){
        overview_at_end=false;
    }
    if (output_choice.find("3timepointsmedian")!=std::string::npos){
        three_timepoint_median=true;
    }
    if (output_choice.find("3timepointsfull")!=std::string::npos){
        three_timepoint_full=true;
    }
    if (output_choice.find("treatdynamics")!=std::string::npos){
        treat_dynamics=true;
    }
}
Print_specifiers::operator bool() const {
    return (per_patient||nolsctime||initialresponse||timetodiagnosis||
            timetoreduction||yearlyburden||three_timepoint_full||
            relapsetime);
}

Stats_Output::Stats_Output(std::string output_choice,unsigned no_stochcomps,Run_modes run_mode):
    _run_mode(run_mode),_print(output_choice){


        _nolsc = 0;
        _diagnosed_nolsc = 0;
        _total_diagnosis_time=0;
        _diagnosed=0;
        _reachedreduction = 0;
        _total_timetoreduction = 0;
        _burden_after_treatment=-1.;
        _avgsize.resize(no_stochcomps+1);

        _treat_dynamics_interval=0.05;

        _redresult.clear();

        _no_recurrence_patients=0;
        _recurrence_count=0;
        _nolsc_recurrence_count=0;

        _patients=0;

        _timer=clock();

        std::cout <<"#output info: ";
        if (_print.nolsctime) std::cout <<"<nolsctime> ";
        if (_print.timetodiagnosis) std::cout <<"<time_to_diag> ";
        if (_print.timetoreduction) std::cout <<"<time to reduction> ";
        if (_print.initialresponse){
            std::cout <<"<init. response>";//<<"<lsc at diag>"<<"<initial c_ratio>";
            // std::cout <<"<burden>";
            // if (_run_mode.treattest)
                // std::cout <<"<relapse>";
            std::cout <<" ";
        }
        if (_run_mode.treattest && _print.relapsetime)
            std::cout <<"<timetorelapse> ";

        if (_print.three_timepoint_full){
            for (int i=0; i<3; ++i){
                std::cout<<"<burden"<<_three_timepoints_measure.t[i]<<"> ";
            }
            std::cout <<" ";
        }

        if (_run_mode.resistance>=0){
            std::cout <<"<resistance_share_treat>";
            if (_run_mode.treattest) std::cout <<"<resistance_share_relapse>";
            std::cout <<" ";
        }

        if (_print.yearlyburden) std::cout <<"<yearlyburden>";
        std::cout << std::endl;
    }

void Stats_Output::initialize_per_patient(int patient){
    _lsc_at_diagnosis=true;
    _nolsc_treattest=false;
    _diagnosis_reached=false;
    _diagnosis_time=-1;
    _burden_after_treatment=-1.;
    _timetorelapse=-1.;
    _yearlyburden.clear();

    _timetoreduction=0.;
    _timebeforerelapserun=0.;
    ++_patients;

    if (_print.per_patient){
        std::cout << "#patient " << patient<<" "<< std::endl;
    }
}


void Stats_Output::save_data_after_diagnosisrun(const Kernel& ker, double time){

    ker.addStochCompSizes(_avgsize);
    if(!ker.hasLSC())
        _nolsc +=1;

    if (ker.doctor().diagnosis_reached()){
        _diagnosis_time=time;
        _diagnosed +=1;
        _diagnosis_reached=true;

        _total_diagnosis_time += time;
        if(!ker.hasLSC()){
            _lsc_at_diagnosis=false;
            _diagnosed_nolsc +=1;
            if (_run_mode.treattest){
                _nolsc_treattest=true;
            }
        }
    }
}

void Stats_Output::save_data_after_treatment(const Kernel &ker, double time){

    _initialburden_alpha=ker.doctor().return_initial_cratio();
    if(ker.doctor().reduction_reached()){
        if (_run_mode.treattest) _no_recurrence_patients++;
        _reachedreduction +=1;
        _timetoreduction=(ker.doctor().reduction_time() - _diagnosis_time);
        _total_timetoreduction +=_timetoreduction;
        _redresult.push_back(_timetoreduction);
    }
    for (unsigned i =0; i<3; ++i){
        _three_timepoints_measure.v[i].push_back(ker.doctor().get_tumor_burden((_diagnosis_time+_three_timepoints_measure.t[i])*365.));
    }
    _timebeforerelapserun=time;
    _burden_after_treatment=ker.doctor().get_tumor_burden();
    _resshare_treat=ker.doctor().get_resistant_share();

    std::vector<double> yearlyburden=ker.doctor().get_yearly_burden();
    std::stringstream strs;
    for (auto x : yearlyburden){
        strs << x<<" ";
    }
    _yearlyburden=strs.str();
    if (_yearlyburden.length()>0){
        _yearlyburden.pop_back();
    }
    if (_print.treat_dynamics){
        burden_record.push_back(ker.doctor().get_burden_at_interval(_treat_dynamics_interval*365));
    }

}

void Stats_Output::save_data_after_relapse(const Kernel &ker, double time){
    _resshare_relapse=ker.doctor().get_resistant_share();
    if(ker.doctor().relapse_reached()) {
        _recurrence_count++;
        _timetorelapse=time-(_timebeforerelapserun);
        if (_nolsc_treattest) 
            _nolsc_recurrence_count++;
    }
    else
        _timetorelapse=-2.;
}

void Stats_Output::print_patient(const Kernel& ker) const{

    if (_print.nolsctime)     std::cout  << ker.get_nolsctime() <<" ";
    if (_print.timetodiagnosis)
        std::cout << _diagnosis_time << "  " ;
    if (_print.timetoreduction)
        std::cout << _timetoreduction << "  ";

    if (_print.initialresponse){
        std::cout <<ker.doctor().calc_response()<<" ";
        // std::cout <<_lsc_at_diagnosis<<" ";
        // std::cout <<_initialburden_alpha<<" ";
        // std::cout <<_burden_after_treatment<<" ";
        // if (_run_mode.treattest)
        //     std::cout <<ker.doctor().diagnosis_reached()<< " ";
        std::cout <<" ";
    }
    if (_run_mode.treattest && _print.relapsetime)
        std::cout <<_timetorelapse<<"  ";

    if (_print.three_timepoint_full){
        if (_diagnosis_reached){
            for (int i=0; i<3; ++i){
                std::cout<<_three_timepoints_measure.v[i].back()<<" ";
            }
        }
        else {
            std::cout <<"-1 -1 -1 ";
        }
        std::cout <<" ";

    }

    if (_print.yearlyburden) 
        std::cout <<_yearlyburden<<"  ";

    if (_run_mode.resistance>=0){
        std::cout <<_resshare_treat<<"  ";
        if (_run_mode.treattest) std::cout <<_resshare_relapse<<" ";
    }

    if (_print.fullburden){
        std::cout <<std::endl<<"#full doctor report"<<std::endl;
        ker.print_full_doctors_report(std::cout);
    }
    if (_print||_run_mode.resistance>=0) std::cout <<std::endl;



}

void Stats_Output::print_at_end() const{

    if (_print.three_timepoint_median){
        std::cout <<"#median burden at three timepoint: ";
        for (int i=0; i<3; ++i) std::cout <<"<"<<_three_timepoints_measure.t[i]<<">";
        std::cout<<std::endl;
        for (int i=0; i<3; ++i){
            std::cout<<_three_timepoints_measure.return_median()[i]<<" ";
        }
        for (int i=0; i<3; ++i){
            std::cout<<_three_timepoints_measure.return_std()[i]<<" ";
        }
        std::cout<<std::endl;
    }

    if (_print.treat_dynamics){
        unsigned int i=0; //row number
        bool stop=false;
        while (!stop){
            unsigned no_columns=0;
            double sum=0.;
            std::vector<double> row;
            for (unsigned j=0; j<burden_record.size(); ++j){
                if (burden_record[j].size()>i){
                    row.push_back(burden_record[j][i]);
                    sum+=burden_record[j][i];
                    ++no_columns;
                }
            }
            if (no_columns==0){
                stop=true;
                break;
            }
            //median
            std::nth_element(row.begin(), row.begin() + row.size() / 2, row.end());
            double med= *std::next(row.begin(), row.size() / 2);

            std::cout <<i*_treat_dynamics_interval<<" "<<sum/double(burden_record.size())<<" "<<med<<std::endl;
            ++i;
        }
    }

    if (!_print.overview_at_end) return;

    if (_run_mode.treattest){
        std::cout <<"#results cancer recurrence: <ratio> <recurrences> <total. diag.> <nolsc_ratio> <nolsc_recurrences> <no_lscdiags>"<<std::endl;
        if (_print) std::cout <<"# ";
        std::cout <<_recurrence_count/double(_no_recurrence_patients)
            <<" "<<_recurrence_count   <<" "  << _no_recurrence_patients<< " "<<_nolsc_recurrence_count/double(_diagnosed_nolsc)<<" "<<_nolsc_recurrence_count <<" "<< _diagnosed_nolsc<< std::endl;
    }

    std::cout << "#Real time elapsed in seconds: " << ((double)clock()-_timer)/CLOCKS_PER_SEC << std::endl;
    std::cout <<"#<av. diagtime>"<<"<diagnosed frac>"<<"<nolsc frac>"<<"<nolsc diagnosed frac>"<<std::endl;
    std::cout <<"# "<< (_diagnosed > 0?(_total_diagnosis_time / (double) _diagnosed):0)<<" "
        << (_diagnosed /(double) _patients) <<" " << (_nolsc / (double) _patients) <<" "
        << (_diagnosed > 0?(_diagnosed_nolsc / (double) _diagnosed):0) << std::endl;

    std::cout << "#<reduction freq.> <#of reductions> <diagnosed> <noscl at dignose>" << std::endl;
    if (_print || !_run_mode.treattest) std::cout <<"# ";
    std::cout << ((_reachedreduction > 0&&_diagnosed>0)?(_reachedreduction / (double) _diagnosed):0)
        << " " << _reachedreduction << " " << _diagnosed << " " << _diagnosed_nolsc<< std::endl;


    double stddev = 0;
    double avg = (_reachedreduction > 0?(_total_timetoreduction / (double) _reachedreduction):-1.) ;
    for(unsigned int i=0; i < _redresult.size(); i++){
        stddev += pow((double)(_redresult[i] - avg), 2.0);
    }
    stddev = stddev / (double)_redresult.size();
    stddev = sqrt(stddev);
    std::cout << "#time to reduction avg="<< avg << " stddev=" << stddev << std::endl;
    std::cout << "#avg size comps. ";
    for (unsigned int i=0; i< _avgsize.size(); i++) {
        std::cout << "<" << i << "> ";
    }
    std::cout<<std::endl<<"# ";
    for (unsigned int i=0; i< _avgsize.size(); i++) {
        std::cout << (_avgsize[i]/(double)_patients) <<" ";
    }
    std::cout << std::endl;
}



std::vector<double> Three_timepoint_measurements::return_av() const{ 
    std::vector<double> av {0.,0.,0.};
    double number=v.back().size();
    for (int i=0; i<3; ++i){
        av[i]=std::accumulate(v[i].begin(),v[i].end(),0.)/number;
    }
    return av;
}

std::vector<double> Three_timepoint_measurements::return_std() const{ 
    std::vector<double> stdev {0.,0.,0.};
    std::vector<double> allmean=return_av();
    for (int i=0; i<3; ++i){
        double mean=allmean[i];
        double sq_sum = std::inner_product(v[i].begin(), v[i].end(), v[i].begin(), 0.0);
        stdev[i] = std::sqrt(sq_sum / v[i].size() - mean * mean);
    }
    return stdev;
}

std::vector<double> Three_timepoint_measurements::return_median() const {
    std::vector<double> returnvalues(3,std::numeric_limits<double>::quiet_NaN());
    for (int i = 0 ; i<3 ;++i){
        std::vector<double>x=v[i];

        if (x.size() == 0) return returnvalues;

        size_t n = 0.;
        if(x.size()%2 == 1)
            n=(x.size()-1)/2;
        else 
            n=x.size()/2;

        std::nth_element(x.begin(), x.begin()+n, x.end());

        double xn = x[n];
        if(x.size()%2 == 1) {
            returnvalues[i]= xn;
        }else {
            // std::nth_element(x.begin(), x.begin()+n-1, x.end());
            returnvalues[i]= 0.5*(xn+x[n-1]);
        }

    }
    return returnvalues;
}
