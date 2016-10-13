
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
}

Stats_Output::Stats_Output(std::string output_choice,unsigned no_stochcomps,bool treattest):
    _print(output_choice),_treattest(treattest){


        _nolsc = 0;
        _diagnosed_nolsc = 0;
        _total_diagnosis_time=0;
        _diagnosed=0;
        _reachedreduction = 0;
        _total_timetoreduction = 0;
        _burden_after_treatment=-1.;
        _avgsize.resize(no_stochcomps+1);

        _redresult.clear();

        _no_recurrence_patients=0;
        _recurrence_count=0;
        _nolsc_recurrence_count=0;

        _patients=0;

        _timer=clock();

        std::cout <<"#output info: ";
        if (_print.nolsctime) std::cout <<"<nolsctime>";
        if (_print.timetodiagnosis) std::cout <<"<time_to_diag>";
        if (_print.timetoreduction) std::cout <<"<time to reduction>";
        if (_print.initialresponse){
            std::cout <<"<init. response>"<<"<lsc at diag>";
            if (!_print.yearlyburden) std::cout <<"<burden>";
            if (_treattest)
                std::cout <<"<relapse>";
        }
        if (_treattest && _print.relapsetime)
            std::cout <<"<timetorelapse>";

        if (_print.yearlyburden) std::cout <<"<yearlyburden>";
        std::cout << std::endl;
}

void Stats_Output::initialize_per_patient(int patient){
    _lsc_at_diagnosis=true;
    _nolsc_treattest=false;
    _diagnosis_reached=false;
    _burden_after_treatment=-1.;
    _timetorelapse=-1.;
    _yearlyburden.clear();

    _timetoreduction=0.;
    ++_patients;

    if (_print.per_patient){
        std::cout << "#patient " << patient<<" "<< std::endl;
    }
}


void Stats_Output::save_data_after_diagnosisrun(const Kernel& ker, double time){

    ker.addStochCompSizes(_avgsize);
    if(!ker.hasLSC())
        _nolsc +=1;

    if (ker.reachedDiagnosis()){
        _diagnosis_time=time;
        _diagnosed +=1;
        _diagnosis_reached=true;
        if (_treattest) _no_recurrence_patients++;

        _total_diagnosis_time += time;
        if(!ker.hasLSC()){
            _lsc_at_diagnosis=false;
            _diagnosed_nolsc +=1;
            if (_treattest){
                _nolsc_treattest=true;
            }
        }
    }
}

void Stats_Output::save_data_after_treatment(const Kernel &ker, double time){

    if(ker.reachedReduction()){
        _reachedreduction +=1;
        _timetoreduction=(ker.whenReduction() - _diagnosis_time);
        _total_timetoreduction +=_timetoreduction;
        _redresult.push_back(_timetoreduction);
    }
    _burden_after_treatment=ker.doctor().get_tumor_burden();
    std::vector<double> yearlyburden=ker.doctor().get_yearly_burden();
    std::stringstream strs;
    for (auto x : yearlyburden){
        strs << x<<" ";
    }
    _yearlyburden=strs.str();
    if (_yearlyburden.length()>0){
        _yearlyburden.pop_back();
    }

}

void Stats_Output::save_data_after_relapse(const Kernel &ker, double time){
    if(ker.reachedDiagnosis()) {
        _recurrence_count++;
        if (_nolsc_treattest) 
            _nolsc_recurrence_count++;
        _timetorelapse=time-(_timetoreduction+_diagnosis_time);
    }
}

void Stats_Output::print_patient(const Kernel& ker) const{
    if (_diagnosis_reached){

        if (_print.nolsctime)     std::cout  << ker.get_nolsctime() <<" ";
        if (_print.timetodiagnosis)
            std::cout << _diagnosis_time << "  " ;
        if (_print.timetoreduction)
            std::cout << _timetoreduction << "  ";

        if (_print.initialresponse){
            std::cout <<ker.doctor().calc_response()<<" ";
            std::cout <<_lsc_at_diagnosis<<" ";
            if (!_print.yearlyburden) std::cout <<_burden_after_treatment<<" ";
            if (_treattest)
                std::cout <<ker.reachedDiagnosis()<< " ";
        }
        if (_treattest && _print.relapsetime)
            std::cout <<_timetorelapse<<" ";

        if (_print.yearlyburden) 
            std::cout <<_yearlyburden<<" ";

        if (_print.fullburden){
            std::cout <<std::endl<<"#full doctor report"<<std::endl;
            ker.print_full_doctors_report(std::cout);
        }
        if (_print) std::cout <<std::endl;

    }


}

void Stats_Output::print_at_end() const{

    if (!_print.overview_at_end) return;
    if (_treattest){
        std::cout <<"#results cancer recurrence: <ratio> <recurrences> <total. diag.> <nolsc_ratio> <nolsc_recurrences> <no_lscdiags>"<<std::endl;
        if (_print) std::cout <<"# ";
        std::cout <<_recurrence_count/double(_no_recurrence_patients)
         <<" "<<_recurrence_count   <<" "  << _no_recurrence_patients<< " "<<_nolsc_recurrence_count/double(_diagnosed_nolsc)<<" "<<_nolsc_recurrence_count <<" "<< _diagnosed_nolsc<< std::endl;
    }

    std::cout << "#Real time elapsed in seconds: " << ((double)clock()-_timer)/CLOCKS_PER_SEC << std::endl;
    std::cout << "#Average time to diagnosis " << (_diagnosed > 0?(_total_diagnosis_time / (double) _diagnosed):0) << std::endl;
    std::cout << "#Fraction diagnosed "<< (_diagnosed /(double) _patients) << std::endl;
    std::cout << "#Fraction with no LSC " << (_nolsc / (double) _patients) << std::endl;
    std::cout << "#Fraction diagnosed with no LSC " << (_diagnosed > 0?(_diagnosed_nolsc / (double) _diagnosed):0) << std::endl;

    std::cout << "#<reduction freq.> <#of reductions> <diagnosed> <noscl at dignose>" << std::endl;
    if (_print || !_treattest) std::cout <<"# ";
    std::cout << ((_reachedreduction > 0&&_diagnosed>0)?(_reachedreduction / (double) _diagnosed):0)
        << " " << _reachedreduction << " " << _diagnosed << " " << _diagnosed_nolsc<< std::endl;


    double stddev = 0;
    double avg = (_reachedreduction > 0?(_total_timetoreduction / (double) _reachedreduction):0) ;
    for(unsigned int i=0; i < _redresult.size(); i++){
        stddev += pow((double)(_redresult[i] - avg), 2.0);
    }
    stddev = stddev / (double)_redresult.size();
    stddev = sqrt(stddev);
    std::cout << "#Average time to reduction "<< avg << "\t" << stddev << std::endl;
    for (int i=0; i< _avgsize.size(); i++) {
        std::cout << "#avg size comp " << i << " = " << (_avgsize[i]/(double)_patients) << std::endl;
    }

}
