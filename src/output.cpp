
#include "output.h"

Stats_Output::Stats_Output(int output_specifier,unsigned no_stochcomps,bool treattest):
    _output_specifier(output_specifier),_treattest(treattest){

        _nolsc = 0;
        _diagnosed_nolsc = 0;
        _total_diagnosis_time=0;
        _diagnosed=0;
        _reachedreduction = 0;
        _total_timetoreduction = 0;
        _avgsize.resize(no_stochcomps+1);

        _redresult.clear();

        _no_recurrence_patients=0;
        _recurrence_count=0;
        _nolsc_recurrence_count=0;

        _patients=0;

        _timer=clock();
}

void Stats_Output::initialize_per_patient(int patient){
    _lsc_at_diagnosis=true;
    _nolsc_treattest=false;
    _diagnosis_reached=false;

    _timetoreduction=0.;
    ++_patients;

        if (_output_specifier==0){
            cout << "#patient " << patient << endl;
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
                if (_output_specifier==3){
                    cout << "#<years to diag.> <years to red.> <total> <nolsctime> "<<std::endl
                        << _diagnosis_time<< "  " 
                        << _timetoreduction << "  "
                        << ker.whenReduction() << " "
                        << ker.get_nolsctime() << endl;
                }

            }

}

void Stats_Output::save_data_after_relapse(const Kernel &ker, double time){
    if(ker.reachedDiagnosis()) {
        _recurrence_count++;
        if (_nolsc_treattest) 
            _nolsc_recurrence_count++;
    }
}

void Stats_Output::print_patient(const Kernel& ker) const{
    if (_diagnosis_reached){

            if (_output_specifier==4){
                std::cout <<ker.doctor().calc_response();
                std::cout <<" "<<_lsc_at_diagnosis;
                if (!_treattest){
                    std::cout <<" "<<ker.doctor().get_tumor_burden();
                }
                if (_treattest)
                    std::cout << " "<<ker.reachedDiagnosis();
                std::cout <<std::endl;
            }
            if (_output_specifier==5){
                ker.print_full_doctors_report(std::cout);
            }
            //			cout << "Reduction is " << ker.getReduction() << endl;
            //
        }//######### everything for case of diagnosis done
        

        if (_output_specifier==1){
            cout  << ker.get_nolsctime() << endl;
        }
        else if (_output_specifier==2){
            std::cout << _diagnosis_time << "  " 
                << _timetoreduction << "  "
                << ker.whenReduction() << " "
                << ker.get_nolsctime() << endl;
        }

}

void Stats_Output::print_at_end() const{

    if (_treattest){
        std::cout <<"#results cancer recurrence: <ratio> <recurrences> <total. diag.> <nolsc_ratio> <nolsc_recurrences> <no_lscdiags>"<<std::endl;
        if (_output_specifier==4) std::cout <<"# ";
        std::cout <<_recurrence_count/double(_no_recurrence_patients)
         <<" "<<_recurrence_count   <<" "  << _no_recurrence_patients<< " "<<_nolsc_recurrence_count/double(_diagnosed_nolsc)<<" "<<_nolsc_recurrence_count <<" "<< _diagnosed_nolsc<< std::endl;
    }

    std::cout << "#Real time elapsed in seconds: " << ((double)clock()-_timer)/CLOCKS_PER_SEC << std::endl;
    std::cout << "#Average time to diagnosis " << (_diagnosed > 0?(_total_diagnosis_time / (double) _diagnosed):0) << std::endl;
    std::cout << "#Fraction diagnosed "<< (_diagnosed /(double) _patients) << std::endl;
    std::cout << "#Fraction with no LSC " << (_nolsc / (double) _patients) << std::endl;
    std::cout << "#Fraction diagnosed with no LSC " << (_diagnosed > 0?(_diagnosed_nolsc / (double) _diagnosed):0) << std::endl;
    std::cout << "#Fraction that reached (default 4.5) "
        << " log reduction : <reduction freq.> <#of reductions> <diagnosed> <noscl at dignose>" << std::endl;
    if (_output_specifier!=3) std::cout <<"# ";
    std::cout << ((_reachedreduction > 0&&_diagnosed>0)?(_reachedreduction / (double) _diagnosed):0)
        << " " << _reachedreduction << " " << _diagnosed << " " << _diagnosed_nolsc<< std::endl;
    double stddev = 0;
    double avg = (_reachedreduction > 0?(_total_timetoreduction / (double) _reachedreduction):0) ;
    for(unsigned int i=0; i < _redresult.size(); i++){
        stddev += pow((double)(_redresult[i] - avg), 2.0);
    }
    stddev = stddev / (double)_redresult.size();
    stddev = sqrt(stddev);
    cout << "#Average time to reduction "<< avg << "\t" << stddev << endl;
    for (int i=0; i< _avgsize.size(); i++) {
        cout << "#avg size comp " << i << " = " << (_avgsize[i]/(double)_patients) << endl;
    }

}
