#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>
#include "data.h"
#include "rangen.h"
#include "kernel.h"
#include "parameter_handler.h"

int main (int argc, char *argv[]) {

    int runid(1);
    int size(7); // 1 means only the stem cell compartment
    float treatmenttime (10);
    float mass(70); //human mass
    float reduction(4.5);
    unsigned patients(1);
    double months(1); //how often output is written
    std::string path("./outinput/");
    std::string inpath(path);
    double ntime(25.);
    int output_specifier(0); //0-> normal;-1 -> only results in the end; 1-> nolsctime data; 2-> time to diag. etc. ; 3-> same as 2, but only if reduction reached; 4-> initial treatment response; 5-> full doctor output of burden
    bool treattest=false;

    ParameterHandler parameters(argc,argv);


    parameters.SetValue("run",	"Run identifier (default 1)",	runid);
    parameters.SetValue("size",	"Number of stochastic compartments (default 7)",	size);
    parameters.SetValue("treattime", "Years of treatment (default 10 year)", treatmenttime);
    parameters.SetValue("mass",	"animal mass (default 70 kg)",	mass);
    parameters.SetValue("reduction", "Required reduction level (default 4.5 logs)", reduction);
    parameters.SetValue("patients", "Number of patients (default 1)", patients);
    parameters.SetValue("path", "Output path (default ./outinput/) ", path);
    parameters.SetValue("inpath", "Input path (default ./outinput/) ", inpath);
    parameters.SetValue("ntime", "Maximum simulation time (years, default 25)", ntime);
    parameters.SetValue("output", "Specifiy kind of output", output_specifier);
    parameters.SetValue("treattest", "test the treatment", treattest);

    parameters.print_help(std::cout);



    double Nbase(16.52861491); // Number of HSCs=Nbase*mass^(0.75)
    double Bbase(2.892507609); // division rate HSC: 1/tau; tau = 365.0/(Bbase*pow(mass,-0.25))
    double Sbase(0.034572078); // deterministic timestep dt=Sbase*pow(mass,0.25)
    double Lbase(8.643019616); // maximum simulation time Tmax=(Lbase*pow(mass,0.25)); elephant 8.54663017, human 8.643019616
    double factor(1.0); // maximum simulation time factor

    RanGen ran;
    Data data;
    data.calcFromMass(mass, Nbase, Bbase, Sbase, (Lbase * factor), months);
    if (ntime > 0.){ //non-default
        data.setTmax(ntime);
    }
    data.set_treatment_rate(0.05); //sets the rate of new bound cell under treatment (per day)
    data.setStop(12); // Diagnosis limit
    data.setReduction(reduction); //treatment stop 
    data.setLimit(size);
    data.setTreatment(treatmenttime);
    	// cout << data << endl;

    double nolsc = 0;
    double diagnosed_nolsc = 0;
    double total_diagnosis_time=0;
    double diagnosed=0;
    double reachedreduction = 0;
    double total_timetoreduction = 0;
    double avgsize[size+1];

    bool recurrence_run=false;
    int no_recurrence_patients=0;
    unsigned recurrence_count=0;
    unsigned nolsc_recurrence_count=0;
    bool nolsc_treattest=false;
    for (int i=0; i< size+1; i++) {
        avgsize[i]=0.0;
    }

    clock_t timer=clock();
    vector<unsigned> redresult;
    for(unsigned i=0 ; i < patients; ++i){
        if (output_specifier==0){
            cout << "#patient " << i << endl;
        }

        bool lsc_at_diagnosis=true;

        Kernel ker(ran, data, size);
        // ker.writeModel(std::cout);

        // //#################read compartment data from file##########
        // stringstream ssin;
        // ssin << path<< "patient-"<< runid << "-"<< i << ".txt";
        // std::ifstream input(ssin.str().c_str());
        // if(!input.is_open()){
        //     if (recurrence_run){
        //         std::cout << "# unable to open input file " << ssin.str() << std::endl;
        //         // std::cout << " exiting program " << std::endl;
        //         // exit(-1);
        //         continue;
        //     }
        // }
        // else{//input is open
        //     ker.readModel(input);
        //     // ker.writeModel(std::cout);
        //     input.close();
        //     recurrence_run=true;
        //     no_recurrence_patients+=1;
        // }
        // //#############end reading model data ######################
        
        //make run without treatment until diagnosis (or time exceeded).
        double time = ker.execute(ran,0.0,false);
        ker.addStochCompSizes(avgsize);
        if(!ker.hasLSC())
            nolsc +=1;

        //#### check if diagnosis is reached
        double timetoreduction=-1.;
        if(ker.reachedDiagnosis()) {
            diagnosed +=1;
            if (treattest) no_recurrence_patients++;

            total_diagnosis_time += time;
            if(!ker.hasLSC()){
                lsc_at_diagnosis=false;
                diagnosed_nolsc +=1;
                if (treattest){
                    nolsc_treattest=true;
                }
            }

            if (recurrence_run){
                if (output_specifier==1){
                    std::cout << ker.getDiagnosisTime() << "  " 
                        << ker.get_nolsctime() << endl;
                }
                continue; // end this if we only check for recurrence
            }

            //start treatment until limit is reached or maxmum time of treatment has passed
            // cout << "#burden is " << ker.burden() << " reduction is " << ker.getReduction() << endl;
            time=ker.execute(ran,time,true);
            // cout << "#burden is " << ker.burden() << " reduction is " << ker.getReduction() << endl;


            if(ker.reachedReduction()){
                reachedreduction +=1;
                timetoreduction=(ker.whenReduction() - ker.getDiagnosisTime());
                total_timetoreduction +=timetoreduction;
                redresult.push_back(ker.whenReduction() - ker.getDiagnosisTime());
                if (output_specifier==3){
                    cout << "#<years to diag.> <years to red.> <total> <nolsctime> "<<std::endl
                        << ker.getDiagnosisTime() << "  " 
                        << timetoreduction << "  "
                        << ker.whenReduction() << " "
                        << ker.get_nolsctime() << endl;
                }

            }

            if (treattest){

                ker.set_ntime(time+10.);
                ker.reset_treatment(ran,time);
                time=ker.execute(ran,time,false); //look for diagnosis again
                if(ker.reachedDiagnosis()) {
                    recurrence_count++;
                    if (nolsc_treattest) 
                        nolsc_recurrence_count++;
                }
                nolsc_treattest=false;
            }

            if (output_specifier==4){
                std::cout <<ker.initial_treatment_response();
                std::cout <<" "<<lsc_at_diagnosis;
                if (!treattest){
                    std::cout <<" "<<ker.doctor().get_tumor_burden();
                }
                if (treattest)
                    std::cout << " "<<ker.reachedDiagnosis();
                std::cout <<std::endl;
            }
            if (output_specifier==5){
                ker.print_full_doctors_report(std::cout);
            }
            //			cout << "Reduction is " << ker.getReduction() << endl;
            //
            // //##########write compartment data to file
            // stringstream ss;
            // ss << path<< "patient-"<< runid << "-"<< i << ".txt";
            // ofstream output(ss.str().c_str());
            // if(!output.is_open()){
            //     cout << " unable to open output file " << ss.str() << endl;
            //     cout << " exiting program " << endl;
            //     exit(-1);
            // }
            // ker.writeModel(output);
            // output.close();
            // //######end writing patient data to file
        }//######### everything for case of diagnosis done
        

        if (output_specifier==1){
            cout  << ker.get_nolsctime() << endl;
        }
        else if (output_specifier==2){
            std::cout << ker.getDiagnosisTime() << "  " 
                << timetoreduction << "  "
                << ker.whenReduction() << " "
                << ker.get_nolsctime() << endl;
        }
        if (recurrence_run&&output_specifier==1){
            std::cout << ker.getDiagnosisTime() << "  " 
                << ker.get_nolsctime() << endl;
        }

    }//end loop over patients

    if (treattest){
        std::cout <<"#results cancer recurrence: <ratio> <recurrences> <total. diag.> <nolsc_ratio> <nolsc_recurrences> <no_lscdiags>"<<std::endl;
        if (output_specifier==4) std::cout <<"# ";
        std::cout <<recurrence_count/double(no_recurrence_patients)
         <<" "<<recurrence_count   <<" "  << no_recurrence_patients<< " "<<nolsc_recurrence_count/double(diagnosed_nolsc)<<" "<<nolsc_recurrence_count <<" "<< diagnosed_nolsc<< std::endl;
    }

    std::cout << "#Real time elapsed in seconds: " << ((double)clock()-timer)/CLOCKS_PER_SEC << std::endl;
    std::cout << "#Average time to diagnosis " << (diagnosed > 0?(total_diagnosis_time / (double) diagnosed):0) << endl;
    std::cout << "#Fraction diagnosed "<< (diagnosed /(double) patients) << std::endl;
    std::cout << "#Fraction with no LSC " << (nolsc / (double) patients) << std::endl;
    std::cout << "#Fraction diagnosed with no LSC " << (diagnosed > 0?(diagnosed_nolsc / (double) diagnosed):0) << endl;
    std::cout << "#Fraction that reached "<< reduction 
        << " log reduction : <reduction freq.> <#of reductions> <diagnosed> <noscl at dignose>" << endl;
    if (output_specifier!=3) std::cout <<"# ";
    std::cout << ((reachedreduction > 0&&diagnosed>0)?(reachedreduction / (double) diagnosed):0)
        << " " << reachedreduction << " " << diagnosed << " " << diagnosed_nolsc<< std::endl;
    double stddev = 0;
    double avg = (reachedreduction > 0?(total_timetoreduction / (double) reachedreduction):0) ;
    for(unsigned int i=0; i < redresult.size(); i++){
        stddev += pow((double)(redresult[i] - avg), 2.0);
    }
    stddev = stddev / (double)redresult.size();
    stddev = sqrt(stddev);
    cout << "#Average time to reduction "<< avg << "\t" << stddev << endl;
    for (int i=0; i< size+1; i++) {
        cout << "#avg size comp " << i << " = " << (avgsize[i]/(double)patients) << endl;
    }


    return 0;
}


