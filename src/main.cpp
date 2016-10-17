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
#include "output.h"
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
    double ntime(25.);// maximum simulation time in years
    std::string output; //"patient nolsctime diagtime initresponse fullburden"
    Run_modes run_mode;

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
    parameters.SetValue("output", "Specifiy kind of output. possible: 'patient,nolsctime,diagtime,initresponse,fullburden,nooverview,yearlyburden,relapsetime'", output);
    parameters.SetValue("treattest", "test the treatment", run_mode.treattest);
    parameters.SetValue("resistance", "introduce resistant cell at diagnosis in specified compartment or in lowest(=100)", run_mode.resistance);

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

    Stats_Output out(output,size,run_mode);

    std::vector<unsigned> redresult;
    for(unsigned i=0 ; i < patients; ++i){
        out.initialize_per_patient(i);
        Kernel ker(ran, data, size);
        
        //make run without treatment until diagnosis (or time exceeded).
        double time = ker.execute(ran,0.0,false);
        out.save_data_after_diagnosisrun(ker,time);

        //#### check if diagnosis is reached
        if(ker.reachedDiagnosis()) {
            if (run_mode.resistance==100) ker.introduce_immunity_inlowest();
            else if (run_mode.resistance>=0 && run_mode.resistance<=32)
                ker.introduce_resistance(run_mode.resistance);

            //start treatment until limit is reached or maxmum time of treatment has passed
            // cout << "#burden is " << ker.burden() << " reduction is " << ker.getReduction() << endl;
            time=ker.execute(ran,time,true);
            out.save_data_after_treatment(ker,time);

            if (run_mode.treattest){
                ker.set_ntime(time+10.);
                ker.reset_treatment(ran,time);
                time=ker.execute(ran,time,false); //look for diagnosis again
                out.save_data_after_relapse(ker,time);
            }

        }//######### everything for case of diagnosis done
        out.print_patient(ker);

    }//end loop over patients

    out.print_at_end();

    return 0;
}


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
