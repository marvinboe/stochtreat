#include "data.h"
#include "rangen.h"
#include "kernel.h"
#include "output.h"
#include "parameter_handler.h"

int main (int argc, char *argv[]) {

    int runid(1);
    int size(7); // 1 means only the stem cell compartment
    float treatmenttime (20);
    float mass(70); //human mass
    float reduction(4.5);
    double relapse_logreduction(3.);
    unsigned patients(1);
    double collectinterval(30.); //how often data is collected
    std::string path("./outinput/");
    std::string inpath(path);
    Diff_probabilities diff_probs;
    double ntime(25.);// maximum simulation time in years
    std::string output; //"patient nolsctime diagtime initresponse fullburden"
    Run_modes run_mode;

    ParameterHandler parameters(argc,argv);


    parameters.SetValue("run",	"Run identifier (1)",	runid);
    parameters.SetValue("size",	"Number of stochastic compartments (7)",	size);
    parameters.SetValue("treattime", "Years of treatment (10 years)", treatmenttime);
    parameters.SetValue("mass",	"animal mass (70 kg)",	mass);
    parameters.SetValue("reduction", "Required reduction level (4.5 logs)", reduction);
    parameters.SetValue("patients", "Number of patients (1)", patients);
    parameters.SetValue("path", "Output path (./outinput/) ", path);
    parameters.SetValue("inpath", "Input path (./outinput/) ", inpath);
    parameters.SetValue("ntime", "Maximum simulation time in years(25)", ntime);
    parameters.SetValue("output", "Specifiy kind of output (). possible: 'patient,nolsctime,diagtime,initresponse,fullburden,nooverview,yearlyburden,relapsetime,3timepointsmedian,3timepointsfull'. Can be combined: 'output=x1;x2;etc'.", output);
    parameters.SetValue("treattest", "test the treatment", run_mode.treattest);
    parameters.SetValue("resistance", "introduce resistant cell at diagnosis in specified compartment or in lowest(100)", run_mode.resistance);
    parameters.SetValue("epsn", "change differentiation probability for healthy cells (0.85)", diff_probs.epsh);
    parameters.SetValue("epsc", "change differentiation probability for cancer cells (0.71)", diff_probs.epsc);
    parameters.SetValue("epsb", "change differentiation probability for bound cells (0.89)", diff_probs.epsb);

    parameters.print_help(std::cout);



    double Nbase(16.52861491); // Number of HSCs=Nbase*mass^(0.75)
    double Bbase(2.892507609); // division rate HSC: 1/tau; tau = 365.0/(Bbase*pow(mass,-0.25))
    double Sbase(0.034572078); // deterministic timestep dt=Sbase*pow(mass,0.25)
    double Lbase(8.643019616); // maximum simulation time Tmax=(Lbase*pow(mass,0.25)); elephant 8.54663017, human 8.643019616
    double factor(1.0); // maximum simulation time factor

    RanGen ran;
    Data data;
    data.initialize(mass, Nbase, Bbase, Sbase, (Lbase * factor), collectinterval, diff_probs);
    if (ntime > 0.){ //non-default
        data.setTmax(ntime);
    }
    data.set_treatment_rate(0.05); //sets the rate of new bound cell under treatment (per day)
    data.set_diagnosis_limit(12); // Diagnosis limit
    data.set_treatment_stop_reduction(reduction); //treatment stop 
    data.set_relapse_reduction(relapse_logreduction); //treatment stop 
    data.set_numstochcomps(size);
    data.set_maximum_treatment_duration(treatmenttime);
    	// cout << data << endl;

    Stats_Output out(output,size,run_mode);

    std::vector<unsigned> redresult;
    for(unsigned i=0 ; i < patients; ++i){
        out.initialize_per_patient(i);
        Kernel ker(ran, data, size);
        
        //make run without treatment until diagnosis (or time exceeded).
        double time = ker.execute(ran,0.0,DIAGNOSISRUN);
        out.save_data_after_diagnosisrun(ker,time);

        //#### check if diagnosis is reached
        if(ker.doctor().diagnosis_reached()) {
            if (run_mode.resistance==100) ker.introduce_immunity_inlowest();
            else if (run_mode.resistance>=0 && run_mode.resistance<=32)
                ker.introduce_resistance(run_mode.resistance);

            //start treatment until limit is reached or maxmum time of treatment has passed
            // cout << "#burden is " << ker.burden() << " reduction is " << ker.getReduction() << endl;
            time=ker.execute(ran,time,TREATMENTRUN);
            out.save_data_after_treatment(ker,time);

            if (run_mode.treattest && ker.doctor().reduction_reached()){
                ker.set_ntime(time+10.);
                ker.reset_treatment(ran,time);
                time=ker.execute(ran,time,RELAPSERUN); //look for diagnosis again
                out.save_data_after_relapse(ker,time);
            }

        }//######### everything for case of diagnosis done
        out.print_patient(ker);

    }//end loop over patients

    out.print_at_end();

    return 0;
}

