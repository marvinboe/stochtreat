#include "data.h"
#include "rangen.h"
#include "kernel.h"
#include "output.h"
#include "parameter_handler.h"

int main (int argc, char *argv[]) {

    Simulation_Parameters simparams;

    ParameterHandler parameters(argc,argv);
    simparams.set_parameters(parameters);
    parameters.print_help(std::cout);



    double Nbase(16.52861491); // Number of HSCs=Nbase*mass^(0.75)
    // double Bbase(2.892507609); // division rate HSC: 1/tau; tau = 365.0/(Bbase*pow(mass,-0.25))
    double Sbase(0.034572078); // deterministic timestep dt=Sbase*pow(mass,0.25)
    double Lbase(8.643019616); // maximum simulation time Tmax=(Lbase*pow(mass,0.25)); elephant 8.54663017, human 8.643019616
    double factor(1.0); // maximum simulation time factor ... WHY?

    RanGen ran;
    Data data;
    data.initialize(simparams, Nbase, Sbase, (Lbase * factor));
    // std::cout << data << std::endl;

    Stats_Output out(simparams.output,simparams.n_stochastic_compartments,simparams.run_mode);

    std::vector<unsigned> redresult;
    for(unsigned i=0 ; i < simparams.patients; ++i){
        out.initialize_per_patient(i);
        Kernel ker(ran, data);
        
        //make run without treatment until diagnosis (or time exceeded).
        double time = ker.execute(ran,0.0,DIAGNOSISRUN);
        out.save_data_after_diagnosisrun(ker,time);

        //#### check if diagnosis is reached
        if(ker.doctor().diagnosis_reached()) {
            if (simparams.run_mode.resistance==100) ker.introduce_immunity_inlowest();
            else if (simparams.run_mode.resistance>=0 && simparams.run_mode.resistance<=32)
                ker.introduce_resistance(simparams.run_mode.resistance);

            //start treatment until limit is reached or maxmum time of treatment has passed
            // cout << "#burden is " << ker.burden() << " reduction is " << ker.getReduction() << endl;
            time=ker.execute(ran,time,TREATMENTRUN);
            out.save_data_after_treatment(ker,time);

            if (simparams.run_mode.treattest && ker.doctor().reduction_reached()){
                ker.set_ntime(time+data.relapse_waiting_time());
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

