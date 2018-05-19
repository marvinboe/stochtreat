/*
 *  data.cpp
 *  HemaMass
 *
 *  Created by Tom Lenaerts on 15/08/06.
 *  Modified by Marvin A. Böttcher 2017
 *
 *
 */

#include "data.h"

Data::Data(){
    _mass = -1.;
    _dt=-1.;
    _N0=400.0;
    _frac_csc=0.0;
    _numlsc = 0;
    //are the same accross mamals
    _diffprobs.epsh=0.85; //in paper 0.85 //0.8476;
    _diffprobs.epsc=0.72;
    _diffprobs.epsr=_diffprobs.epsc;
    _diffprobs.epsb=0.89;
    _relapse_waiting_time=10.;
    _tmax=25.;
    _ncompartments=32;
    _diagnosis_amplification=2.875; //10.39;//DP RAT
    _reduction = 3;
    _numstochcomps=7;
    _additional=0;
    _threshold = 0.2;
    //for storing data
    _outputstep=100;

    //not relevant for the moment
    // _p_csc=0;
    // _p_imm=0;
    _treatment_duration=10.0;
}


void Data::initialize(const Simulation_Parameters & simparams){ 

    double Nbase(16.52861491); // Number of HSCs=Nbase*mass^(0.75)
    // double Bbase(2.892507609); // division rate HSC: 1/tau; tau = 365.0/(Bbase*pow(mass,-0.25))
    double Sbase(0.034572078); // deterministic timestep dt=Sbase*pow(mass,0.25)
    double Lbase(8.643019616); // maximum simulation time Tmax=(Lbase*pow(mass,0.25)); elephant 8.54663017, human 8.643019616

    _simparams=simparams;
    _mass = _simparams.mass;

    if (_mass > 0.){
        double hsc=Nbase * pow(_mass, 0.75);
        _N0=std::round(hsc);//round to integer number of HSC
    }
    else {
        _N0=std::round(_simparams.n_hsc);
    }

    _numlsc = _simparams.inital_lsc; //
    _frac_csc=_numlsc/_N0;

    _prolif=_simparams.prolif;
    _prolif.set_unset_params();


    _dt=_simparams.dt;
    if (_dt<0.){
        _dt=Sbase*pow(_mass,0.25);	
    }

    // _age=(L*pow(mass,0.25));
    _tmax=(Lbase*pow(_mass,0.25));

    _outputstep=int(_simparams.collectinterval/_dt); //output_steps

    //are the same accross mammals or changd by command line
    _diffprobs=_simparams.diff_probs;
    _diffprobs.epsr=_diffprobs.epsc;
    if (_diffprobs.beta<0) _diffprobs.beta=1-_diffprobs.eps_asym;

    _ncompartments=_simparams.n_compartments;
    _n_neutral_compartments=_simparams.n_neutral;
    _diagnosis_amplification=_simparams.diagnosis_amplification;
    _reduction = _simparams.reduction; //required log reduction
    _required_redtime = _simparams.required_reduction_time*365.; //required log reduction
    _treatment_rate=_simparams.treatment_rate;
    _numstochcomps=_simparams.n_stochastic_compartments;


    set_relapse_reduction(_simparams.relapse_logreduction); //treatment stop 
    set_maximum_treatment_duration(_simparams.treatmenttime);
    _relapse_waiting_time=_simparams.relapse_waiting_time;
    if (_simparams.ntime > 0.){ //non-default
        setTmax(_simparams.ntime);
    }
}

std::ostream& Data::display(std::ostream& os){
    os << "#Inputdata_for_hematopoietic_model" << std::endl;
    os << "    N0 " << _N0 << std::endl;
    os << "    dt " << _dt << std::endl;
    os << "    Tmax " << _tmax<<" years"<<std::endl;
    os << "    mass " << _mass << std::endl;
    os << "    #stoch. comp. " << _numstochcomps << std::endl;
    os << "    step " << _outputstep << std::endl;
    os << "    threshold " << _threshold << std::endl;
    os << "    frac_csc " << _frac_csc << std::endl;
    os << "    numlsc " << _numlsc << std::endl;
    os << "    epsh " << _diffprobs.epsh << std::endl;
    os << "    espc " << _diffprobs.epsc << std::endl;
    os << "    espb " << _diffprobs.epsb << std::endl;
    os << "    espr " << _diffprobs.epsr << std::endl;
    os << "    esp_asym " << _diffprobs.eps_asym << std::endl;
    os << "    beta " << _diffprobs.beta << std::endl;
    os << "    p_csc " << _p_csc << std::endl;
    os << "    p_imm " << _p_imm << std::endl;
    os << "    treatment_rate " << _treatment_rate << std::endl;
    os << "    ncompartment " << _ncompartments << std::endl;
    os << "    additional " << _additional << std::endl;
    os << "    treatment duration " << _treatment_duration << std::endl;
    os << "    stop " << _diagnosis_amplification << std::endl;
    os << "    reduction " << _reduction << std::endl;

    return os;
}


void Simulation_Parameters::set_parameters(ParameterHandler & parameters){
    parameters.SetValue("run",	"Run identifier (1)",	runid);
    parameters.SetValue("size",	"Number of stochastic compartments (7)",	n_stochastic_compartments);
    parameters.SetValue("stochcomps",	"Number of stochastic compartments (7)",	n_stochastic_compartments);
    parameters.SetValue("n_hsc", "Total size of hematopoietic stem cell compartment (400)", n_hsc);
    parameters.SetValue("compartments",	"Total number of compartments (32)",	n_compartments);
    parameters.SetValue("n_neutral","Number of compartments where bcr/abl is neutral (1)",	n_neutral);
    parameters.SetValue("lsc", "Initial number of leukemic stem cells (1)", inital_lsc);
    parameters.SetValue("treattime", "Maximum years of treatment (10 years)", treatmenttime);
    parameters.SetValue("treatrate", "rate at which cells are bound to drug (0.05/day)", treatment_rate);
    parameters.SetValue("mass",	"obsolete: animal mass (-1)",	mass);
    parameters.SetValue("reduction", "Required reduction level (4.5 log)", reduction);
    parameters.SetValue("diagnosis", "Amplification of cell count in last compartment for diagnosis (2.857x)", diagnosis_amplification);
    parameters.SetValue("reductiontime", "Required time in reduction (2 years)", required_reduction_time);
    parameters.SetValue("relapse_reduction", "Required reduction level (3 log)", relapse_logreduction);
    parameters.SetValue("patients", "Number of patients (1)", patients);

    parameters.SetValue("ntime", "Maximum simulation time (25 years)", ntime);
    parameters.SetValue("timestep", "Timestep for integration of deterministic equations (0.1 days)", dt);
    parameters.SetValue("relapse_waiting", "Waiting time for relapse after treatment (10 years)", relapse_waiting_time);

    parameters.SetValue("resistance", "introduce resistant cell at diagnosis in specified compartment or in lowest=100 (-1)", run_mode.resistance);
    parameters.SetValue("epsn", "change differentiation probability for healthy cells (0.85)", diff_probs.epsh);
    parameters.SetValue("epsc", "change differentiation probability for cancer cells (0.71)", diff_probs.epsc);
    parameters.SetValue("epsb", "change differentiation probability for bound cells (0.89)", diff_probs.epsb);
    parameters.SetValue("epsr", "change differentiation probability for resistant cells (0.71)", diff_probs.epsr);
    parameters.SetValue("eps_asym", "asymmetric division probability stem cells (0.0)", diff_probs.eps_asym);
    parameters.SetValue("beta", "stem cell replacement probability (1-eps_asym)", diff_probs.beta);

    parameters.SetValue("r0n", "base proliferation rate of stem cells (1/365 per day)", prolif.r0n());
    parameters.SetValue("gamman", "proliferation rate expansion between comps (1.263)", prolif.gamman());
    parameters.SetValue("r0c", "base proliferation rate of cancer cells (r0n)", prolif.r0c());
    parameters.SetValue("gammac", "proliferation rate expansion cancer cells(gamman)", prolif.gammac());
    parameters.SetValue("r0b", "base proliferation rate bound cells (r0n)", prolif.r0b());
    parameters.SetValue("gammab", "proliferation rate expansion bound cells (gamman)", prolif.gammab());

    parameters.SetValue("output", "Specifiy kind of output (). possible: 'patient,nolsctime,\
            diagtime,reductiontime,initresponse,fullburden,nooverview,yearlyburden,lastburden,\
            relapsetime,timepointburdenmedian,timepointburdenfull,treatdynamics,relapsedynamics,\
            treatdynall,relapsedynall,print_simparams'. Can be combined: 'output=x1;x2;etc'.", output);
    parameters.SetValue("treattest", "test the treatment", run_mode.treattest);

}
