/*
 *  data.cpp
 *  HemaMass
 *
 *  Created by Tom Lenaerts on 15/08/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "data.h"

Data::Data(){
	_mass = 70;
	_dt=0.1;
	_N0=400.0;
	_frac_csc=0.0;
	_numlsc = 0;
	//are the same accross mamals
	_diffprobs.epsh=0.85; //in paper 0.85 //0.8476;
	_diffprobs.epsc=0.72;
	_diffprobs.epsr=_diffprobs.epsc;
	_diffprobs.epsb=0.89;
        _tmax=25.;
	_ncompartments=32;
	_diagnosis_level=10.39;//DP RAT
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


void Data::initialize(const Simulation_Parameters & simparams, double N,double B, double T, double L){ 
	_mass = simparams.mass;
	
	double hsc=N * pow(_mass, 0.75);
        _N0=std::round(hsc);//round to integer number of HSC

	_numlsc = 1; //always start with 1 LSC
	_frac_csc=_numlsc/_N0;
	
        _prolif=simparams.prolif;

	_dt=T*pow(_mass,0.25);	
	
	// _age=(L*pow(mass,0.25));
        _tmax=(L*pow(_mass,0.25));
	
	_outputstep=int(simparams.collectinterval/_dt); //output_steps
	
	//are the same accross mamals or changd by command line
	_diffprobs.epsh=simparams.diff_probs.epsh;//in paper 0.85// 0.8476;
	_diffprobs.epsc=simparams.diff_probs.epsc;
	_diffprobs.epsb=simparams.diff_probs.epsb; //imatinib;
	_diffprobs.epsr=_diffprobs.epsc;

	_ncompartments=simparams.n_compartments;
	_diagnosis_level=simparams.diagnosis_level;//log value at diagnosis
	_reduction = simparams.reduction; //required log reduction
	_required_redtime = simparams.required_reduction_time*365.; //required log reduction
	_treatment_rate=simparams.treatment_rate;
	_numstochcomps=simparams.n_stochastic_compartments;


        set_relapse_reduction(simparams.relapse_logreduction); //treatment stop 
        set_maximum_treatment_duration(simparams.treatmenttime);
        if (simparams.ntime > 0.){ //non-default
            setTmax(simparams.ntime);
        }
}

Data::Data(const Data& other){
	_dt=other.dt();
	_diffprobs.epsh=other.epsh();
	_diffprobs.epsc=other.epsc();
	_diffprobs.epsb=other.epsb();
	_diffprobs.epsr=other.epsi();
        _prolif=other.return_prolif_params();
	_p_csc=other.p_csc();
	_p_imm=other.p_imm();
	_frac_csc=other.frac_csc();
	_numlsc = other.numlsc();
	_treatment_rate=other.treatment_rate();
        _tmax=other.getTmax_in_years();
	_ncompartments=other.ncompartments();
	_diagnosis_level=other.diagnosis_level();
	_reduction=other.reduction();
        _required_redtime=other.required_reduction_time();
	_N0=other.N0();
	_numstochcomps=other.nstochcomp();
	_threshold = other.threshold();
	_additional=other.additional();
	_outputstep=other.step();
	_treatment_duration=other.treatment_dur();
	_mass = other.mass();
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
	os << "    p_csc " << _p_csc << std::endl;
	os << "    p_imm " << _p_imm << std::endl;
	os << "    treatment_rate " << _treatment_rate << std::endl;
	os << "    ncompartment " << _ncompartments << std::endl;
	os << "    additional " << _additional << std::endl;
	os << "    treatment duration " << _treatment_duration << std::endl;
	os << "    stop " << _diagnosis_level << std::endl;
	os << "    reduction " << _reduction << std::endl;
	
	return os;
}


void Simulation_Parameters::set_parameters(ParameterHandler & parameters){
    parameters.SetValue("run",	"Run identifier (1)",	runid);
    parameters.SetValue("size",	"Number of stochastic compartments (7)",	n_stochastic_compartments);
    parameters.SetValue("stochcomps",	"Number of stochastic compartments (7)",	n_stochastic_compartments);
    parameters.SetValue("n_compartments",	"Total number of compartments (32)",	n_compartments);
    parameters.SetValue("treattime", "Maximum years of treatment (10 years)", treatmenttime);
    parameters.SetValue("treatrate", "rate at which cells are bound to drug (0.05/day)", treatment_rate);
    parameters.SetValue("mass",	"animal mass (70 kg)",	mass);
    parameters.SetValue("reduction", "Required reduction level (4.5 log)", reduction);
    parameters.SetValue("reductiontime", "Required time in reduction (0 years)", required_reduction_time);
    parameters.SetValue("relapse_reduction", "Required reduction level (3 log)", relapse_logreduction);
    parameters.SetValue("patients", "Number of patients (1)", patients);
    parameters.SetValue("ntime", "Maximum simulation time (25 years)", ntime);

    parameters.SetValue("resistance", "introduce resistant cell at diagnosis in specified compartment or in lowest=100 (-1)", run_mode.resistance);
    parameters.SetValue("epsn", "change differentiation probability for healthy cells (0.85)", diff_probs.epsh);
    parameters.SetValue("epsc", "change differentiation probability for cancer cells (0.71)", diff_probs.epsc);
    parameters.SetValue("epsb", "change differentiation probability for bound cells (0.89)", diff_probs.epsb);

    parameters.SetValue("kn", "base proliferation rate of stem cells (1/365 per day)", prolif.kn);
    parameters.SetValue("gamman", "proliferation rate expansion between comps (1.263)", prolif.gamman);
    parameters.SetValue("kc", "base proliferation rate of cancer cells(1/365 per day)", prolif.kc);
    parameters.SetValue("gammac", "proliferation rate expansion cancer cells(1.263)", prolif.gammac);
    parameters.SetValue("kb", "base proliferation rate bound cells (1/365 per day)", prolif.kb);
    parameters.SetValue("gammab", "proliferation rate expansion bound cells (1.263)", prolif.gammab);

    parameters.SetValue("output", "Specifiy kind of output (). possible: 'patient,nolsctime,diagtime,reductiontime,initresponse,fullburden,nooverview,yearlyburden,relapsetime,3timepointsmedian,3timepointsfull,treatdynamics'. Can be combined: 'output=x1;x2;etc'.", output);
    parameters.SetValue("treattest", "test the treatment", run_mode.treattest);

}
