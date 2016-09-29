/*
 *  data.cpp
 *  HemaMass
 *
 *  Created by Tom Lenaerts on 15/08/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "data.h"
#include <boost/lexical_cast.hpp> 
#include <fstream>
#include <sstream>
#include <cmath>

#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>


double Data::mylog(double p1, double base){
	if(p1<=0)
		return 0.0;
	else {
		gsl_sf_result logp1,logbase;
		int status_p1=gsl_sf_log_e(p1,&logp1);
		if(status_p1 != GSL_SUCCESS){
			cout << "Kernel::mylog  GSL Error "<<  status_p1 << " for p1 =" << p1 << endl;
			exit(-1);
		}
		int status_base=gsl_sf_log_e(base,&logbase);  // this division will normalize the entropy value between 0 and 1.
		if(status_base != GSL_SUCCESS){
			cout << "Kernel::mylog GSL Error "<< status_base<< " for base =" << base << endl;
			exit(-1);
		}
		return logp1.val/logbase.val;
	}
}

Data::Data(){
	_mass = 70;
	_rbase=1.263;
	_dt=0.1;
	_N0=400.0;
	_tau=365.0;	
	_frac_csc=0.0;
	_numlsc = 0;
	//are the same accross mamals
	_epsh=0.85; //TODO in paper 0.85 0.8476;
	_epsc=0.72;
	_epsi=_epsc;
	_epsb=0.89;
        _tmax=25.;
	_age=25;
	_ncompartments=32;
	_stop=10.39;//DP RAT
	_reduction = 3;
	_limit=7;
	_additional=0;
	_threshold = 0.2;
	//for storing data
	_step=100;
	_ofcompartment="compartment.txt";
	_offinal="mammal.txt";
	_ofname = "result.txt";
	_hdlocation ="./";

	//not relevant for the moment
	// _p_csc=0;
	// _p_imm=0;
	// _treatment_rate=0.0;
	// _treatment=3.0;
	// _rcancer = 1.0;
}


void Data::calcFromMass(double mass, double N,double B, double T, double L, double nummonth){ 
	_mass = mass;
	
	double hsc=N * pow(mass, 0.75);
	double temp;
	modf(hsc,&temp); 	
	double diff=hsc-temp;
		
	if(diff >= 0.5)
		_N0 = (temp+1); // always round up to next complete HSC
	else _N0=temp;
	_frac_csc=1/_N0; //always start with 1 LSC
	_numlsc = 1;
	
	_tau = 365.0/(B*pow(mass,-0.25));	
	_dt=T*pow(mass,0.25);	
	
	_age=(L*pow(mass,0.25));
        _tmax=(L*pow(mass,0.25));
	
	_step=int(nummonth*(30./_dt)); //output per x months
	
	//are the same accross mamals or changd by command line
	_rbase=1.26;//TODO in paper 1.26// 1.263;
	_epsh=0.85;//TODO in paper 0.85// 0.8476;
	_epsc=0.72;
	_epsb=0.89; //imatinib;
	_epsi=_epsc;
	_ncompartments=32;
	_stop=12;//log value at diagnosis
	_reduction = 3; //required log reduction
	_treatment_rate=0.05;
	_limit=7;
	_additional=0;
	_threshold = 0.2;
	

	//not relevant for the moment TODO???
	// _p_csc=0;
	// _p_imm=0;
	_treatment=3.0;
	_rcancer = 1.0;
}

Data::Data(const Data& other){
	_dt=other.dt();
	_epsh=other.epsh();
	_epsc=other.epsc();
	_epsb=other.epsb();
	_epsi=other.epsi();
	_p_csc=other.p_csc();
	_p_imm=other.p_imm();
	_frac_csc=other.frac_csc();
	_numlsc = other.numlsc();
	_treatment_rate=other.treatment_rate();
        _tmax=other.getTmax_in_years();
	_age=other.age();
	_rbase=other.rbase();
	_ncompartments=other.ncompartments();
	_stop=other.stop();
	_reduction=other.reduction();
	_N0=other.N0();
	_tau=other.tau();
	_limit=other.nstochcomp();
	_threshold = other.threshold();
	_additional=other.additional();
	_step=other.step();
	_treatment=other.treatment_dur();
	_ofcompartment=other.ofcompartment();
	_offinal=other.offinal();
	_ofname = other.ofname();
	_hdlocation =other.storage();
	_rcancer = other.rcancer();
	_mass = other.mass();
}


std::ostream& Data::display(std::ostream& os){
	os << "#Inputdata_for_hematopoietic_model" << endl;
	os << "  calculated :: " << endl;
	os << "    N0 " << _N0 << endl;
	os << "    tau " << _tau << endl;
	os << "    dt " << _dt << endl;
	os << "    Tmax " << _tmax<<" iterations=" << _age << " years" << endl;
	os << "  command line :: " << endl;
	os << "    mass " << _mass << endl;
	os << "    limit " << _limit << endl;
	os << "    step " << _step << endl;
	os << "    mammal " << _offinal << endl;
	os << "    threshold " << _threshold << endl;
	os << "    frac_csc " << _frac_csc << endl;
	os << "    numlsc " << _numlsc << endl;
	os << "  fixed :: " << endl;
	os << "    rbase " << _rbase << endl;	
	os << "    epsh " << _epsh << endl;
	os << "    espc " << _epsc << endl;
	os << "    espb " << _epsb << endl;
	os << "    espi " << _epsi << endl;
	os << "    p_csc " << _p_csc << endl;
	os << "    p_imm " << _p_imm << endl;
	os << "    treatment_rate " << _treatment_rate << endl;
	os << "    ncompartment " << _ncompartments << endl;
	os << "    additional " << _additional << endl;
	os << "    treatment " << _treatment << endl;
	os << "    stop " << _stop << endl;
	os << "    reduction " << _reduction << endl;
	os << "    rcancer " << _rcancer << endl;
	os << "  data collection :: " << endl;	
	os << "    output 1  " << _ofcompartment << endl;
	os << "    output 3  " << _ofname << endl;
	os << "    storage  " << _hdlocation << endl;
	
	return os;
}

