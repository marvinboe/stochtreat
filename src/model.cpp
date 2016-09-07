/*
 *  Model.cpp
 *  stochmut
 *
 *  Created by Tom Lenaerts on 1/31/11.
 *  Copyright 2011 Universite Libre de Bruxelles. All rights reserved.
 *
 */

#include <cassert>
#include <cmath>
#include <iomanip>

#include "model.h"
#include "dependency.h"

#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>


double Model::myround(double val){
	double temp;
	modf(val,&temp); 	
	double diff=val - temp;
	if(diff >=0.5)
		return temp+1;
	return temp;
}

double Model::mylog(double p1, double base){
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



Model::Model(Data data, unsigned int ns):_numstoch(ns),_diagnosis(0),_nolsctime(-1.){ //ns = 1 is alleen de stem cell Model
	assert(_numstoch > 0);
	_numcomp = data.ncompartments()+1;
	
	unsigned int tlen = ((_numcomp-1)*4) + 3;
	_endstoch = ((_numstoch-1)*4) + 3;
	_compartments=new double[tlen];
	_previous=new double[tlen];
	
	_rates = new double[_numcomp];
	//initialize Model: 1. Stem Cell compartment
	setH(0,data.N0());
	setC(0,0);
	setI(0,0);//initially no resistant mutants
	setRate(0, 1.0/data.tau());
	storeH(0,0);
	storeC(0,0);
	storeI(0,0);
	
	
	//initialize Model: 2. other compartments
	double gamma =( (2.0*data.epsh()) / (2.0*data.epsh() - 1.0) ) / data.rbase();
	for(unsigned k=1; k < _numcomp; ++k){
		double tempHk=getH(k-1)*gamma;
		if(k == 1)
			tempHk*=(1.0/ (2.0*data.epsh()));
		double Hk=myround(tempHk);
		setH(k,Hk);
		setC(k,0);
		setI(k,0);
		setB(k,0);
		storeH(k,0);
		storeC(k,0);
		storeI(k,0);
		storeB(k,0);
		
		setRate(k, getRate(k-1)*data.rbase());
	}
	//initialize CML
	setH(0, getH(0)-data.numlsc());
	setC(0, data.numlsc());
}

Model::Model(const Model& other){
	_numstoch = other.numStoch();
	_numcomp = other.numComp();
	unsigned int tlen = ((_numcomp-1)*4) + 3;
	_endstoch = ((_numstoch-1)*4) + 3;
	_compartments=new double[tlen];
	_previous=new double[tlen];

	_rates = new double[_numcomp];
	
	for(unsigned k=0; k < _numcomp; ++k){
		setRate(k, other.getRate(k));
		setH(k, other.getH(k));
		setC(k, other.getC(k));
		setI(k, other.getI(k));
		storeH(k,0);
		storeC(k,0);
		storeI(k,0);
		if(k>0){
			setB(k, other.getB(k));
			storeB(k,0);
		}
	}
}

Model& Model::operator=(const Model& other){
	_numstoch = other.numStoch();
	_numcomp = other.numComp();
	unsigned int tlen = ((_numcomp-1)*4) + 3;
	_endstoch = ((_numstoch-1)*4) + 3;
	_compartments=new double[tlen];
	_rates = new double[_numcomp];
	
	for(unsigned k=0; k < _numcomp; ++k){
		setRate(k, other.getRate(k));
		setH(k, other.getH(k));
		setC(k, other.getC(k));
		setI(k, other.getI(k));
		storeH(k,0);
		storeC(k,0);
		storeI(k,0);
		if(k>0){
			setB(k, other.getB(k));
			storeB(k,0);
		}
	}
	return *this;
}


void Model::setRate(unsigned int k, double v){
	assert (k>=0 && k < _numcomp);
	_rates[k]=v;
}

double Model::getRate(unsigned k) const{
	assert (k>=0 && k < _numcomp);
	return _rates[k];
}

double Model::getN(unsigned int k) const{
	assert (k>=0 && k < _numcomp);
	double total=0;
	if(k>0){
		unsigned int tmp = ((k-1)*4)+3;
		total = _compartments[tmp]+_compartments[tmp+1]+_compartments[tmp+2]+_compartments[tmp+3];
	}
	else {
		total = _compartments[0]+_compartments[1]+_compartments[2];
	}
	return total;
}


void Model::setH(unsigned int k, double v){
	assert (k>=0 && k < _numcomp);
	unsigned int tmp(0);
	if(k>0){		
		tmp = ((k-1)*4)+3;
	}
	_compartments[tmp]=v;
}

double Model::getH(unsigned int k) const{
	assert (k>=0 && k < _numcomp);
	unsigned int tmp(0);
	if(k>0)
		tmp = ((k-1)*4)+3;
	return _compartments[tmp];
}

void Model::setC(unsigned int k, double v){
	assert (k>=0 && k < _numcomp);
	unsigned int tmp(0);
	if(k>0)
		tmp = ((k-1)*4)+3;
	_compartments[tmp+1]=v;
}
 
double Model::getC(unsigned int k) const{
	assert (k>=0 && k < _numcomp);
	unsigned int tmp(0);
	if(k>0)
		tmp = ((k-1)*4)+3;
	return _compartments[tmp+1];
}
 
void Model::setI(unsigned int k, double v){
	assert (k>=0 && k < _numcomp);
	unsigned int tmp(0);
	if(k>0)
		tmp = ((k-1)*4)+3;
	_compartments[tmp+2]=v;
}
 
double Model::getI(unsigned int k) const{
	assert (k>=0 && k < _numcomp);
	unsigned int tmp(0);
	if(k>0)
		tmp = ((k-1)*4)+3;
	return _compartments[tmp+2];
}

void Model::setB(unsigned int k, double v){ //no cured cells in stem-cell compartment
	assert (k>0 && k < _numcomp);
	unsigned int tmp = ((k-1)*4)+3;
	_compartments[tmp+3]=v;
}

double Model::getB(unsigned int k) const {
	assert (k>0 && k < _numcomp);
	unsigned int tmp= ((k-1)*4)+3;
	return _compartments[tmp+3];
}

double Model::get(unsigned k, unsigned t){
	switch(t){
		case H:return getH(k);
		case C:return getC(k);
		case I:return getI(k);
		case B:return getB(k);
	}
	return -1.0;
}

void Model::set(unsigned k, unsigned t, double v){
	switch(t){
		case H:
			setH(k,v);
			break;
		case C:
			setC(k,v);
			break;
		case I:
			setI(k,v);
			break;
		case B:
			setB(k,v);
			break;
	}
}

void Model::incr(unsigned k, unsigned t, double v){
	switch(t){
		case H:
			setH(k,getH(k)+v);
			break;
		case C:
			setC(k,getC(k)+v);
			break;
		case I:
			setI(k,getI(k)+v);
			break;
		case B:
			setB(k,getB(k)+v);
			break;
	}
}

void Model::decr(unsigned k, unsigned t, double v){
	switch(t){
		case H:
			if(getH(k)>0)
				setH(k,getH(k)-v);
			break;
		case C:
			if(getC(k)>0)
				setC(k,getC(k)-v);
			break;
		case I:
			if(getI(k)>0)
				setI(k,getI(k)-v);
			break;
		case B:
			if(getB(k)>0)
				setB(k,getB(k)-v);
			break;
	}
}

double Model::retrieveN(unsigned int k) const{
	assert (k>=0 && k < _numcomp);
	double total=0;
	if(k>0){
		unsigned int tmp = ((k-1)*4)+3;
		total = _previous[tmp]+_previous[tmp+1]+_previous[tmp+2]+_previous[tmp+3];
	}
	else {
		total = _previous[0]+_previous[1]+_previous[2];
	}
	return total;
}


void Model::storeH(unsigned int k, double v){
	assert (k>=0 && k < _numcomp);
	unsigned int tmp(0);
	if(k>0){		
		tmp = ((k-1)*4)+3;
	}
  	_previous[tmp]=v;
}

double Model::retrieveH(unsigned int k) const{
	assert (k>=0 && k < _numcomp);
	unsigned int tmp(0);
	if(k>0)
		tmp = ((k-1)*4)+3;
	return _previous[tmp];
}

void Model::storeC(unsigned int k, double v){
	assert (k>=0 && k < _numcomp);
	unsigned int tmp(0);
	if(k>0)
		tmp = ((k-1)*4)+3;
	_previous[tmp+1]=v;
}

double Model::retrieveC(unsigned int k) const{
	assert (k>=0 && k < _numcomp);
	unsigned int tmp(0);
	if(k>0)
		tmp = ((k-1)*4)+3;
	return _previous[tmp+1];
}

void Model::storeI(unsigned int k, double v){
	assert (k>=0 && k < _numcomp);
	unsigned int tmp(0);
	if(k>0)
		tmp = ((k-1)*4)+3;
	_previous[tmp+2]=v;
}

double Model::retrieveI(unsigned int k) const{
	assert (k>=0 && k < _numcomp);
	unsigned int tmp(0);
	if(k>0)
		tmp = ((k-1)*4)+3;
	return _previous[tmp+2];
}

void Model::storeB(unsigned int k, double v){ //no cured cells in stem-cell compartment
	assert (k>0 && k < _numcomp);
	unsigned int tmp = ((k-1)*4)+3;
	_previous[tmp+3]=v;
}

double Model::retrieveB(unsigned int k) const {
	assert (k>0 && k < _numcomp);
	unsigned int tmp = ((k-1)*4)+3;
	return _previous[tmp+3];
}

double Model::retrieve(unsigned k, unsigned t){
	switch(t){
		case H:return retrieveH(k);
		case C:return retrieveC(k);
		case I:return retrieveI(k);
		case B:return retrieveB(k);
	}
	return -1.0;
}

void Model::store(unsigned k, unsigned t, double v){
	switch(t){
		case H:
			storeH(k,v);
			break;
		case C:
			storeC(k,v);
			break;
		case I:
			storeI(k,v);
			break;
		case B:
			storeB(k,v);
			break;
	}
}

ostream& Model::display(ostream& os){
    os <<_numcomp<<std::endl;
    for(unsigned k=0; k < _numcomp; ++k){
        os << k <<" " << setprecision(6) << getRate(k) <<" "; 
        os << getN(k) << " " << getH(k) <<" " << getC(k) <<" " << getI(k);
        if (k>0) os  <<" "<< getB(k);
        os << std::endl; 
    }
    os << "# " << _numstoch << endl;
    os << "# " << _numcomp << endl;
    os << "# " << _alpha << endl;
    os << "# " << _diagnosis << endl;
    os << "# " << getReduction() << endl;
    os << "# " << when() << endl;
    return os;
}

std::istream& Model::read(std::istream& is){
    is >> _numcomp;
    for(unsigned i=0; i < _numcomp; ++i){
        unsigned int k;
        double N,H,C,I,B;
        double rate;
        is >> k >> rate ;
        is >> N >>  H >> C >> I;
        if (k > 0) {
            is >> B; 
            setB(k,B);
        }
        setRate(k,rate);
        setH(k,H);
        setC(k,C);
        setI(k,I);
    }
    calcAlpha();
    return is;
}

void Model::memorize(){
	for(unsigned k=0 ; k < _numcomp; ++k){
		storeH(k,getH(k));
		storeC(k,getC(k));
		storeI(k,getI(k));
		if(k>0)
			storeB(k,getB(k));
		if(k >= _numstoch){
			setH(k, 0.0);
			setC(k, 0.0);
			setI(k, 0.0);
			if(k>0)
				setB(k,0.0);
		}
	}
}


bool Model::updateDet(unsigned k, Data& data){	
	double p = data.dt() *getRate(k);
	double prev_p = data.dt() *getRate(k-1); 
	double prev_epsh = data.epsh(); 

//	cout << k-1 << " : " << getN(k-1) << "\t" << getH(k-1) << "\t" << getC(k-1) << endl;
//	cout << k << " : " << getN(k) << "\t" << getH(k) << "\t" << getC(k) << "\t" << retrieveN(k) << "\t" << retrieveH(k) << "\t" << retrieveC(k) << endl;
//	cout << p << "\t" << data.epsh() << "\t" << data.epsc() << endl;
	
	if(k<=_numstoch){
//		cout << "#influx " << (2.0 * prev_epsh * prev_p * retrieveH(k-1)) << endl;
		prev_epsh = 0.0;
	}
	double tempH= retrieveH(k) + (2.0 * prev_epsh * prev_p * retrieveH(k-1)) + 
				(p * (1.0 - data.epsh()) * retrieveH(k)) - 
				(p * data.epsh() * retrieveH(k));
//	cout << "#current H("<<k<<")=" << getH(k) << endl; 
	incr(k, H,tempH);
//	cout << "#new H " << getH(k) << endl; 
	
	double prev_epsc = data.epsc();
	if(k<=_numstoch)
		prev_epsc = 0.0;
    double tempC= retrieveC(k) + (2.0 * prev_epsc * prev_p * retrieveC(k-1)) +
			(p * (1.0 - data.epsc()) * retrieveC(k)) -
			(p * data.epsc() * retrieveC(k));
	incr(k, C,tempC);
	
	double prev_epsi = data.epsi();
	if(k<=_numstoch)
		prev_epsi = 0.0;
	double tempI= retrieveI(k) + (2.0 * prev_epsi  * prev_p * retrieveI(k-1)) +
			(p * (1.0 - data.epsi()) * retrieveI(k-1)) -
			(p * data.epsi() * retrieveI(k-1));
	incr(k, I, tempI);
	
	double prev_epsb = data.epsb();
	if(k<=_numstoch || k == 1 )
		prev_epsb = 0.0;
	double tempB= retrieveB(k) + (2.0 * prev_epsb * prev_p * ((k-1)>0?retrieveB(k-1):0.0)) +
			(p * (1.0 - data.epsb()) * retrieveB(k)) -
			(p * data.epsb() * retrieveB(k));
	incr(k, B,tempB);
//	cout << k << " : " << getN(k) << "\t" << getH(k) << "\t" << getC(k) << "\t" << getI(k) << "\t" << getB(k) << endl;
	return true;
}

bool Model::treatDeterministically(unsigned k, double amount){	
	double ccells = getC(k);
	double bcells = getB(k);
	double tmp = ccells * amount;
	setC(k, ccells - tmp);
	setB(k, bcells + tmp);
	return true;
}

bool Model::treatStochastically(unsigned k, double rate, RanGen& ran){	
	double ccells = getC(k);
	double bcells = getB(k);

	double tmp (0);
	for(unsigned i = 0; i < ccells; ++i){
		if(ran.randouble() < rate) ++tmp;
	}
	
	bool changed = (tmp>0)?true:false;
	if (changed){
		setC(k, ccells - tmp);
		setB(k, bcells + tmp);
//		cout << "#Fraction changed in " << k << " : " << setprecision(3)<< (tmp/ccells)*100 << "% (" << ccells << " , " << getC(k)<< " , " << tmp << " , " << bcells << " , " << getB(k)<< ")" << endl;
	}
	return changed;	
}


bool Model::diagnosis(Data& data){
	// double res = mylog(getN(_numcomp-1),10);
	// cout <<  res<<" "<<mylog(lastN(),10) << "\t" << data.stop() << endl;
	return mylog(lastN(),10)>= data.stop();
}

bool Model::reduction(Data& data){
	return getReduction() >= data.reduction();
}

float Model::getReduction(){
	double b = diseaseBurden();
	return (2.0 - mylog(b,10));
}


bool Model::containsLSC(){
	return (getC(0) > 0);
}

void Model::calcAlpha(){
	double NC=lastC();
	double NB=lastB();
	double NH=lastH();
	_alpha = (NC + NB + (2.0 * NH)) / (NC + NB);
//	cout << "Alpha = " << _alpha << endl;
}

double Model::diseaseBurden(){
	double NC=lastC();
	double NB=lastB();
	double NH=lastH();
	double burden = (_alpha*(100.0 * ((NC + NB) / (NC + NB + (2.0 * NH)))));
//	cout << burden << endl;
	return burden;
}


void Model::print_cells(std::ostream & os,double _time){
    os <<_time/365.0<<" ";
    os <<getH(0)<<" ";
    os <<getC(0)<<" ";
    os <<getI(0)<<" ";
    os <<std::endl;
}

void Model::check_LSCvanished(double t){
    if (_nolsctime>0.)
        return;
    if (getC(0)==0)
        _nolsctime=t/365.;
}




