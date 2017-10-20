/*
 *  Model.cpp
 *  stochmut
 *
 *  Created by Tom Lenaerts on 1/31/11.
 *  Copyright 2011 Universite Libre de Bruxelles. All rights reserved.
 *
 */


#include "model.h"



double Model::myround(double val){
    double temp;
    modf(val,&temp); 	
    double diff=val - temp;
    if(diff >=0.5)
        return temp+1;
    return temp;
}

double Model::mylog(double p1, double base) const{
    if(p1<=0)
        return 0.0;
    else {
        return std::log(p1)/std::log(base);
    }
}



Model::Model(Data data):_diagnosis(0){
    _numstoch=data.nstochcomp();
    assert(_numstoch > 0);
    _numcomp = data.ncompartments()+1;
    _numtypes=4;

    unsigned int tlen = ((_numcomp-1)*4) + 3;
    _endstoch = ((_numstoch-1)*4) + 3;
    _compartments=new double[tlen];
    _previous=new double[tlen];

    for (unsigned k=0 ; k < _numcomp; ++k){ //TODO that can be be done nicely
        std::vector<double> bufvec(4,0.);
        _rates.push_back(bufvec);
    }
    //initialize Model: 1. Stem Cell compartment
    setH(0,data.N0());
    setC(0,0);
    setI(0,0);//initially no resistant mutants
    for (unsigned t = 0; t< _numtypes; ++t){
        setRate(0,t, data.base_proliferation(t));
    }
    storeH(0,0);
    storeC(0,0);
    storeI(0,0);

    setTreatRate(data.treatment_rate());


    //initialize Model: 2. other compartments
    double gamma =( (2.0*data.epsh()) / (2.0*data.epsh() - 1.0) ) / data.prolif_exp(0); //to calculate steady state compartment sizes
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

        for (unsigned t = 0; t< _numtypes; ++t){
            // if (k<data.n_neutral_compartments()){
            //     setRate(k,t, getRate(k-1,0)*data.prolif_exp(t));//neutral compartments
            // }
            // else{ FIXME
                setRate(k,t, getRate(k-1,t)*data.prolif_exp(t));
            // }
        }
    }
    //initialize CML
    setH(0, getH(0)-data.numlsc());
    setC(0, data.numlsc());
}

Model::Model(Data data,std::istream & is):_diagnosis(0){ 
    _numcomp = data.ncompartments()+1;
    _numtypes=4;
    unsigned int tlen = ((_numcomp-1)*4) + 3;
    _compartments=new double[tlen];
    _previous=new double[tlen];
    for (unsigned k=0 ; k < _numcomp; ++k){ //TODO that can be be done nicely
        std::vector<double> bufvec(4,0.);
        _rates.push_back(bufvec);
    }

    read(is);
}

Model::Model(const Model& other){
    _numstoch = other.numStoch();
    _numcomp = other.numComp();
    _numtypes=4;
    unsigned int tlen = ((_numcomp-1)*4) + 3;
    _endstoch = ((_numstoch-1)*4) + 3;
    _compartments=new double[tlen];
    _previous=new double[tlen];

    for (unsigned k=0 ; k < _numcomp; ++k){ //TODO that can be be done nicely
        std::vector<double> bufvec(4,0.);
        _rates.push_back(bufvec);
    }

    for(unsigned k=0; k < _numcomp; ++k){
        for (unsigned t = 0; t< _numtypes; ++t){
            setRate(k,t, other.getRate(k,t));
            set(k,t,other.get(k,t));
        }
        // setH(k, other.getH(k));
        // setC(k, other.getC(k));
        // setI(k, other.getI(k));
        storeH(k,0);
        storeC(k,0);
        storeI(k,0);
        if(k>0){
        //     setB(k, other.getB(k));
            storeB(k,0);
        }
    }
}

Model& Model::operator=(const Model& other){
    _numstoch = other.numStoch();
    _numcomp = other.numComp();
    _numtypes=4;
    unsigned int tlen = ((_numcomp-1)*4) + 3;
    _endstoch = ((_numstoch-1)*4) + 3;
    _compartments=new double[tlen];
    for (unsigned k=0 ; k < _numcomp; ++k){ //TODO that can be be done nicely
        std::vector<double> bufvec(4,0.);
        _rates.push_back(bufvec);
    }

    for(unsigned k=0; k < _numcomp; ++k){
        for (unsigned t = 0; t< _numtypes; ++t){
            setRate(k,t, other.getRate(k,t));
        }
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


void Model::setRate(unsigned int k, unsigned t, double v){
    assert (k>=0 && k < _numcomp);
    assert (t>=0 && t < 4);
    _rates[k][t]=v;
}

double Model::getRate(unsigned k, unsigned t) const{
    assert (k>=0 && k < _numcomp);
    assert (t>=0 && t < 4);
    return _rates[k][t];
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

double Model::get(unsigned k, unsigned t) const{
    switch(t){
        case H:return getH(k);
        case C:return getC(k);
        case I:return getI(k);
        case B:return getB(k);
    }
    return -1.0;
}

void Model::set(unsigned k, unsigned t, double v){
    if (k==0 && t==4)
        return;
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

std::ostream& Model::display(std::ostream& os) const{
    os <<_numcomp<<" "<<_numstoch<<std::endl;
    for(unsigned k=0; k < _numcomp; ++k){
        os << k <<" " << getRate(k,0) <<" ";//<< setprecision(6) 
        os << getN(k) << " " << getH(k) <<" " << getC(k) <<" " << getI(k);
        if (k>0) os  <<" "<< getB(k);
        os << std::endl; 
    }
    os << "# " << _numstoch << std::endl;
    os << "# " << _numcomp << std::endl;
    os << "# " << _alpha << std::endl;
    os << "# " << _diagnosis << std::endl;
    return os;
}

std::istream& Model::read(std::istream& is){ //FIXME not working!
    is >> _numcomp;
    is >> _numstoch;
    _endstoch = ((_numstoch-1)*4) + 3;
    for(unsigned i=0; i < _numcomp; ++i){
        unsigned int k;
        double N,H,C,I,B;
        double rate;
        is >> k >> rate ;
        is >> N >>  H >> C >> I;
        if (k > 0) {
            is >> B; 
            setB(k,B);
            storeB(k,B);
        }
        setRate(k,0,rate);
        setH(k,H);
        storeH(k,H);
        setC(k,C);
        storeC(k,C);
        setI(k,I);
        storeI(k,I);
    }
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

void Model::reset_treatment(){

    for(unsigned k=1 ; k < _numcomp; ++k){
        // std::cout <<"before: "<<k<<" "<<getC(k)<<" "<<getB(k)<<std::endl;
        incr(k,C,getB(k));
        setB(k,0);
        // std::cout <<"after: "<<k<<" "<<getC(k)<<" "<<getB(k)<<std::endl;
    }
}


bool Model::updateDet(unsigned k, Data& data){	//FIXME this has to be checked!!!
    assert(k>=_numstoch);



    double epsh = data.epsh(); 
    double epsc = data.epsc();
    double epsr = data.epsr();
    double epsb = data.epsb();
    double prev_epsh = data.epsh(); 
    double prev_epsc = data.epsc();
    double prev_epsr = data.epsr();
    double prev_epsb = data.epsb();
    if (k<=data.n_neutral_compartments()){
        prev_epsc = prev_epsh;
        prev_epsr = prev_epsh;
        prev_epsb = prev_epsh;
    }
    if (k<data.n_neutral_compartments()){
        epsc = epsh;
        epsr = epsh;
        epsb = epsh;
    }
    if(k==_numstoch){
        prev_epsh = 0.0;
        prev_epsc = 0.0;
        prev_epsr = 0.0;
        prev_epsb = 0.0;
    }

    // std::cout << k-1 << " : " << getN(k-1) << " " << getH(k-1) << " " << getC(k-1) << " " << retrieveN(k-1) << " " << retrieveH(k-1) << " " << retrieveC(k-1) << std::endl;
    // std::cout << k << " : " << getN(k) << " " << getH(k) << " " << getC(k) << " " << retrieveN(k) << " " << retrieveH(k) << " " << retrieveC(k) << std::endl;
    // std::cout <<p << " " << data.epsh() << " " << data.epsc() << std::endl;

    double p = data.dt() *getRate(k,0);
    double prev_comp_p = data.dt() *getRate(k-1,0); 

    double tempH= retrieveH(k) + (2.0 * prev_epsh * prev_comp_p * retrieveH(k-1)) 
        + (p * (1.0 - epsh) * retrieveH(k)) 
        - (p * epsh * retrieveH(k));
    //	std::cout << "#current H("<<k<<")=" << getH(k) << std::endl;
    incr(k, H,tempH);
    //	std::cout << "#new H " << getH(k) << std::endl; 

    p = data.dt() *getRate(k,1);
    prev_comp_p = data.dt() *getRate(k-1,1); 
    double tempC= retrieveC(k) + (2.0 * prev_epsc * prev_comp_p * retrieveC(k-1))
        + (p * (1.0 - epsc) * retrieveC(k))
        - (p * epsc * retrieveC(k));
    incr(k, C,tempC);

    p = data.dt() *getRate(k,2);
    prev_comp_p = data.dt() *getRate(k-1,2); 
    double tempI= retrieveI(k) + (2.0 * prev_epsr  * prev_comp_p * retrieveI(k-1)) 
        + (p * (1.0 - epsr) * retrieveI(k))
        - (p * epsr * retrieveI(k));
    // std::cout << "#current I("<<k<<")=" << retrieveI(k) << std::endl; 
    incr(k, I, tempI);
    // std::cout << "#new I("<<k<<")=" << getI(k) << std::endl; 

    p = data.dt() *getRate(k,3);
    prev_comp_p = data.dt() *getRate(k-1,3); 
    double tempB= retrieveB(k) + (2.0 * prev_epsb * prev_comp_p * ((k-1)>0?retrieveB(k-1):0.0))
        + (p * (1.0 - epsb) * retrieveB(k))
        - (p * epsb * retrieveB(k));
    incr(k, B,tempB);
    //	cout << k << " : " << getN(k) << "\t" << getH(k) << "\t" << getC(k) << "\t" << getI(k) << "\t" << getB(k) << std::endl;
    return true;
}

bool Model::treatDeterministically(unsigned k, double amount){	

    double ccells=retrieveC(k);
    double tmp = ccells * amount;
    // if (amount >0.)
    // std::cout <<"#treatment: "<<k<<" "<<tmp<<" "<<ccells<<" "<<amount<<" "<<getC(k);
    incr(k,C,-tmp);
    incr(k,B,tmp);
    // std::cout <<" "<<getC(k)<<std::endl;
    return true;
}

bool Model::containsLSC() const{
    return (getC(0) > 0);
}

void Model::print_cells(std::ostream & os,double _time){
    os <<_time/365.0<<" ";
    os <<getH(0)<<" ";
    os <<getC(0)<<" ";
    os <<getI(0)<<" ";
    os <<std::endl;
}

bool Model::manual_mutation(unsigned int k,unsigned int celltype_from,unsigned int celltype_to){
    if ( get(k,celltype_from)>=1.){
        //do fun stuff
        incr(k,celltype_to,1.);
        decr(k,celltype_from,1.);
        return true;
    }
    else {
        return false;
    }


}


