/*
 *  patientdata.cpp
 *  stochestimate
 *
 *  Created by Tom Lenaerts on 12/07/11.
 *  Copyright 2011 Université Libre de Bruxelles. All rights reserved.
 *
 */

#include "patientdata.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <boost/lexical_cast.hpp> 
#include <cmath>
#include <cassert>

int PatientData::posTime(float t){
    unsigned i(0);
    for(; i < _timepoints.size() && _timepoints[i] != t; ++i);
    return i;
}

void PatientData::tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters){
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the std::vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

void PatientData::processTimepoints(std::vector<std::string>& tokens, int p){
    for (unsigned i=1; i<tokens.size()-1; ++i) {
        _timepoints.push_back(boost::lexical_cast<float>(tokens[i]));
    }
}

void PatientData::processPatientdata(std::vector<std::string>& tokens, int p){
    for (unsigned i=1; i<tokens.size()-1; ++i) {
        _data.push_back(boost::lexical_cast<float>(tokens[i]));
    }
}

void PatientData::processData(std::vector<std::string>& tokens, int p){
    std::vector<float> tmp;
    for (unsigned i=1; i<tokens.size()-1; ++i) {
        float value = boost::lexical_cast<float>(tokens[i]);
        tmp.push_back(value);
    }	
    _alldata.push_back(tmp);
}

void PatientData::parseBuffer(std::stringstream& ss){
    std::string line;
    int i=0; //line in file = column in matrix datapoints
    bool patient_found = false;
    while(getline(ss,line,'\n')){
        std::vector<std::string> tokens;
        tokenize(line,tokens,"\t");
        if(tokens.size() > 0){
            unsigned which = boost::lexical_cast<unsigned>(tokens[0]);
            if(i == 0) { //time points
                processTimepoints(tokens,i);			
            }
            if(i == _patient && _patient > 0 && !patient_found){
                processPatientdata(tokens,i);
                patient_found = true;
                _patient = which;
            }
            if(i > 0){
                processData(tokens,i);
            }
        }
        ++i;
    }
    if(!patient_found && _patient >0){
        std::cout << "!! Requested patient " << _patient << " NOT FOUND !! " << std::endl;
        exit(-1);
    }
}


void PatientData::readBinaryFile(std::istream& is){
    is.seekg (0, std::ios::end);
    int length = (int)is.tellg();
    is.seekg (0, std::ios::beg);
    _alldata.clear();
    // allocate memory:
    char* buffer = new char [length];
    // read data as a block:
    is.read (buffer,length);
    //parse data
    std::stringstream ss;
    ss << buffer;
    parseBuffer(ss);

    delete[] buffer;
}

void PatientData::transformAllData(std::vector< std::vector<float> >& transformed){
    unsigned numpatients = (unsigned)_alldata.size(); // number of patients
    unsigned numpoints = (unsigned)_alldata[0].size(); // number of datapoints per patietn

    for(int column=0; column <  numpoints; column++){
        std::vector<float> tmp;
        for(int row=0; row <  numpatients; row++){
            if(_alldata[row][column] >= 0)
                tmp.push_back(_alldata[row][column]);
        }
        transformed.push_back(tmp);
    }
}

void PatientData::calculateStatistics(){
    _min.clear();
    _max.clear();
    _fquartile.clear();
    _tquartile.clear();		
    _median.clear();
    std::vector< std::vector<float> > tdata;
    transformAllData(tdata);
    unsigned numpoints = (unsigned)tdata.size(); // number of datapoints per patient

    for(int column=0; column <  numpoints; column++){
        //calculate statistics for every column
        //median, first quartile and thrid quartile
        std::vector< float> tmp(tdata[column]);
        sort(tmp.begin(), tmp.end());

        //median
        unsigned median_pos = (unsigned)tmp.size()/2.0;
        _median.push_back(tmp[median_pos]);

        //first quartile
        float firstQ_tmp = (median_pos)/2.0;
        float intpart;
        std::modf(firstQ_tmp,&intpart); 	
        unsigned firstQ_pos = (median_pos)/2.0;
        if( (firstQ_tmp-intpart) >0.0 ){
            //firstQ is average of 2 positions
            _fquartile.push_back( (tmp[firstQ_pos]+tmp[firstQ_pos+1])/2.0 );
        }
        else {
            _fquartile.push_back(tmp[firstQ_pos]);
        }

        //third quartile
        float thirdQ_tmp = median_pos + ((tmp.size()- (median_pos + 1))/2.0);
        std::modf(thirdQ_tmp,&intpart); 	
        unsigned thirdQ_pos = median_pos + ((tmp.size()- (median_pos + 1))/2.0);
        if( (thirdQ_tmp-intpart) > 0.0 ){
            _tquartile.push_back( (tmp[thirdQ_pos] + tmp[thirdQ_pos+1])/2.0 );
        }
        else {
            _tquartile.push_back(tmp[thirdQ_pos]);
        }

        _min.push_back(tmp[0]);
        _max.push_back(tmp[tmp.size()-1]);

    }
    if(_patient == 0){
        _data = _median;
    }
}


PatientData::PatientData(std::string fname, int p):_patient(p),_lasttime(0), _lastmonth(0.0){
    std::ifstream patientstream(fname.c_str(), std::ios::binary);
    if( patientstream.is_open()){
        readBinaryFile(patientstream);
        _estimate.clear();
        _monthly.clear();
        calculateStatistics();
    }
}

PatientData::PatientData(const PatientData& other){
    if(other.sizeEstimate() == other.sizeData()	)	
        _lasttime = -1;
    else _lasttime = 0;
    _lastmonth = other.lastMonth();
    _patient = other.patientID();
    _timepoints.clear();
    _data.clear();
    _fquartile.clear();
    _tquartile.clear();
    _min.clear();
    _max.clear();		
    _estimate.clear();		
    _median.clear();
    _monthly.clear();

    for(unsigned i = 0; i < other.sizeData(); ++i){
        _timepoints.push_back(other.time(i));
        _data.push_back(other.data(i));
        _fquartile.push_back(other.fqvalue(i));
        _tquartile.push_back(other.tqvalue(i));
        _min.push_back(other.minvalue(i));
        _max.push_back(other.maxvalue(i));		
        _median.push_back(other.median(i));		
    }
    for(unsigned i = 0; i < other.sizeEstimate(); ++i){
        _estimate.push_back(other.estimate(i));
    }
    for(unsigned i = 0; i < other.sizeMonthly(); ++i){
        _monthly.push_back(other.monthly(i));
    }	
}

PatientData& PatientData::operator=(const PatientData& other){
    if(other.sizeEstimate() == other.sizeData()	)	
        _lasttime = -1;
    else _lasttime = 0;
    _lastmonth = other.lastMonth();
    _patient = other.patientID();
    _timepoints.clear();
    _data.clear();
    _fquartile.clear();
    _tquartile.clear();
    _min.clear();
    _max.clear();		
    _estimate.clear();		
    _median.clear();
    _monthly.clear();

    for(unsigned i = 0; i < other.sizeData(); ++i){
        _timepoints.push_back(other.time(i));
        _data.push_back(other.data(i));
        _fquartile.push_back(other.fqvalue(i));
        _tquartile.push_back(other.tqvalue(i));
        _min.push_back(other.minvalue(i));
        _max.push_back(other.maxvalue(i));		
        _median.push_back(other.median(i));		
    }

    for(unsigned i = 0; i < other.sizeEstimate(); ++i){
        _estimate.push_back(other.estimate(i));
    }		
    for(unsigned i = 0; i < other.sizeMonthly(); ++i){
        _monthly.push_back(other.monthly(i));
    }	
    return *this;
}


std::ostream& PatientData::display(std::ostream& os){
    os << "[" << _patient << ":";
    for(unsigned i = 0; i < _data.size(); ++i){
        os << "(" << _timepoints[i] << ";" << _data[i] << ")";
        if(i < (_data.size() -1))
            os << ", "; 
    }
    os << "],[ESTIMATE:";
    for(unsigned i = 0; i < _estimate.size(); ++i){
        os << "(" << _timepoints[i] << ";" << _estimate[i] << ")";
        if(i < (_estimate.size() -1))
            os << ", "; 
    }
    os << "],[ STATS:";
    for(unsigned i = 0; i < _fquartile.size(); ++i){
        os << "(" << _fquartile[i] << ";"<< _median[i] << ";" << _tquartile[i]<< ";" << _min[i]<< ";" << _max[i] << ")";
        if(i < (_data.size() -1))
            os << ", "; 
    }
    os << "]";
    return os;
}

bool PatientData::estimateTime(double time){
    bool needtorecord = false;
    if(_lasttime != -1 && time >= _timepoints[_lasttime]) { //_timepoints expressed in years !
        needtorecord = true;
        _lasttime = (_lasttime < (_timepoints.size()-1)? _lasttime + 1: -1);
    }
    return needtorecord;	
}

void PatientData::recordEstimate(Model& _pool){
    //store disease burden in estimate folder.
    _estimate.push_back(_pool.diseaseBurden());
    //	cout << _data[_estimate.size()-1] << "\t" << _pool.diseaseBurden() << "\t" << _estimate.size() << std::endl;;
}

bool PatientData::monthlyTime(double time){
    bool needtorecord = false;
    if(time >= _lastmonth) { 
        needtorecord = true;
        _lastmonth += ((365.0/12.0)/365.0) ;
    }
    return needtorecord;	
}

void PatientData::recordMonthly(Model& _pool){
    //store disease burden in _monthly std::vector
    _monthly.push_back(_pool.diseaseBurden());
}


double PatientData::estimateError(){
    double err(0.0);
    //	cout << _estimate.size() << " \t " <<  _timepoints.size() << std::endl;
    assert(_estimate.size() == _timepoints.size());
    //	cout << "[ESTIMATE : " << std::endl;
    for(unsigned i = 1; i < _estimate.size(); ++i){ //start from one size values in position zero are both 100.
        //		cout << _data[i] << "\t" << _estimate[i] << "\t" << _tquartile[i] << "\t" << _fquartile[i] << std::endl;
        //		cout << _estimate[i] << " \t";
        if(_data[i] != -1)
            err += (pow(_data[i] - _estimate[i],2) / pow(_tquartile[i] - _fquartile[i],2));
    }
    //	cout << std::endl;
    return err;
}

