/*
 *  patientdata.h
 *  stochestimate
 *
 *  Created by Tom Lenaerts on 12/07/11.
 *  Copyright 2011 Université Libre de Bruxelles. All rights reserved.
 *
 */

#ifndef __PATIENTDATA_H
#define __PATIENTDATA_H

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "model.h"



class PatientData {
public:
	PatientData(std::string str, int p);
	PatientData(const PatientData&);
	
	unsigned sizeData() const {return (unsigned)_data.size();} // should be the same as the timepoints size
	unsigned sizeEstimate() const {return (unsigned)_estimate.size();}
	unsigned sizefQuartile() const {return (unsigned)_fquartile.size();}
	unsigned sizetQuartile() const {return (unsigned)_tquartile.size();}
	unsigned sizeMin() const {return (unsigned)_min.size();}
	unsigned sizeMax() const {return (unsigned)_max.size();}
	unsigned sizeMedian() const {return (unsigned)_median.size();}
	unsigned sizeMonthly() const {return (unsigned)_monthly.size();}
	
	float time(unsigned pos) const {return _timepoints[pos];}
	int posTime(float t);
	float data(unsigned pos) const {return _data[pos];}
	float fqvalue(unsigned pos) const {return _fquartile[pos];}
	float tqvalue(unsigned pos) const {return _tquartile[pos];}
	float minvalue(unsigned pos) const {return _min[pos];}
	float maxvalue(unsigned pos) const {return _max[pos];}
	float median(unsigned pos) const {return _median[pos];} 
	float monthly(unsigned pos) const {return _monthly[pos];} 

	bool clearEstimates() {_estimate.clear(); _lasttime=0; return (_estimate.size() == 0);}
	bool clearMonthly() {_monthly.clear();  return (_monthly.size() == 0);}
	
	bool estimateTime(double time);
	void recordEstimate(Model& _pool);
	bool monthlyTime(double time);
	void recordMonthly(Model& _pool);
	float lastMonth() const {return _lastmonth;}
	
	float estimate(unsigned pos) const {return _estimate[pos];}
	void addEstimate(float e) {_data.push_back(e);}
	double estimateError();
	
	unsigned patientID() const { return _patient;}
	
	PatientData& operator=(const PatientData& other);
	
	friend std::ostream & operator<<(std::ostream &o, PatientData& pd){return pd.display(o);}
	
private:
	std::ostream& display(std::ostream&);
	void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters="\t ");
	void processTimepoints(std::vector<std::string>& tokens, int p);
	void processPatientdata(std::vector<std::string>& tokens, int p);
	void processData(std::vector<std::string>& tokens, int p);
	void transformAllData(std::vector< std::vector<float> >& transformed);
	void calculateStatistics();
	void parseBuffer(std::stringstream& ss);
	void readBinaryFile(std::istream& is);
	
	unsigned _patient;
	std::vector<float> _timepoints;
	std::vector<float> _data;
	std::vector<float> _estimate;
	std::vector<float> _monthly;
	
	std::vector<float> _fquartile;
	std::vector<float> _tquartile;
	std::vector<float> _min;
	std::vector<float> _max;
	std::vector<float> _median;
	unsigned _lasttime;
	float _lastmonth;
	//temporary variable
	std::vector< std::vector<float> > _alldata;
	
};

#endif
