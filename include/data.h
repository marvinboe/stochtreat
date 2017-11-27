/*
 *  data.h
 *  HemaMass
 *
 *  Created by Tom Lenaerts on 15/08/06.
 *  Copyright 2006 SWITCH. All rights reserved.
 *
 */


#ifndef __DATA_H
#define __DATA_H

#include <string>
#include <istream>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include "parameter_handler.h"

struct Diff_probabilities{
    double epsh=0.85;
    double epsc=0.72;
    double epsb=0.89;
    double epsr=epsc; //differentation probability immmune cell
    double eps_asym=0.0; //asymmetric division probability for stem cells

    /* write differentiation probabilities output os. */
    void write(std::ostream & os){ 
        os <<"#diffprobs: "<<epsh<<" "<<epsc<<" "<<epsb<<" "<<epsr<<" "<<eps_asym<<std::endl;
    }

};

struct Proliferation_parameters{
    double kn=1/365.;
    double kc=kn;
    double kb=kn;
    double kr=kn;

    double gamman=1.263;//in paper 1.26// 1.263;
    double gammac=1.263;
    double gammab=1.263;
    double gammar=gammac;
    void write(std::ostream & os){ 
        os <<"#prolifs (k,gamma): "<<kn<<" "<<kc<<" "<<kr<<" "<<kb<<" "<<gamman<<" "<<gammac<<" "<<gammar<<" "<<gammab<<" "<<std::endl;
    }
};

struct Run_modes{
    int resistance=-1;
    bool treattest=false;
    bool fixed_time_treatment=true;
    operator bool() const { return (resistance>=0||treattest);}
};

struct Simulation_Parameters{
    Diff_probabilities diff_probs;
    Proliferation_parameters prolif;
    Run_modes run_mode;
    std::string output; //"patient nolsctime diagtime initresponse fullburden"
    int runid = 1;
    int n_stochastic_compartments = 7; // 1 means only the stem cell compartment
    int n_neutral = 1; //bcr/abl neutral up to this compartment
    int n_compartments = 32;
    unsigned inital_lsc = 1;
    double diagnosis_level = 12;
    float treatmenttime  = 20;
    float mass = 70; //human mass
    float reduction = 4.5;
    double relapse_logreduction = 3.;
    double required_reduction_time = 0;
    double relapse_waiting_time = 5;
    double treatment_rate = 0.05;
    unsigned patients = 1;
    double collectinterval=30.; //how often data is collected
    double ntime=25.;// maximum simulation time in years
    double dt=0.1;
    void set_parameters(ParameterHandler & inputparams);
};

/** Stores and handles all data for the model and input/ouput.
 * Recieves and handels user input and simulation defaults.
 * Main input:
 * Simulation_Parameters simparams - all user editable simulation parameters
 * Note: When adding new parameters to constructor, also add to copy constructor!!!
 * */
class Data {
    public:
        Data(); 
        Data(const Data&); 
        ~Data(){};

        double dt() const {return _dt;}
        void setdt(double v) {_dt=v;}
        double frac_csc() const {return _frac_csc;}
        void setFrac_csc (double v) {
            _frac_csc = v;
            _numlsc = v*N0();
            double temp;
            modf(_numlsc,&temp); 	
            double diff=_numlsc-temp;

            if(diff >= 0.5)
                _numlsc = (temp+1); // always round up to next complete lsc
            else _numlsc=temp;			
        }
        double numlsc() const {return _numlsc;}
        void setNumLsc(double v) {
            _numlsc = v;
            _frac_csc = v/N0();
        } 
        double N0() const {return _N0;}
        void setN0(double v) {_N0 = v;}

        /* Returns the base proliferation rate for each cell type,
         * which is the stem cell proliferation.*/
        double base_proliferation(unsigned type){
            switch(type){
                case 0:
                    return _prolif.kn;
                case 1:
                    return _prolif.kc;
                case 2:
                    return _prolif.kr;
                case 3:
                    return _prolif.kb;
                default :
                    return -1.0;
            }
        }

        /* Returns the base proliferation rate for each cell type,
         * which is the stem cell proliferation.*/
        double prolif_exp(unsigned type){
            switch(type){
                case 0:
                    return _prolif.gamman;
                case 1:
                    return _prolif.gammac;
                case 2:
                    return _prolif.gammar;
                case 3:
                    return _prolif.gammab;
                default :
                    return -1.0;
            }
        }

        Proliferation_parameters return_prolif_params() const{ return _prolif;}

        double mass() const {return _mass;}
        void setMass(double v) {_mass = v;}

        double epsh() const {return _diffprobs.epsh;}
        double epsc() const {return _diffprobs.epsc;}
        double epsb() const {return _diffprobs.epsb;}
        double epsr() const {return _diffprobs.epsr;}
        double eps_asym() const {return _diffprobs.eps_asym;}

        /** return self-renewal probability epsilon for type.
         * 0: healthy, 1: cancerous, 2: resistant, 3: bound (treated).
         * Special case 4: asymmetric division probability for stem cells.*/
        double eps(unsigned type) const {
            switch(type){
                case 0:
                    return _diffprobs.epsh;
                case 1:
                    return _diffprobs.epsc;
                case 2:
                    return _diffprobs.epsr;
                case 3:
                    return _diffprobs.epsb;
                case 4:
                    return _diffprobs.eps_asym;
                default :
                    return -1.0;
            }
        }

        void setEpsh (double v) {_diffprobs.epsh = v;} 
        void setEpsc (double v) {_diffprobs.epsc = v;} 
        void setEpsb (double v) {_diffprobs.epsb = v;} 
        void setEpsi (double v) {_diffprobs.epsr = v;} 
        void set_eps_asym (double v) {_diffprobs.eps_asym = v;} 
        void setEps(unsigned type, double v)  {
            switch(type){
                case 0:
                    _diffprobs.epsh = v;
                    break;
                case 1:
                    _diffprobs.epsc = v;
                    break;
                case 2:
                    _diffprobs.epsr = v;
                    break;
                case 3:
                    _diffprobs.epsb = v;
                    break;
                case 4:
                    _diffprobs.eps_asym = v;
                    break;
            }
        }

        double p_csc() const {return _p_csc;}
        void setPcsc (double v) {_p_csc = v;} 
        double p_imm() const {return _p_imm;}
        void setPimm (double v) {_p_imm = v;} 


        /** returns the treatment rate of all cancer cells. */
        double treatment_rate() const {return _treatment_rate;}

        /** sets the percentage of cells that is affected by treatment. */
        void set_treatment_rate(double v) {_treatment_rate = v;}


        /** returns the total number of compartments in the model. */
        int ncompartments() const {return _ncompartments;}

        /** returns the number of bcrabl-neutral compartments. */
        unsigned int n_neutral_compartments() const {return _n_neutral_compartments;}

        /** Sets the total number of compartments in the model. */
        void setNCompartments(int v) {_ncompartments = v;}

        /** Returns the threshold power of 10 when diagnosis is reached:
         * When 10^(stop) cells are in the compartment -> diagnosis. */
        double diagnosis_level() const {return _diagnosis_level;}

        /** Sets the threshold power of 10 when diagnosis is reached:
         * When 10^(stop) cells are in the compartment -> diagnosis. */
        void set_diagnosis_limit(double v) {_diagnosis_level = v;}

        /** Returns the recuction level that is required to stop treatment.*/
        double reduction() const {return _reduction;}
        void set_treatment_stop_reduction(double v) {_reduction = v;}

        /** Returns required time in days for reduction to be 
         * maintained before treatment is stopped.*/
        double required_reduction_time() const {return _required_redtime;}
        void set_required_reduction_time(double v) {_required_redtime=v;}

        void set_relapse_reduction(double v) {_relapse_reduction=v;}
        double relapse_reduction() {return _relapse_reduction;}

        double additional() const {return _additional;}
        void setAdditional(double a)  {_additional=a;}

        /** returns the number of stochstic compartments */
        int nstochcomp() const {return _numstochcomps;}

        /** Sets the number of stochastic compartments.*/
        void set_numstochcomps(int l) {_numstochcomps = l;}

        /** threshold of cellnumber for diagnosis. */
        double threshold() const {return _threshold;}
        void setThreshold(double v) {_threshold = v;}


        int step() const {return _outputstep;}
        void setStep(int v) { _outputstep = v;}

        /** returns treatment time in years.*/
        double treatment_dur() const {return _treatment_duration;}

        /** Returns the maximum time a simulation runs in years.*/
        double getTmax_in_years() const {return _tmax;}

        /** returns waiting time for cancer relapse.*/
        double relapse_waiting_time() const {return _relapse_waiting_time;}

        /** Sets maximum simulation time. */
        void setTmax(double v){_tmax=v;}

        /** Sets treatment time in years that will be used 
         * in the treatment phase of the simulation. */
        void set_maximum_treatment_duration(double t) {_treatment_duration=t;}

        /** Calculates patient parameters from given input.
         * parameters:
         * mass     - mass of the modeled animal
         * Nbase    - log base for the number of hematopeotic stem cells
         * Bbase    - log base for the average cell cycle time of hematopeotic stem cells
         * Sbase    - log base for the deterministic timestep of simulation
         * Lbase    - log base for maximum simulation time
         * c_interv - interval for virtual doctor visits (data collection interval)
         * diffprobs- differentiation probabilies. */
        void initialize(double,double,double, double, double, double,Diff_probabilities);

        /** Calculates patient parameters from given input.
         * parameters:
         * simparams- simulation parameters (from default and user input).
         * Nbase    - log base for the number of hematopeotic stem cells
         * Bbase    - log base for the average cell cycle time of hematopeotic stem cells
         * Sbase    - log base for the deterministic timestep of simulation
         * Lbase    - log base for maximum simulation time */
        void initialize(const Simulation_Parameters & ,double, double,double);

        friend std::ostream & operator<<(std::ostream &o, Data& c){return c.display(o);}

    private:
        std::ostream& display(std::ostream&);
        std::istream& read_from_file(std::istream&);

        double _dt; //time step relative to days
        Diff_probabilities _diffprobs;//differentiation probabilities of all cell types
        Proliferation_parameters _prolif; //proliferation parameters of all cell types
        double _p_csc; //probability that a normal cell turns into a cancer cell
        double _p_imm; //probabilty that a cancer cell turns into an immune cell
        double _frac_csc; //fraction of cancer cells in the stem cell compartment
        double _numlsc; //number of LSC
        double _treatment_rate; //percentage of cells bound to imatinib per day
        double _tmax; // maximum simulation time in years
        int _ncompartments;  //number of compartmens in the hematopoeitic system
        int _n_neutral_compartments; //number of compartments where bcr/abl is neutral
        double _N0; //Numbe of cells in the stem cell compartment
        int _numstochcomps; //index of first deterministic compartment
        double _additional; //additional number of years to continue simulation after X
        double _treatment_duration; //number of years of treatment
        double _diagnosis_level; //stop value = diagnosis level
        double _reduction; //stop value (required log reduction in bcr-abl transcript level)
        double _required_redtime; //time stop value to be maintained
        double _relapse_reduction; //stop value relapse
        double _relapse_waiting_time; //time to wait for relapse after stopping of treatment
        double _mass;  //mammal mass
        double _threshold; //percentage increase in number of cells for diagnosis

        int _outputstep;//steps after which output is saved. unused
};

#endif
