#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <limits>
#include <boost/lexical_cast.hpp>
#include "../include/getopts.h"
#include "../include/data.h"
#include "../include/rangen.h"
#include "../include/dependency.h"
#include "../include/PriorityQueue.h"
#include "../include/kernel.h"
#include <algorithm>

int main (int argc, char *argv[]) {

    enum {RUNID, SIZE,TREATMENTTIME, MASS, REDUCTION, PATIENTS, PATH,NTIME,OUTSPEC};
    Options options;

    options.addOption("r", "run",	"Run identifier (default 1)",	true);
    options.addOption("s", "size",	"Number of stochastic compartments (default 7)",	true);
    options.addOption("t", "treattime", "Years of treatment (default 10 year)", true);
    options.addOption("m", "mass",	"animal mass (default 70 kg)",	true);
    options.addOption("e", "reduction", "Required reduction level (default 4.5 logs)", true);
    options.addOption("p", "patients", "Number of patients (default 2)", true);
    options.addOption("h", "path", "Output path (default ./) ", true);
    options.addOption("n", "ntime", "Maximum simulation time", true);
    options.addOption("n", "ntime", "Maximum simulation time", true);
    options.addOption("o", "output", "Specifiy kind of output", true);


    options.parse(argc, argv);

    int runid(1);
    int size(7); // 1 means only the stem cell compartment
    float treatmenttime (10);
    float mass(70); //human mass
    float reduction(4.5);
    unsigned patients(2);
    double months(1);
    double factor(1.0);
    string path("./");
    double ntime(-1.0);
    int output_specifier(0);

    int opt;
    while ((opt = options.cycle()) >= 0) {
        switch(opt) {
            case RUNID:
                {
                    runid = boost::lexical_cast<int>(options.getArgs(opt));
                    cout << "#Run identifier is: " << runid << endl;
                    break;
                }
            case SIZE:
                {
                    size = boost::lexical_cast<int>(options.getArgs(opt));
                    cout << "#Number of stochastic compartments is: " << size << endl;
                    break;
                }
            case TREATMENTTIME:
                {
                    treatmenttime = boost::lexical_cast<float>(options.getArgs(opt));
                    cout << "#Treatmenttime set to: " << treatmenttime << endl;
                    break;
                }
            case MASS:
                {
                    mass = boost::lexical_cast<float>(options.getArgs(opt));
                    cout << "#Animal mass is set to: " << mass << endl;
                    break;
                }
            case REDUCTION:
                {
                    reduction = boost::lexical_cast<float>(options.getArgs(opt));
                    cout << "#Requested log reduction is: " << reduction << endl;
                    break;
                }
            case PATIENTS:
                {
                    patients = boost::lexical_cast<unsigned>(options.getArgs(opt));
                    cout << "#Number of patients is: " << patients << endl;
                    break;
                }
            case PATH:
                {
                    path = options.getArgs(opt);
                    cout << "#Path is set to: " << path << endl;
                    break;
                }
            case NTIME:
                {
                    ntime = boost::lexical_cast<float>(options.getArgs(opt));
                    cout << "#set maximum time to: " << ntime << endl;
                    break;
                }
            case OUTSPEC:
                {
                    output_specifier = boost::lexical_cast<int>(options.getArgs(opt));
                    cout << "#set output specifier to: " << output_specifier << endl;
                    break;
                }
            default:
                break;
        }
    }

    double Nbase(16.52861491);
    double Bbase(2.892507609);
    double Sbase(0.034572078);
    double Lbase(8.643019616); //elephant 8.54663017, human 8.643019616

    RanGen ran;
    Data data;
    data.calcFromMass(mass, Nbase, Bbase, Sbase, (Lbase * factor), months);
    if (ntime > 0.){ //non-default
        data.setNtimes(ntime);
    }
    data.setPercBound(0.05);
    data.setStop(12);
    data.setReduction(reduction);
    data.setLimit(size);
    data.setTreatment(treatmenttime);
    //	cout << data << endl;

    double nolsc = 0;
    double diagnosed_nolsc = 0;
    double diagnosis=0;
    double diagnosed=0;
    double reachedreduction = 0;
    double total_timetoreduction = 0;
    double avgsize[size+1];
    for (int i=0; i< size+1; i++) {
        avgsize[i]=0.0;
    }

    //	cout << "#Examine " << patients << " patients" << endl;
    //	cout << data << endl;
    clock_t timer=clock();
    vector<unsigned> redresult;
    for(unsigned i=0 ; i < patients; ++i){
        if (output_specifier==0){
            cout << "#patient " << i << endl;
        }
        Kernel ker(ran, data, size);
        //	    ker.printAll(); cout << endl << endl;
        double time = ker.execute(ran,0.0,false);
        if(!ker.hasLSC())
            nolsc +=1;
        if(ker.reachedDiagnosis()) {
            diagnosed +=1;
            diagnosis += time;
            if(!ker.hasLSC())
                diagnosed_nolsc +=1;
            //start treatment until limit is reached or maxmum time of treatment has passed
            //			cout << "#burden is " << ker.burden() << " reduction is " << ker.getReduction() << endl;
            ker.execute(ran,time,true);
            //			cout << "#burden is " << ker.burden() << " reduction is " << ker.getReduction() << endl;


            if(ker.reachedReduction()){
                reachedreduction +=1;
                double timetoreduction=(ker.whenReduction() - ker.getDiagnosis());
                total_timetoreduction +=timetoreduction;
                redresult.push_back(ker.whenReduction() - ker.getDiagnosis());
                if (output_specifier==0){
                    cout << "#<years to diag.> <years to red.> <total> <nolsctime>"
                    << ker.getDiagnosis() << "  " 
                    << timetoreduction << "  "
                    << ker.whenReduction() << " "
                    << ker.get_nolsctime() << endl;
                }

                // //write compartment data to file
                // stringstream ss;
                // ss << path<< "patient-"<< runid << "-"<< i << ".txt";
                // ofstream output(ss.str().c_str());
                // if(!output.is_open()){
                //     cout << " unable to open output file " << ss.str() << endl;
                //     cout << " exiting program " << endl;
                //     exit(-1);
                // }
                // ker.writeModel(output);
                // output.close();
            }
            //			cout << "Reduction is " << ker.getReduction() << endl;
        }

        if (output_specifier==1){
            cout  << ker.get_nolsctime() << endl;
        }
        ker.addStochCompSizes(avgsize);
    }

    cout << "#Real time elapsed in seconds: " << ((double)clock()-timer)/CLOCKS_PER_SEC << endl;
    cout << "#Fraction diagnosed "<< (diagnosed /(double) patients)*100 << endl;
    cout << "#Average time to diagnosis " << (diagnosed > 0?(diagnosis / (double) diagnosed):0) << endl;
    cout << "#Fraction with no LSC" << (nolsc / (double) patients) << endl;
    cout << "#Fraction diagnosed with no LSC" << (diagnosed > 0?(diagnosed_nolsc / (double) diagnosed):0) << endl;
    cout << "#Fraction that reached "<< reduction 
        << " log reduction : <reduction freq.> <#of reductions> <diagnosed> <noscl at dignose>" << endl <<"# "
        << (reachedreduction > 0?(reachedreduction / (double) diagnosed):0)
        << " " << reachedreduction << " " << diagnosed << " " << diagnosed_nolsc<< endl;
    double stddev = 0;
    double avg = (reachedreduction > 0?(total_timetoreduction / (double) reachedreduction):0) ;
    for(int i=0; i < redresult.size(); i++){
        stddev += pow((double)(redresult[i] - avg), 2.0);
    }
    stddev = stddev / (double)redresult.size();
    stddev = sqrt(stddev);
    cout << "#Average time to reduction "<< avg << "\t" << stddev << endl;
    for (int i=0; i< size+1; i++) {
        cout << "#avg size comp " << i << " = " << (avgsize[i]/(double)patients) << endl;
    }


    return 0;
}


