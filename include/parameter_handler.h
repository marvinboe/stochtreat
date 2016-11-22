/** 
 * File:   ParameterHandler.h
 * Author: Marvin Boettcher
 *
 * Reads an input file and stores all parameters found in the input file for later use
 * To avoid mixup, there should never be a constant defined in Const.h and the parameter file!
 * Arguments given to the programm overide existing parameters
 * Checks if an argument is used during run and gives out warning if not 
 * (Requires to call CheckIfAllUsed at end)
 *
 * The parameter file must be handed over as first argument to the programm
 * If a request for an undefined parameter is found the program terminates
 *
 * Comments are started by #
 **/

#ifndef __PARAMETER_HANDLER__
#define	__PARAMETER_HANDLER__





#include <vector>
#include <string>
// #include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <sstream>
#include <map>



const char* const DEFAULT_PARAMETER_FILE = "parameters.dat";

class ParameterHandler {
public:
    ParameterHandler();
    ParameterHandler(int argc, char** argv);
    ParameterHandler(std::istringstream& instream);
    ~ParameterHandler();

    void HandleInputArguments(int argc, char* argv[]);

    template<typename T> bool SetValue(const char* parameter_name,T &parameter){
        auto it= parameters.find(parameter_name);
        if (it != parameters.end())
        {
            std::stringstream stovar(it->second);
            if (stovar >> parameter)
                ;
            else {
                std::cout <<"writing of  '"<<it->second<<"' to parameter "<<parameter_name<<" failed"<<std::endl;
                exit(1);
            }
        }
        else
        {
            std::cout <<"key '"<<parameter_name<<"' not found in input"<<std::endl;
            parameter=T();
            return false;
        }
        return true;
    }

    template<typename T> void SetValue_werror(const char* parameter_name,T &parameter){
        auto it= parameters.find(parameter_name);
        if (it != parameters.end())
        {
            std::stringstream stovar(it->second);
            if (stovar >> parameter)
                ;
            else {
                std::cout <<"writing of  '"<<it->second<<"' to parameter "<<parameter_name<<" failed"<<std::endl;
                exit(1);
            }
        }
        else
        {
            std::cout <<"key '"<<parameter_name<<"' not found in input"<<std::endl;
            exit(1); 
        }
    }

    template<typename T> void SetValue(const char* parameter_name,const char * helptext,T &parameter){
        auto it= parameters.find(parameter_name);
        keys_and_help[parameter_name]=helptext;
        if (it != parameters.end())
        {
            std::stringstream stovar(it->second);
            if (stovar >> parameter)
                ;
            else {
                std::cout <<"writing of  '"<<it->second<<"' to parameter "<<parameter_name<<" failed"<<std::endl;
                exit(1);
            }
        }
    }

    std::string ReturnString(const char* parameter_name);
    int ReturnInt(const char* parameter_name);
    double ReturnDouble(const char* parameter_name);
    bool ReturnBool(const char* parameter_name);
    void PrintAll();

    void print_help(std::ostream & os);
   
private:

    
    struct s_parameter
    {
       std::string name;
       std::string value;  //stored as char an converted later
       bool used;
    };

    bool debug;
    std::map<std::string,std::string> parameters;

    std::map<std::string,std::string> keys_and_help;

    bool help;
    
    
    bool ParseString(std::string input,s_parameter &new_parameter);

    s_parameter& ParseString(std::string input);

    bool file_loaded;


    

};

#endif	/* ParameterHandler  */

