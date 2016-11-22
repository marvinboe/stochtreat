//
/*
 * File:   ParameterHandler.cpp
 */


#include "parameter_handler.h"



ParameterHandler::ParameterHandler()
{
    debug=false;
    std::cout <<"no input arguments passed: Exit(1)"<<std::endl;
    exit(1);
}

ParameterHandler::ParameterHandler(std::istringstream& DataFile)
{
    debug=false;
    if (debug)
        std::cout<<"----------parsing parameter file and arguments--------------"<<std::endl;
   

    std::string line_in_file;

    s_parameter new_parameter;
    while (std::getline(DataFile,line_in_file))      //until end of file is reached
    {
//         std::cout <<line_in_file<<std::endl;
        if (ParseString(line_in_file, new_parameter))
        {
//             std::cout <<" " << new_parameter.name<<" "<<new_parameter.value<<" "<<parameters.size()<<std::endl;
            parameters.at(new_parameter.name)=new_parameter.value;
            
        }


    }
    
    if (debug)
        std::cout<<"----------parsing finished----------------------------------"<<std::endl;
}


ParameterHandler::ParameterHandler(int argc, char** argv):file_loaded(false)
{

    std::ifstream DataFile;
    std::string file_name;
    help=false;
    
    if (argc>=2) // second argument should be the input file
    {
        DataFile.open(argv[1]);     //then try to open file with the corresponding name
        if (DataFile.fail())
        {
            //second argument might be command line input
            // std::cout<<"Could not open parameter file '"<<argv[1]<<"'"<<std::endl;
            // exit(1);
        }
        else {
            file_loaded=true;
            file_name=argv[1];
        }
    }
    else //open default file... has to be provided
    {
        return;
        DataFile.open(DEFAULT_PARAMETER_FILE);
        if (DataFile.fail())
        {
            std::cout<<"Could not open default parameter file '"<<DEFAULT_PARAMETER_FILE<<"'"<<std::endl;
            exit(1);
        }
        else
        {
            std::cout<<"Using default parameter file"<<std::endl;
            file_name=DEFAULT_PARAMETER_FILE;
        }
    }
    
    std::string line_in_file;//buffer for line based reading
    s_parameter new_parameter;
    
    while (getline(DataFile,line_in_file))      //until end of file is reached
    {
        if (ParseString(line_in_file, new_parameter))
        {
            parameters[new_parameter.name]=new_parameter.value;
        }
        
    }

    HandleInputArguments(argc,argv);
    
    DataFile.close();
}





ParameterHandler::~ParameterHandler()
{
}


bool ParameterHandler::ParseString(std::string input, s_parameter &new_parameter)
{

    //check if comment or empty
    if (input[0]=='#'  || input.empty()||(input.find_first_not_of(" \t")==std::string::npos))
    {
        //nothing
        return false;
    }

    //truncating string after # for comments
    size_t found_at=input.find("#");
    if (found_at!=std::string::npos)
    {
        input=input.substr(0,found_at);
    }

    found_at=input.find("help");
    if (found_at!=std::string::npos){
        help=true;
        return false;
    }
    found_at=input.find("-h");
    if (found_at!=std::string::npos){
        help=true;
        return false;
    }



    found_at=input.find("=");
    if (found_at==input.length()-1)
    {
        std::cout<<"parameter '"<<input<<"' has no value"<<std::endl;
        return false;
    }
    else
    {
        new_parameter.name=input.substr(0,found_at);
        new_parameter.value=input.substr(found_at+1,std::string::npos);



        //remove trailing and leading whitespaces or tabs from both strings
        size_t found_at_end;
        found_at=new_parameter.name.find_first_not_of(" \t");
        found_at_end=new_parameter.name.find_last_not_of(" \t");
        new_parameter.name=new_parameter.name.substr(found_at,found_at_end-found_at+1);

        found_at=new_parameter.value.find_first_not_of(" \t");
        found_at_end=new_parameter.value.find_last_not_of(" \t");
        new_parameter.value=new_parameter.value.substr(found_at,found_at_end-found_at+1);


        new_parameter.used=false;
        //             std::cout <<new_parameter.name<<" "<<new_parameter.value<<std::endl;
        return true;
    }

}


void ParameterHandler::HandleInputArguments(int argc, char* argv[])
{
    int checkno=(file_loaded?2:1);
    if (argc<=checkno)
    {
        //nothing to be done
        return;
    }
    else
    {
        for (int i=checkno;i<argc;i++)
        {
            s_parameter new_parameter;
            //string argument=argv[i];
            if (ParseString(argv[i], new_parameter))
            {

                std::map<std::string,std::string>::iterator it;

                if ((it= parameters.find(new_parameter.name))!=parameters.end()) {
                    if (debug) 
                    {
                        std::cout<<"override "<<it->first<<"=";//<<(*it).second<<endl;
                        std::cout<<"by       "<<new_parameter.name<<"="<<new_parameter.value<<std::endl;
                    }
                    (*it).second=new_parameter.value;
                }
                else
                {
                    parameters.insert(std::pair<std::string,std::string>(new_parameter.name,new_parameter.value));
                }


            }
        }
    }
}




void ParameterHandler::PrintAll()
{
    for (std::map<std::string,std::string>::iterator it=parameters.begin();it != parameters.end() ;it++)
    {
        std::cout<<"name="<<std::setw(20)<<it->first<<", value="<<it->second<<std::endl;
    }

}



template<> bool ParameterHandler::SetValue<bool>(const char* name, bool& b_parameter)
{
    std::map<std::string,std::string>::iterator it= parameters.find(name);
    if (it != parameters.end())
    {
        if (     (*it).second.find("true")!=std::string::npos
                 || atoi((*it).second.c_str())==1)
        {
            b_parameter=true;
            //             v_parameters[index].used=true;
            return true;
        }
        else if ((*it).second.find("false")!=std::string::npos
                 || atoi((*it).second.c_str())==0)
        {
            b_parameter=false;
            //             v_parameters[index].used=true;
            return true;
        }
        else
        {
            std::cout<<"Cannot set boolean parameter '"<<(*it).first<<"'"<<std::endl;
            std::cout<<"the value '"<<(*it).second<<"' is not boolean"<<std::endl;
            std::cout<<"'"<<(*it).first<<"' is set to false"<<std::endl;
            b_parameter=false;
            //             v_parameters[index].used=true;
            return false;
        }


    }

    return false;
}


template <> void ParameterHandler::SetValue_werror<bool>(const char* name, bool& b_parameter)
{
    std::map<std::string,std::string>::iterator it= parameters.find(name);
    if (it != parameters.end())
    {
        if (     (*it).second.find("true")!=std::string::npos
                 || atoi((*it).second.c_str())==1)
        {
            b_parameter=true;
            //             v_parameters[index].used=true;
        }
        else if ((*it).second.find("false")!=std::string::npos
                 || atoi((*it).second.c_str())==0)
        {
            b_parameter=false;
            //             v_parameters[index].used=true;
        }
        else
        {
            std::cout<<"Cannot set boolean parameter '"<<(*it).first<<"'"<<std::endl;
            std::cout<<"the value '"<<(*it).second<<"' is not boolean (1,0,true,false)"<<std::endl;
            exit(1); 
        }


    }
    else
    {
            std::cout<<"Cannot set boolean parameter '"<<name<<"'"<<std::endl;
	    std::cout <<"key '"<<name<<"' not found in input"<<std::endl;
            exit(1); 
    }

}

template<> void ParameterHandler::SetValue<bool>(const char* name,const char* helptext, bool& b_parameter)
{
    std::map<std::string,std::string>::iterator it= parameters.find(name);
    keys_and_help[name]=helptext;

    if (it != parameters.end())
    {
        if (     (*it).second.find("true")!=std::string::npos
                 || atoi((*it).second.c_str())==1)
        {
            b_parameter=true;
            //             v_parameters[index].used=true;
            return;
        }
        else if ((*it).second.find("false")!=std::string::npos
                 || atoi((*it).second.c_str())==0)
        {
            b_parameter=false;
            //             v_parameters[index].used=true;
            return;
        }
        // else
        // {
        //     std::cout<<"Cannot set boolean parameter '"<<(*it).first<<"'"<<std::endl;
        //     std::cout<<"the value '"<<(*it).second<<"' is not boolean"<<std::endl;
        //     std::cout<<"'"<<(*it).first<<"' is set to false"<<std::endl;
        //     b_parameter=false;
        //     //             v_parameters[index].used=true;
        //     return false;
        // }


    }

    return;
}

void ParameterHandler::print_help(std::ostream & os){

    if (help){
        os <<" usage: './executable argument1=x argument2=y ...'"<<std::endl;
        os <<" or './executable parameterfile argument1=x argument2=y ...'"<<std::endl;
        os << "<argument name> \t description"<<std::endl;
        os << "-----------------------------------"<<std::endl;
        for (const auto keyhelp_pair: keys_and_help){
            os <<'<'<<keyhelp_pair.first<<'>'<<"\t\t"<<keyhelp_pair.second<<std::endl;
        }
        exit(0);
    }

}



int ParameterHandler::ReturnInt(const char* name)
{
    int to_be_returned=0;
    SetValue(name,to_be_returned);
    return to_be_returned;
}


double ParameterHandler::ReturnDouble(const char* name)
{
    double to_be_returned=0.;
    SetValue(name,to_be_returned);
    return to_be_returned;
}




bool ParameterHandler::ReturnBool(const char* name)
{
    bool to_be_returned=false;
    SetValue(name,to_be_returned);
    return to_be_returned;
}




std::string ParameterHandler::ReturnString(const char* name)
{
    std::string to_be_returned;
    SetValue(name,to_be_returned);
    return to_be_returned;
}



// void ParameterHandler::CheckIfAllUsed()
// {
//     std::cout<<"-------------------------------------"<<endl;
//     std::cout<<"The following parameters are not used"<<endl;
//     int count=0;
//     for (int i=0;i<v_parameters.size();i++)
//     {
//         if (v_parameters[i].used==false)
//         {
//             std::cout<<"name="<<setw(20)<<v_parameters[i].name<<", value="<<v_parameters[i].value<<endl;
//             //printf("name=%20s \t value=%10s \n",v_parameters[i].name.c_str(),v_parameters[i].value.c_str());
//             count++;
//
//         }
//     }
//     if (count==0)
//     {
//         std::cout<<"all parameters are used"<<endl;
//     }
//     std::cout<<"-------------------------------------"<<endl;
// }
// kate: indent-mode cstyle; space-indent on; indent-width 0;
