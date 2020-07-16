#ifndef OUTPUT_H
#define OUTPUT_H
//#include "observer.h"
#include "vector3d.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <vector>
#include "cell.h"
//#include "events.h"

#ifdef _WIN32
#include <windows.h>
#endif

#ifdef __linux__
#include <sys/stat.h>
#endif

#ifdef __APPLE__
#include <sys/param.h>
#include <sys/stat.h>
#endif


using namespace std;



class simulation;
class parameters;
class cell;

class output
{
public:
    string parfname;
    output(string _parfname = string());
    ~output();
    void record_output_time_step(double currentTime, simulation& currentSim, parameters &p);
    void recordEvent(double currentTime  , double value);
//    short Output_ID;
    void initialize_fileds();
    void clear_fileds();
    void write2file_time(double currentTime, simulation& currentSim, parameters&p);
    void write_event( cell* Cellx, stringstream &sim_output);
    void close_event(B_cell* Cellx, stringstream &sim_output,double time);

    void write_event_2file(stringstream &sim_output);

    void Plasma_output(double currentTime, simulation& currentSim, parameters&p);
    void Memory_output(double currentTime, simulation& currentSim, parameters&p);
//    void createFolder(string folderName,parameters &p,short simID);
//    string output_path;
    string currentDateTime();
    vector<string> storage;

};


struct ministat
{
    
    ministat()
    {
        clear_ministat();
    }
    void clear_ministat()
    {
        N=0;
        sum=0;
        sumsq=0;
        values.clear();
    }
    void add(double v)
    {
        if(not(isnan(v))){
            sum += v;
            sumsq += v*v;
            N++;}
        else {
            cout<<"NAN value in output!";
        }
        values.push_back(v);
    }
    int N;
    vector<double> values;
    double sum;
    double sumsq;
    double average()
    {
        // if(N==NAN){return 0;}
        if (N<=0){return 0;}
        if(isnan(sum)){return 0;} //ELENA
        else{return double(sum / (double) N);}
    }
    double stddev()
    {

        if (N<=1){return 0;}
        if(sum == 0){return 0;} //ELENA
        if (isnan(sum)) {return 0;}
        else
        {   if (N>1)
        {
            double std = (double) sqrt( (double) (sumsq / (double) N) - (sum / (double) N) * (sum / (double) double (N - 1)));
            if (isnan(std))
            {return 0.0;}
            return std;
        }
        else {
            return 0;
        }
            
        }
    }
    string print()
    {
        stringstream res;
        res << "[" << N << "]";
        for(int i = 0; i < N; ++i)
        {
            res << "\t" << values[i];
        }
        return res.str();
    }
};




//ministats
extern vector<ministat> Bcell_counts;   //number of Bcells
extern ministat Plasma_counts;

//Mutations
extern ministat Bcell_mutation;
extern ministat Plasma_mutations;

//Affinities
extern ministat Bcell_affinity;
extern ministat Plasma_affinity;

//Antigens
extern vector<ministat> Bcell_Antigen;

//Antibodies
extern ministat Antibody;

//CentroBlasts cycle and divisions
extern ministat Bcell_cycle;
extern ministat Bcell_divisions;
extern ministat Plasma_antigen;
extern ministat Plasma_divisions;
#endif // OUTPUT_H

