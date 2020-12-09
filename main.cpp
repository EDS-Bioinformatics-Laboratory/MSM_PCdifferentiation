#include "mafalda.h"
#include "random.h"
//#include "GC3D.h"
#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>
//#include <sys/stat.h>
using namespace std;
string outputFolder=string();
#include <time.h>
#ifdef _WIN32
#include <windows.h>
#endif
#ifdef __linux__
#include <sys/stat.h>
#endif
#ifdef __APPLE__
#include <sys/param.h>
//#include <boost/filesystem.hpp>
#endif


void createFolder(string folderName) {


#ifdef _WIN32
  const char *p = tmp.str();
  const WCHAR *pwcsName;
  int nChars = MultiByteToWideChar(CP_ACP, 0, p, -1, NULL, 0);
  pwcsName = new WCHAR[nChars];
  MultiByteToWideChar(CP_ACP, 0, p, -1, (LPWSTR)pwcsName, nChars);
  CreateDirectory(pwcsName, NULL);
  delete[] pwcsName;
#endif

//#Recheck
#ifdef __linux__
  const int dir_err = mkdir(tmp.str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if ((-1 == dir_err) && (errno != EEXIST)) {
    cerr << "Error creating directory : " << tmp << endl;
  }
#endif

#ifdef __APPLE__
      char cstr [folderName.size()+1];
      strcpy(cstr,folderName.c_str());
      mkdir(cstr, 0777);
#endif
}



int main(int argc, char** argv){

    cout << "Starting Mafalda " << endl;
    cout << "------------------------------------- How to use: --------------------------------------------------------\n";
    cout << "Mafalda hyphasmaparameterfile.par -h\n";
    cout << "  ... additional options that can be used:\n";
    cout << "  -s seedNumber \n";
    cout << "  -o outputFolder    or -o auto    to create and output in a folder with the parameter file name\n";
    cout << "----------------------------------------------------------------------------------------------------------\n";


//    srand(1233244567887); //Elena: If you want to fix seed from inside code do inside "no argument detected option!"
//    srand(time(NULL));


    ///////////////////////////////////// Parsing arguments given from command line ///////////////////////////////7

    // See https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html#Getopt-Long-Option-Example

    string parfname = string();
    //Elena: we can define here another parameter file for extensions. How can we do for different files. Eg: network, kon/of, others...
    string parfname2 = string();

    int takeHyphasmaFile = false; // Using hyohasma parameter file

    int requested_seed = -1; // Using specific seed for random number generator

    int c = 0;  // will be the remaining number of arguments during parsing
    while (c != -1){

        // Definition of the list of possible arguments: either "a_word" or "-X argument".
        // Other arguments (not from the predefined list), like parameter files, will also be retrieved at the end.
        static struct option long_options[] = {
            {"hypster",     no_argument,       0, 'h'},   // These options don’t set a flag.
            {"seed",    required_argument, 0, 's'},
            {"outputFolder",    required_argument, 0, 'o'},
            {0, 0, 0, 0}
        };

        // Parses the arguments one by one
        int option_index = 0;                     // the index of the identified argument will be put inside.
        // Don't forget to reput the list of -x allowed inside the third argument of getopt
        c = getopt_long (argc, argv, "hs:o:", long_options, &option_index);  // "as:" means, expects -a without argument or -s with argument
        switch (c)
        {
        case 0:
            {if (long_options[option_index].flag != 0)  // i.e. if there was a flag associated. Nothing to do
                cout << long_options[option_index].name << endl;
            if (optarg)
                cout << " with arg " << optarg << endl;
                break;}
        case 'h':
            {
            parfname = string(argv[optind]);
            cout << "   -h detected -> will read hyphasma parameter file:" << parfname<<endl;
            takeHyphasmaFile = true;
                break;}
        case 's':
            { requested_seed = atoi(optarg); //turn string into int
            cout << "   -s detecetd -> using seed: " << requested_seed << endl;
            //printf ("option -s with value `%s'\n", optarg);
                break;}
        case 'o':
            {outputFolder = optarg; //turn string into int

            cout << "   -o detected -> Using output folder: " << outputFolder << endl;
            //printf ("option -s with value `%s'\n", optarg);
                break;}
        case '?':
            { // getopt_long should have printed an error message.
                break;}
        default:
            {// no arguments within the list. Other arguments might be given (see next loop)
                break;}
        }
    }


    // Removing extension from par file if present

    if((parfname.size() > 0) && (!parfname.substr(parfname.size()-4,4).compare(string(".par")))){
        parfname = parfname.substr(0, parfname.size()-4);
        cout << "   ... The parameter file contained .par at the end => cutted it into " << parfname << endl;
    }

    //Elena: create folder using parameter file name.
    if((outputFolder.size() == 0) || (!outputFolder.compare(string("auto")))){
        outputFolder = parfname;
    }

    createFolder(outputFolder);


//// Print all parameters that is used to analyze file
//// Create simulation instance
////    if (Nof_simulations > 1) //Elena: Why not define this from script?
////    {
////    }
////    else {
        parameters currentParameterSet;
        currentParameterSet.convert_parameters();
       currentParameterSet.writeparameters(outputFolder+"/params.txt"); //Elena: Write parameters in output file since they might change per simulation.
        simulation Sim(currentParameterSet);
//        Sim.Simulation_ID=Nof_simulations;
//        initGC3D(argc, argv);
        Sim.simulate(*Sim.currentLattice, currentParameterSet);
////    }
    
    // Print all parameters that is used to analyze file
    // Create simulation instance
    //    if (Nof_simulations > 1) //Elena: Why not define this from script?
    //    {
    //    }
    //    else {
//            parameters currentParameterSet(takeHyphasmaFile,parfname);
//            currentParameterSet.convert_parameters();
//           currentParameterSet.writeparameters(outputFolder+"/params.txt"); //Elena: Write parameters in output file since they might change per simulation.
//            simulation Sim(currentParameterSet);
//    //        Sim.Simulation_ID=Nof_simulations;
//    //        initGC3D(argc, argv);
//            Sim.simulate(*Sim.currentLattice, currentParameterSet);
    //    }
    
    cout << "Starting OpenGL " << endl;
    
    //#check
    ofstream analyze;
    analyze.open(outputFolder+"/ana_ini.out");
    
// Initialize the visualization, later should put all the visualization together.
    analyze.close();
    return 0;

}
