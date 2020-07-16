//#ifndef MAIN_H
//#define MAIN_H
//#include "lattice.h"
//#include "parameters.h"
//#include "chemokines3d.h"
//#include "linkHyphasma/setparam.h"
//#include <iostream>
//#include <list>
//#include <getopt.h> /// Library to filter arguments from the command line to the main
//
//using namespace std;
//
//int main(int argc, char** argv){
//
//    cout << "Starting Mafalda " << endl;
//
//    ///////////////////////////////////// Parsing arguments given from command line ///////////////////////////////7
//    // See https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html#Getopt-Long-Option-Example
//    int takeHyphasmaFile = false;
//    int requested_seed = -1;              // one argument that can be given: a seed
//
//    int c = 0;                            // will be the remaining number of arguments during parsing
//    while (c != -1){
//
//        // Definition of the list of possible arguments: either "a_word" or "-X argument".
//        // Other arguments (not from the predefined list), like parameter files, will also be retrieved at the end.
//        static struct option long_options[] = {
//            {"hyphasma", no_argument,       &takeHyphasmaFile, 1}, // Options to set a flag.
//            {"hypster",     no_argument,       0, 'h'},   // These options don’t set a flag.
//            {"seed",    required_argument, 0, 's'},
//            {0, 0, 0, 0}
//        };
//
//        // Parses the arguments one by one
//        int option_index = 0;                     // the index of the identified argument will be put inside.
//        c = getopt_long (argc, argv, "as:", long_options, &option_index);  // "as:" means, expects -a without argument or -s with argument
//        switch (c){
//        case 0:
//            if (long_options[option_index].flag != 0)  // i.e. if there was a flag associated. Nothing to do
//                cout << long_options[option_index].name << endl;
//            break;
//            if (optarg)
//                cout << " with arg " << optarg << endl;
//            break;
//            case 'h':
//                cout << "-h detected for HYPHASMA MOUAHAHAHAHA" << endl;
//                takeHyphasmaFile = true;
//                break;
//        case 's':
//            requested_seed = atoi(optarg);
//            cout << "Using seed " << requested_seed << endl;
//            //printf ("option -s with value `%s'\n", optarg);
//            break;
//        case '?':
//            // getopt_long should have printed an error message.
//            break;
//        default:
//            // no arguments within the list. Other arguments might be given (see next loop)
//            //cout << "no argument given" << endl;
//            break;
//        }
//    }
//    // Additional arguments, that are not in the 'official' list options (the parameter file for instance).
//    int nArguments = argc - optind;
//    if (optind < argc){
//        cout << "Detected additional parameters: ";
//        for(int i = optind; i < argc; ++i){
//            cout << "\t" << argv[i];
//        }
//        cout << endl;
//    }
//    if(nArguments != 1) {
//        cerr << "You should give one (and only one) parameter file when you run Mafalda" << endl;
//    }
//
//    takeHyphasmaFile = true; //#Recheck Danial: Should be based on arguments of the main function. For the moment I set it to true, so mafalda reads hyphasma parameter file and has the built-in parameter set and it is possible to use either of them. This should be changed to a 3-way option where user can decide to use hyphasma-file mafalda-file or built-in set. Elena, if you don't need to import parameters from a file, we can remove this part and have only hyphasma and built-in.
//
//    parameters* currentParameterSet = new parameters();
//    if(!takeHyphasmaFile){
//        cout << "Reading parameter file" << argv[optind] << "with Mafalda" << endl;
//        currentParameterSet->readparameters(argv[optind]);
//    } else {
//        cout << "Reading parameter file" << argv[optind] << "with hyphasma" << endl;
//        hyphasmaParameter par; // note: this is a Hyphasma parameter set
//        bool done = par.wahl(argv[optind],true,true); // the command for hyphasma to read a parameter file
//        if(!done) cerr << "ERR: hyphasma refused to load the parameter file" << argv[optind] << endl;
//        currentParameterSet->matchFromHyphasma(par); // converts the parameters from hyphasma into parameters here
//
//        // Philipe: we should put all values in minutes
//    }
//    currentParameterSet->writeparameters(string("UsedParameterSet.txt"));
//    simulation Sim(*currentParameterSet);
//    Sim.simulate(*Sim.currentLattice, *currentParameterSet);
//
//    return 0;
//}
//
//#endif // MAIN_H

