#ifndef ODESOLVER_H
#define ODESOLVER_H
//#define USE_BOOST_SOLVER
// Philippe Robert, 19-09-2014 - philippe-robert.com
// odesolver.h and odesolver.cpp  are using the TableCourse structure

#define typeTime double
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>      // std::setprecision
#include <cmath>
#include <fstream>
#include <sstream>
#include <utility>

using namespace std;

/// Data structure to store the discrete kinetics of several variables (ideal for plotting simulations)
/// as a table (meaning, for a list of time-points, store the values of all these selected variables)


struct TableCourse {
    vector<double> attribut;            /// List of time points (left column)
    vector< vector<double> * > storage; /// data (2D) storage[i][j] = value at time point i (t = attribut[i]) of variable j (whose name is headers[j])
    vector<string> headers;             /// Names of the variables (columns), size = nbVars + 1. The first header is the name of the table (or title of the left column (ex: time))
    int nbVar;                          /// nb of variables
    int nbLignes;                       /// nb of time points

    /// Reading a table from another table or from a file
    TableCourse(TableCourse* toCopy);   /// copy from another one
    TableCourse(string fileToRead);     /// reads from a file.
            /// SYNTAX of a TableCourse file (read or exported):
            ///
            /// nbRows(time-pts)     nbColumns(nbVars)
            /// headerLeft      headerVar1      headerVar2      ...
            /// time1           valVar1         valVar2         ...
            /// time2           valVar1         valVar2         ...
            /// ...

    /// Manually making and writing into a table
    TableCourse(int _nbVar);            /// creates empty, number of columns
    void setHeader(int i, string title);                    /// Danger : starts at index 1 for variables. Header[0] = titre of the table
    void addSet(double attr, vector<double> &toCopy);       /// Add a line (i.e. the value of each variable at the new time (attr)
    void reset();

    /// Additional I/O :
    void read(string fileName);
    void save(string fileName, string title = string(""));
    vector<double> getTimeCourse(int var);
    vector<double> getTimePoints(int var = -1);
    double operator()(int vari, typeTime timej);
    TableCourse subKinetics(vector<int> timePoints, vector<string> namesVariables = vector<string>()/* same names than are in the tablecourse headers */);
    //void print();
    string print(bool fileExportVersion = true);
};



#define STOP_VAL 1e8
#define MIN_VAL 1e-10
#define DEFAULT_BACKGROUND_VALUE    1

/// ================= A general mother class for ODE simulations. =========================

struct odesolver {
    virtual ~odesolver(){}
    string name;                                            /// Name of the model (to separate different subclasses ...)
    odesolver(int _nbVars,  int _nbParams);
    int nbVars;
    int nbParams;
    vector<double> init;                                    /// for initial values of variables
    vector<double> params;                                  /// for parameter values
    vector<string> names;                                   /// for the names of the variables
    vector<string> paramNames;

    ///     for running a simulation, these variables will be used
    double t;                                               /// updated during simulation / advised to do 't=0;' in the subclass function initialise
    vector<double> val;                                     /// for variables at time t

public:
    /// Time evolution for dt at time t, initialise and base parameter values
    virtual void derivatives(const vector<double> &x, vector<double> &dxdt, const double t);/// computes the derivatives (at t) from position x
    virtual void setBaseParameters();                                                      /// gives a correct set of parameters, if you don't load other set.
    virtual void initialise(long long _background = DEFAULT_BACKGROUND_VALUE);            /// initialise, (parameters should be set before initializing) - sets t = 0 if necessary.
                                                                                         /// the background parameter allows to give options of simulation (such as deficient mice)
    virtual void loadParameters(string file_name);          /// reads parameters from a text file with format : "NB_PARAMS\tparam1\tparam2\tparam3 ..."
    virtual void saveParameters(string file_name);          /// writes parameters from a text file with format : "NB_PARAMS\tparam1\tparam2\tparam3 ..."

    /// A function to perform personnalized actions by the model, without a specific function, but by a name of action, (ex : adding a cytokine at day 2)
    virtual void action(string name, double parameter);
    virtual void action(string name, vector<double> parameters);

protected:
    void initialiseDone();
    void setBaseParametersDone();
    void loadParametersDone();


public:
/// Now, the interests of using the odesolver mother class comes :
/// whatever the sub-class model, the following functions are implemented and can be used to simulate, or do interesting manipulations of the model.
    typedef void Evaluator;
/// Functions for simulating (automatically calling the sub-class derivatives function, using the solver)
    void simulate(double _sec_max, Evaluator* E = NULL);    /// Simulates from t (not changed), to t + _sec_max, by calling one_step
                                                            /// the evaluator, if given, is a structure that says when to store data from the simulation, and avoids to store everything.
                                                            /// if the simulation diverges, it is stopped and a penalty is computed in the field penalities (of the mother class)
    double dt;                                              /// minimum dt used for simulations. The tunable dt will start at dt*10 and be limited between dt and dt*100
    double penalities;                                      /// automatically updated by the mother class : put to zero when initialiseDone() called and increased when the simulation diverges.
/// To get simulations as tables of values (kinetics) of simulation(before a simulation activate the kinetic mode).
/// Then, for every simulation, the simulation data will be put in a new table, that can be retrieved by getCinetique.
/// Note: that 'initialiseDone' clears the current table.
    bool saveKinetics;                                      /// sets the mode : simulates a kinetics to record or just simulate what you ask without stopping
    double print_all_secs;                                  /// frequency of saving/displaying data for kinetics
    void setPrintMode(bool z, double _print_all_secs);      /// to set the 'recording mode' to ON/OFF, and the frequency of recording
    TableCourse getCinetique();                             /// then, each time initialiseDone() is called, the kinetics is cleared.
    vector<double> getInit();
    double getT();
    string print();
    string printVarNames(double _t);
    string printVarValues(double _t);
    string printParNames();
    string printParValues();

/// Accessing variables from an external name, and defining backgrounds (even if the model has its internal variable names, you might want to say that this vairable represents this cytokine).
protected:
    vector<int> extNames;             /// for each variable, the model can give a global ID
    vector<long long> backSimulated; /// the model can list in here the IDs of backgrounds it is able to simulate. This info can be used inside the sub-class functions (network...)

/// Variables can be accessed with the 'external ID' with the following functions
public:
    virtual void   setValue(int idGlobalVariable, double val);      /// to modify the value of a variable from the global ID of it
    virtual void   addValue(int idGlobalVariable, double val);      /// to modify the value of a variable from the global ID of it
    virtual double getValue(int idGlobalVariable);                  /// to get the value of a variable from its global ID
    virtual bool   isSimuBack(long long idGlobalBack);              /// to know if a background can be simulated by the model
    virtual bool   isVarKnown(int idGlobalVariable);                /// to know if a variable can be simulated by the model (from its global ID)
    virtual int    internValName(int idGlobalVariable);             /// to know the index of a variable inside the tables

/// Working directly with the parameters of a model
    int getNbParams();
    void setParam(int i, double v);
    double getParam(int i);
    vector<double> getParameters();
    void setParameters(vector<double> &newParamSet);

/// Get information on the (internal) variables.
    int getNbVars();
    string getName(int internalID);
    int getGlobalID(int internalID);

/// Managing kinetics (automatically used)
private:
    TableCourse* cinetique;
    void deleteCinetique();
    void newCinetiqueIfRequired();
    void save_state(double _t);         /// the simule function call this to record data of a time point into the Cinetique.
protected:
    double max(double, double);
    bool parametersLoaded;
    bool checkDivergence();
    bool stopCriterionReached;
};


#endif // ODESOLVER_H


