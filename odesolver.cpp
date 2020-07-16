#include "odesolver.h"
#include <algorithm>
#include <sstream>
using namespace std;


void TableCourse::reset(){
    int nL = storage.size();
    for(int i = 0; i < nL; ++i) delete storage[i];
    storage.clear();
    attribut.clear();
    nbLignes = 0;
    headers.clear();
    headers.resize(nbVar+1);
}
void TableCourse::setHeader(int i, string title){
    // Header (0) = titre des attributs !!, commence à 1 !!
    if((i < 0) || (i > nbVar)) {cerr << "ERR:setHeades, out of bounds " << i << endl; return;}
    headers[i] = title;
}

TableCourse::TableCourse(int _nbVar) : nbVar(_nbVar), nbLignes(0) {reset();}

TableCourse::TableCourse(string fileToRead){
    read(fileToRead);
}

TableCourse::TableCourse(TableCourse* toCopy){
    nbVar = toCopy->nbVar;
    reset();
    nbLignes = toCopy->nbLignes;
    //if((int) toCopy->attribut.size() !=  nbLignes) {cerr << "ERR TableCourse::TableCourse(other TableCourse), attribut[] has wrong size" << endl; return;}
    //if((int) toCopy->storage.size() !=  nbLignes) {cerr << "ERR TableCourse::TableCourse(other TableCourse), storage[] has wrong size" << endl; return;}
    if((int) toCopy->headers.size() !=  nbVar + 1) {cerr << "ERR TableCourse::TableCourse(other TableCourse), headers[] has wrong size" << endl; return;}
    attribut.resize(nbLignes, 0);
    storage.resize(nbLignes, NULL);
    for(int i = 0; i < nbLignes; ++i){
        storage[i] = new vector<double>(nbVar, 0.0);
        attribut[i] = toCopy->attribut[i];
        if((int) toCopy->storage[i]->size() != nbVar) {cerr << "ERR TableCourse::TableCourse(other TableCourse), storage[" << i << "] has wrong size" << endl; return;}
        for(int j = 0; j < nbVar; ++j){
            (* (storage[i]))[j] = (* (toCopy->storage[i]))[j];
        }
    }
    for(int i = 0; i < nbVar + 1; ++i){
        headers[i] = toCopy->headers[i];
    }
}

double TableCourse::operator()(int vari, typeTime timej){
    if((vari < 0) || (vari >= nbVar)) return 0;
    if((timej < 0) || (timej >= nbLignes)) return 0;
    return (* (storage[timej]))[vari];
}

void TableCourse::addSet(double attr, vector<double> &toCopy){
    //if(!toCopy) {cerr << "ERR: TableCourse, given table doesn't exist (NULL)" << endl; return;}
    if((int) toCopy.size() != nbVar) {cerr << "ERR: TableCourse, given table has wrong size (" << toCopy.size() << ") instead of NbVar = " << nbVar << endl; return;}
    attribut.push_back(attr);
    vector<double> * vres = new vector<double>(nbVar, 0.0);
    for(int i = 0; i < nbVar; ++i){
        (*vres)[i] = toCopy[i];
    }
    storage.push_back(vres);
    nbLignes++;
}

string TableCourse::print(bool fileExportVersion){
    stringstream res;
    if(fileExportVersion) res << attribut.size() << "\t" << headers.size() - 1<< "\n";
    else res << "sizes : header " << headers.size() << " Attr" << attribut.size() << " storage" << storage.size() << " & nbLignes = " << nbLignes << endl;

    res << headers[0];
    for(int i = 1; i <= nbVar; ++i){
        //if(fileExportVersion) if(headers[i-1].size() < 8) res << "\t"; // to align the headers whti the following data (when >= 8, the tab is shifted)
        res << "\t" << headers[i];
    }
    res << "\n"; // << std::fixed << setprecision(6);
    for(int i = 0; i < nbLignes; ++i){
        res << attribut[i];
        for(int j = 0; j < nbVar; ++j){
            res << "\t" << (*(storage[i]))[j];
        }
        res << "\n";
    }
    return res.str();
}

void TableCourse::save(string fileName, string title){
    ofstream fichier(fileName.c_str(), ios::out);
    if(fichier){
        //fichier << "#   Time course for " << title << "\n";
        fichier << storage.size() << "\t" << headers.size()-1 << "\n";
        fichier << headers[0];
        for(int i = 1; i <= nbVar; ++i){
            fichier << "\t" << headers[i];
        }
        fichier << "\n";
        for(int i = 0; i < nbLignes; ++i){
            fichier << attribut[i];
            for(int j = 0; j < nbVar; ++j){
                double currentValue = (*(storage[i]))[j];
                if(std::isnan(currentValue)) {
                    fichier << "\tNaN";
                } else if(std::isinf(currentValue)) {
                    fichier << "\tinf";
                } else {
                    fichier << "\t" << currentValue;
                }
            }
            fichier << "\n";
        }
        fichier.close();
    }
}

void TableCourse::read(string fileName){

    nbLignes = 0;
    nbVar = 0;
    reset();
    ifstream fichier(fileName.c_str(), ios::out);
    if(fichier){
        fichier >> nbLignes >> nbVar;
        cout << "      Reading kinetics from file " << fileName << ", got NbLignes : " << nbLignes << ", nbVar :" << nbVar << endl;
        int nbLignesToRemember = nbLignes;//because of reset
        reset();    // to put everything to the good size
        string push;
        fichier >> push;
        setHeader(0, push);
        for(int i = 1; i <= nbVar; ++i){
            fichier >> push;
            setHeader(i, push);
        }
        stringstream tempBuffer;
        for(int i = 0; i < nbLignesToRemember; ++i){
            double attr;
            fichier >> attr;
            vector<double> values;

            // Shit happens with iostream : NaN and inf are not recognized by iostream => read string, and then converts to number via a stringstream
            double db;
            string trash;
            for(int j = 0; j < nbVar; ++j){
                fichier >> trash;
                if((!trash.compare("NaN")) || (!trash.compare("NAN")) || (!trash.compare("nan")) || (!trash.compare("+inf")) || (!trash.compare("-inf")) || (!trash.compare("inf"))){
                    trash.clear();
                    values.push_back(NAN);
                    //cout << "NaN" << "\t";
                } else {
                    tempBuffer << trash;
                    tempBuffer >> db;
                    tempBuffer.clear();
                    values.push_back(db);
                    //cout << db << "\t";
                }
                /*
                if(fichier >> db){
                    values.push_back(db);
                } else {
                    fichier.ignore(3);
                    string potentialNaN;
                    fichier >> potentialNaN;
                    fichier.unget();
                    //cout << "Success " << potentialNaN;
                }
                cout << db << "\t"; */
            }
            //cout << endl;
            addSet(attr, values);
        }
        fichier.close();
    } else cerr << "ERR:TableCourse::read(" << fileName << ", file not found - note that the tablecourse was cleared\n";
}


vector<double> TableCourse::getTimeCourse(int var){
    vector<double> res = vector<double>();
    if((var < 0) || (var > nbVar)){cerr << "ERR::TableCourse, getTimeCourse(" << var << "), out of bounds. nbVars = " << nbVar << endl; return res;}
    //res.resize(nbLignes);
    for(int i = 0; i < nbLignes; ++i){
        double currentVal = (*(storage[i]))[var];
        if(!(std::isnan(currentVal) || std::isinf(currentVal))){
            res.push_back((*(storage[i]))[var]);
        }
    }
    return res;
}
vector<double> TableCourse::getTimePoints(int var){
    if(var < 0) {
        vector<double> res1 = vector<double>(attribut);
        return res1;
    }
    vector<double> res = vector<double>();
    if(var > nbVar){cerr << "ERR::TableCourse, getTimePoints(" << var << "), out of bounds. nbVars = " << nbVar << endl; return res;}
    for(int i = 0; i < nbLignes; ++i){
        double currentVal = (*(storage[i]))[var];
        if(!(std::isnan(currentVal) || std::isinf(currentVal))){
            res.push_back(attribut[i]);
        } //else cout << "!";
    }
    return res;
}


// timePoints should be ordered !!!
TableCourse TableCourse::subKinetics(vector<int> timePoints, vector<string> namesVariables /* same names than in the tablecourse */){
    cout << "  --> Extracting subKinetics ..." << endl;

    int nbTP = timePoints.size();
    int nbWantedVars = namesVariables.size();
    if(nbTP == 0) {cerr << "ERR: TableCourse::subKinetics, no time point given " << endl; return TableCourse(0);}
    sort(timePoints.begin(), timePoints.end());     // the time points should be sorted !!

    // tables to extract the variables we want to keep in the subKinetics:
    vector<int> placeHeaderInNewKin;    //STARTS AT ZERO for each variable placeHeaderInNewKin[header variable position] -> will go to position x in the subkinetics.-1 if not taken
    placeHeaderInNewKin.resize(headers.size(), -1);
    vector<string> listeNamesNewTable;     // list of headers of the subkinetics   listeIDsNewTable[header position in subkinetics] -> name

    if(nbWantedVars == 0){
        for(int i = 1; i < (int) headers.size(); ++i){
           placeHeaderInNewKin[i] = i-1;    // STARTS AT ZERO, (indice for data vector later)
           listeNamesNewTable.push_back(headers[i]);
        }
    }
    else {
        stringstream errs;
        for(int i = 0; i < nbWantedVars; ++i){
            string wantedVar = namesVariables[i];
            bool done = false;
            for(int j = 1; j < (int) headers.size(); ++j){
                if(! headers[j].compare(wantedVar)) {   // if the jth variable
                    //cout << " compare " << headers[j] << " and " << wantedVar << endl;
                    if(!done) {listeNamesNewTable.push_back(wantedVar);
                        if(placeHeaderInNewKin[j] >= 0) cerr << "ERR : TableCourse::subKinetics, you're asking twice the same variable !!\n";
                        placeHeaderInNewKin[j] = listeNamesNewTable.size() - 1;} // i.e. the last position
                    done = true;
                }
            }
            if(!done) errs  << "\t" << namesVariables[i];
        }
        if(errs.str().size() > 0) cerr << "WRN : TableCourse::subKinetics, variables not found (probably not simulated) in kinetics" << errs.str() << endl;
    }

    int nbVarReduced = listeNamesNewTable.size();
    TableCourse Tbres(nbVarReduced);
    Tbres.setHeader(0,"time");
    for(int i = 0; i < nbVarReduced; ++i){
        Tbres.setHeader(i+1, listeNamesNewTable[i]);
    }

    vector<double> data;
    data.resize(nbVarReduced, 0.0);

    int pointerTimePoint = 0; // already checked NBTP > 0
    for(int i = 0; (i < nbLignes) && (pointerTimePoint < nbTP); ++i){
        if(attribut[i] == timePoints[pointerTimePoint]){ // if timepoint is required
            pointerTimePoint++;
            for(int j = 0; j < nbVar; ++j){
                if(placeHeaderInNewKin[j+1] >= 0)   // le +1 qui tue ...
                    data[placeHeaderInNewKin[j+1]] = (*(storage[i]))[j];
            }
            if((int) data.size() != nbVarReduced) cerr << "ERR : TableCourse::subKinetics, for time point " << attribut[i] << " , some variables where not found. Should not happen !!\n";
            Tbres.addSet(attribut[i], data);
        }
    }
    if(pointerTimePoint < nbTP) cerr << "WRN :  TableCourse::subKinetics, some time points were not found (first not found is : " << timePoints[pointerTimePoint] << ")\n";

    return Tbres;
}


double boundedLog(double x){
    return 6.90776 + max(log(max(1e-8,x)), -6.90776);
}

double fitnessFunction(double valueToCompare, double refExpectedValue, double StdDev){
#ifdef SQUARE_COST
    return (refExpectedValue- valueToCompare) * (refExpectedValue - valueToCompare);
#endif
#ifdef SQUARE_COST_STD
    if(fabs(StdDev < 1e-8)) return (refExpectedValue - valueToCompare) * (refExpectedValue - valueToCompare);
    return (refExpectedValue- valueToCompare) * (refExpectedValue - valueToCompare) / StdDev;
#endif
#ifdef LOG_COST
    return fabs(boundedLog(fabs(refExpectedValue)) - boundedLog(fabs(valueToCompare)));
#endif
#ifdef LOG_COST_STD
    //if(fabs(StdDev < 1e-8)) return boundedLog(fabs(refExpectedValue- valueToCompare));
    return fabs(boundedLog(fabs(refExpectedValue)) - boundedLog(fabs(valueToCompare)));
#endif
#ifdef PROPORTION_COST
    if(fabs(fabs(refExpectedValue) < 1e-3)) return min(2.0, fabs(valueToCompare) * 1000);
    return min(2.0, fabs((refExpectedValue- valueToCompare) / refExpectedValue));
#endif
    cerr << "No cost function policy defined. Please put #define SQUARE_COST, SQUARE_COST_STD, LOG_COST, LOG_COST_STD, PROPORTION_COST in a file included by evaluator.cpp\n";
    return -1;
}

void testFitnessFunction(){
    vector<double> exp = {0, 0, 0.1, 0.1, 0.1, 1, 1, 1,      1.5,   1.5, 2,     5, 5, 5, 5,     10, 0.00001, 0.00001, 0.001};
    vector<double> mes = {0, 1, 0.1, 0.3, 1  , 1, 2, 0.0001, 1000,  0,   1.5,   3, 4, 5, 10,    100,0.0001,  0.001,    1};
    int L = exp.size();
    for(int i = 0; i < L; ++i){
        cout << std::setprecision(8) << "Exp= " << exp[i] << ",   \tmes= " << mes[i] << "   \tcost= " << fitnessFunction(mes[i], exp[i], 0) << "\t\t" << boundedLog(exp[i]) << "," << boundedLog(mes[i]) << endl;
    }
}





#ifdef USE_BOOST_SOLVER
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;
#endif

/// ======================= ODE Solver (defined with the same interface than BOOST) ===============================

/// This class is the interface to use boost::integrate functions (a class where the operator () gives the derivatives)
/// I don't use boost, but I keep this formalism so it's easy to use boost when we want

class mySim {
public:
    odesolver* mm;
    mySim( odesolver* _mm ) : mm(_mm) {}
    void operator()(const vector<double> &_x, vector<double> &_dxdt, const double t){
        mm->derivatives(_x,_dxdt,t);

    }
};

/// a function to evalute the percent of error between two vectors (for comparing dt time step versus two times dt/2)
double evalErr(vector<double> &v1ref, vector<double> &v2, int size){
    double res = 0;
    for(int i = 0; i < size; ++i)
        if(v1ref[i] != 0) res = max(res, fabs((v2[i] - v1ref[i]) / (v1ref[i])));
    return res/size;
}


/// Integrates the ODEs (in the model in mySim), updating the start vector, from t to nxt, with the minimum dt dtmin (and maximum dt dtmin * 500, starting from dtmin * 10)
size_t myintegrate( mySim BS , vector<double> &start , double t , double nxt , double dtmin ){
    int size = start.size();
    vector<double> oneshotdxdt  = vector<double>(size, 0.0);
    vector<double> oneshotpoint = vector<double>(size, 0.0);
    vector<double> halfdxdt     = vector<double>(size, 0.0);
    size_t nbSteps = 0;

    double dtlocal = min(dtmin*10, (nxt-t) / 10); //hrs.
    if(dtlocal == 0) return 0;  // simulation finished
    while(t < nxt){
        // if integrating one step of dtlocal
        BS(start, oneshotdxdt, t); //Elena: (x, dx/dt, t)
        for(int i = 0; i < size; ++i){
            oneshotpoint[i] = start[i] + dtlocal * oneshotdxdt[i];
            //if(isnan(oneshotpoint[i]) || (isinf(oneshotpoint[i])) || (fabs(oneshotpoint[i]) > 1e12)) {/*cerr << "Divergence ! " << endl; */ return 0; }
        }

        // versus if integrating two steps of dtlocal/2 (the result from this technique is the one kept (more accurate than one step), but the next dt will be chosen depending if one step is accurate enough (close to two steps) or not
        //BS(start, halfdxdt, t);
        for(int i = 0; i < size; ++i)
            {start[i] = start[i] + 0.5* dtlocal * oneshotdxdt[i];}
        BS(start, halfdxdt, t + 0.5* dtlocal);
        for(int i = 0; i < size; ++i){
            start[i] = start[i] + 0.5* dtlocal * halfdxdt[i];
            if((start[i] < MIN_VAL) && (start[i] > (-1.0 * MIN_VAL))) start[i] = 0;
            }  // maybe better to do it every 10 or 20 loops ...
        t += dtlocal;
    }
    return nbSteps;
}


#ifdef USE_BOOST_SOLVER
runge_kutta4 <vector<double>> rk4 = runge_kutta4<vector<double>> ();

/// Integrates the ODEs (in the model in mySim), updating the start vector, from t to nxt, with the minimum dt dtmin (and maximum dt dtmin * 500, starting from dtmin * 10)
size_t RK4integrate( mySim BS , vector<double> &start , double t , double nxt , double dtmin ){

    bool diverge = false;
    int size = start.size();
    vector<double> oneshotdxdt  = vector<double>(size, 0.0);
    vector<double> oneshotpoint = vector<double>(size, 0.0);
    vector<double> halfdxdt     = vector<double>(size, 0.0);
    size_t nbSteps = 0;

    double dtlocal = min(dtmin*10, (nxt-t) / 10);
    if(dtlocal == 0) return 0;  // simulation finished
    while(t < nxt){
        // if integrating one step of dtlocal

        rk4.do_step( BS , oneshotpoint , t , dtlocal );
        /*BS(start, oneshotdxdt, t);
        for(int i = 0; i < size; ++i){
            oneshotpoint[i] = start[i] + dtlocal * oneshotdxdt[i];
            //if(isnan(oneshotpoint[i]) || (isinf(oneshotpoint[i])) || (fabs(oneshotpoint[i]) > 1e12)) { return 0; }
        }*/

        // versus if integrating two steps of dtlocal/2 (the result from this technique is the one kept (more accurate than one step), but the next dt will be chosen depending if one step is accurate enough (close to two steps) or not

        //implicitely (already called before) BS(start, halfdxdt, t + 0.5* dtlocal);
        /*for(int i = 0; i < size; ++i)
            {start[i] = start[i] + 0.5* dtlocal * oneshotdxdt[i];}*/
        rk4.do_step( BS , start , t , 0.5* dtlocal );               // In-place transformation of inout

        //BS(start, halfdxdt, t + 0.5* dtlocal);
        rk4.do_step( BS , start , t , 0.5* dtlocal );


        for(int i = 0; i < size; ++i){
            // start[i] = start[i] + 0.5* dtlocal * halfdxdt[i];
            if((start[i] < MIN_VAL) && (start[i] > (-1.0 * MIN_VAL))) start[i] = 0;
            if(fabs(start[i]) > STOP_VAL) diverge = true;
        }                                       // maybe better to do it every 10 or 20 loops ...
        t += dtlocal;

        double compare = evalErr(start, oneshotpoint, size);
        if(compare < 0.002) dtlocal *= 1.15; else dtlocal *= 0.5;
        if(dtlocal < dtmin) {dtlocal = dtmin;}          // minimum dt
        if(dtlocal > (nxt-t)) dtlocal = (nxt-t);        // should not exceed the stopping time point
        if(dtlocal > (500*dtmin)) dtlocal = 500*dtmin;  // should not exceed 500 * dtmin
        nbSteps += 2;                                   // for information, reminds the number of steps

        if(diverge) return nbSteps;
    }
    return nbSteps;
}
#endif

/// function to simulate from t to t+sec_max. The value of odesolver::t is updated
void odesolver::simulate(double sec_max, Evaluator* E){
//    if(!E) cerr << "No Evaluator\n";
    if(saveKinetics) newCinetiqueIfRequired(); else deleteCinetique();
    mySim BS = mySim(this);

    double nxt = 0;         // nxt will be the next stopping point each time (for recording kinetics or for getting specific value)
    double t_init = t;
    size_t steps = 0;

    while(nxt < t_init + sec_max){       // when finished, 1e8

        // finds the next wanted time point
        nxt = t_init + sec_max;
        if(saveKinetics){nxt = min(nxt, print_all_secs * (double) ((int) (t / print_all_secs + 1e-8) + 1));} // the 1e-8 is due to compensate the errors. Indeed, (0.03 / 0.01) = 3 but (int) (0.03 / 0.01) was giving 2 ...
        #ifdef USE_EVALUATORS
        if(E) nxt = min(nxt,E->nextPoint());
        #endif


#ifdef USE_BOOST_SOLVER
// Version 1 : boost adaptive runge kuutta => freezes
//typedef controlled_runge_kutta< runge_kutta_dopri5<vector<double>> > ctrl_rkck54;
//steps += integrate_adaptive( ctrl_rkck54(), BS, val , (double) t , (double) nxt , (double) dt /*, my_observer */ );

// Version 2 : same adaptive as aeuler but with rk4 inside
steps += RK4integrate( BS , val , (double) t , (double) nxt , (double) dt );     // this line can be directly replaced by the boost solver ...

#else
        steps += myintegrate( BS , val , (double) t , (double) nxt , (double) dt );     // this line can be directly replaced by the boost solver ...
#endif

        t = nxt;    // this might not be right if the simulation stopped before ...

        if(saveKinetics){if( fabs((double) ((int) (t / print_all_secs + 1e-8)) - t / print_all_secs) < 1e-8) save_state(t);} // to account for division errors
    }
    //{ // in case, finishes the simulation (but in theory, nxt should have reached t+sec_max)
    //    steps += myintegrate( BS , val , (double) t , (double) t_init +sec_max , dt );
    //    t = t_init + sec_max;
    //}
}








// ======================================= 2 Basic definitions ==========================

odesolver::odesolver(int _nbVars,  int _nbParams) : nbVars(_nbVars), nbParams(_nbParams), penalities(0), saveKinetics(false), cinetique(NULL) {
    t = 0; stopCriterionReached = false; dt = 0.02; print_all_secs = 1200;	// basal values for simulation, can be changed by the subclass
    params.clear();         params.resize(nbParams, 0.0);
    val.clear();            val.resize(nbVars, 0.0);
    names.clear();          names.resize(nbVars, string(""));
    init.clear();           init.resize(nbVars, 0.0);
    extNames.clear();       extNames.resize(nbVars, 0);
    paramNames.clear();     paramNames.resize(nbParams, string(""));
    backSimulated.clear();}
void odesolver::loadParameters(string file_name){
    cout << "      -> loading paramter set from " << file_name << endl;
    ifstream fichier(file_name.c_str(), ios::in);
        if(fichier){int nb_p = 0; fichier >> nb_p;
        for(int i = 0; i < min(nb_p, (int) nbParams); ++i){
            double tampon = 0.0;
            fichier >> tampon;
            if(tampon != 0) params[i] = tampon;
            //cout << "Loading p" << i << " = " << params[i] << endl;
        }}
    loadParametersDone();}
// Virtual functions that SHOULD be reimplanted
void odesolver::saveParameters(string file_name){
    ofstream fichier(file_name.c_str(), ios::out);
    if(fichier){
        fichier << nbParams << endl;
        for(int i = 0; i < nbParams; ++i){
            if(i > 0) fichier << "\t";
            fichier << getParam(i);
        }
        fichier << endl;
        fichier.close();
    } else {cerr << "odesolver::saveParameters(), file not found : " << file_name << endl;}
}
// Virtual functions that SHOULD be reimplanted
void odesolver::derivatives(const vector<double> &x, vector<double> &dxdt, const double t){
    cerr << "ERR : the function derivatives has not been implemented properly in the daughter class\n";}
void odesolver::initialise(long long ){
    cerr << "WRN : Nothing specified for initialization\n";}
void odesolver::setBaseParameters(){
    cerr << "ERR : the function setBaseParameters has been called but not implemented in the daughter class\n";}				// gives a correct set of parameters, if you don't load other set.
void odesolver::initialiseDone(){
    if(! parametersLoaded) cerr << "ERR : Please provide a full parameter set to the odesolver before simulating. Details : you called initialise() without loading/setting parameters before initializing ! - you can not perform simulation ! Solutions : call 'loadParameters(file)', 'setParameters(full set)' or 'setBaseParameters()\n";
    deleteCinetique(); penalities = 0;stopCriterionReached=false;}
void odesolver::setBaseParametersDone(){parametersLoaded = true;}
void odesolver::loadParametersDone(){parametersLoaded = true;}
double odesolver::max(double a, double b){if(a > b) return a; else return b;}
string odesolver::print(){
    stringstream res;
    // I preferred to remove any use depending on the namespaces N:: and Back:: so that the odesolver class is general and independent of them.
    res << "   Nb Params " << nbParams << endl;
    res << "   Nb Vars   " << nbVars << endl;
    res << "   Internal Variables : " << endl;
    if((int) names.size() == nbVars){
        for(int i = 0; i < nbVars; ++i){
            res << "\t" << i << "\t" << names[i] << "\t";
            if(extNames[i] >= 0) res << "Globally Known As Var ID=: (" << extNames[i] << ")"; //<< GlobalName(extNames[i]) <<
            else res << "Internal only ";
            res << endl;
        }
     } else res << "ERR: odesolver::print, wrong size for names[] - there are not enough names found for species" << endl;
    res << "    Current parameter values\n";
    for(int i = 0; i < nbParams; ++i){res << "   " << i << "\t" << params[i] << endl;}
    return res.str();
}
string odesolver::printVarNames(double _t){
    stringstream res;
    res << "\n" << _t <<" h\t";
    for(int i = 0; i < nbVars; ++i){
        res << names[i] << "\t";
    } res << endl;
    return res.str();
    }
string odesolver::printVarValues(double _t){
    stringstream res;
    res << _t << "h" <<  "'\t" << std::fixed;
    for(int i = 0; i < nbVars; ++i){
        res << setprecision(4) << val[i] << "\t";
    } res << "\n";
    return res.str();}
string odesolver::printParValues(){
    stringstream res;
    res << "Parameter values (\t" << nbParams << " params)\n";
    for(int i = 0; i < nbParams; ++i){
        res << i << "\t" << params[i] << "\t";
    } res << endl;
    return res.str();}



// ============================ 3 - Managing the names of variables and I/O to variables ========================

// virtual functions that don't need to be overrided (provided extNames and backSimulated are filled).
void odesolver::setValue(int idGlobalVariable, double value){
    int id = internValName(idGlobalVariable); if((id >= 0) && (id < nbVars)) val[id] = value; else cerr << "odesolver::setValue(" << idGlobalVariable << "," << value << "), this species ID is not implanted !! \n";}
void odesolver::addValue(int idGlobalVariable, double value){
    int id = internValName(idGlobalVariable); if((id >= 0) && (id < nbVars)) val[id] += value; else cerr << "odesolver::addValue(" << idGlobalVariable << "," << value << "), this species ID is not implanted !! \n";}
double odesolver::getValue(int idGlobalVariable){
    int id = internValName(idGlobalVariable); if((id >= 0) && (id < nbVars)) return val[id]; else {cerr << "odesolver::getValue(" << idGlobalVariable << "), this species ID is not implanted !! \n"; return -1;}}
bool odesolver::isSimuBack(long long idGlobalBack){
    long long nb = backSimulated.size();
    for(int i = 0; i < nb; ++i){if(backSimulated[i] == idGlobalBack) return true;}
    return false;}
bool odesolver::isVarKnown(int idGlobalVariable){return (internValName(idGlobalVariable) >= 0);}
int odesolver::internValName(int idGlobalVariable){
    int nb = extNames.size();
    for(int i = 0; i < nb; ++i){if(extNames[i] == idGlobalVariable) return i;}
    return -1;}
int odesolver::getNbParams(){return nbParams;}
void odesolver::setParam(int i, double v){
    if((i < 0) || (i >= nbParams)) {cerr << "odesolver::setParam(" << i << "," << v << "), index out of bounds. NbParams = " << nbParams << endl; return;}
    params[i] = (v > 0? v : -v);}   // fuck the system !!
double odesolver::getParam(int i){
    if((i < 0) || (i >= nbParams)) {cerr << "odesolver::getParam(" << i << "), index out of bounds. NbParams = " << nbParams << endl; return 0.0;}
    return params[i];}              // fuck the system !!
string odesolver::getName(int internalID){
    if((internalID < 0) || (internalID >= nbVars)) {return string("");}
    return names[internalID];}
int odesolver::getGlobalID(int internalID){
    if((internalID < 0) || (internalID >= nbVars)) {return -1;}
    return extNames[internalID];}
int odesolver::getNbVars(){
    return nbVars;}
vector<double> odesolver::getInit(){
    return init;}
double odesolver::getT(){return t;}
vector<double> odesolver::getParameters(){
    return params;}
void odesolver::setParameters(vector<double> &newParamSet){
    if((int) newParamSet.size() != nbParams) {cerr << "ERR- odesolver::setParameters(vector), wrong size for the given vector (" <<  newParamSet.size() << ") instead of " << nbParams << " parameters\n"; return;}
    for(int i = 0; i < nbParams; ++i){params[i] = newParamSet[i];}
    parametersLoaded = true;
}

void odesolver::action(string name, double parameter){cerr << "ERR- odesolver::action(string, double), this function should and was not implemented in the odesolver subclass\n";}
void odesolver::action(string name, vector<double> parameters){cerr << "ERR- odesolver::action(string, vector<double>), this function should and was not implemented in the odesolver subclass\n";}
// ============================ 4 - Managing the kinetic mode : records every xx secs ==============================

void odesolver::setPrintMode(bool z, double _print_all_secs){
    saveKinetics = z;
    if(_print_all_secs > 0) print_all_secs = _print_all_secs;}
void odesolver::newCinetiqueIfRequired(){
    if(cinetique != NULL) return;
    cinetique = new TableCourse(nbVars);
    cinetique->setHeader(0, string("Time"));
    for(int i = 0; i < nbVars; ++i){
        cinetique->setHeader(i+1, names[i]);}}
void odesolver::deleteCinetique(){
    if(cinetique){delete cinetique; cinetique = NULL;}}
TableCourse odesolver::getCinetique(){
    if(cinetique) return *cinetique; else return TableCourse(nbVars);}
void odesolver::save_state(double _t){
    if(saveKinetics && cinetique){cinetique->addSet(t, val);}}
