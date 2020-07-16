#ifndef MAIN_H
#define MAIN_H
#include "lattice.h"
#include "parameters.h"
#include "chemokines3d.h"
#include "events.h"
#include "output.h"
#include "setparam.h"
#include <iostream>
#include <list>
using namespace std;
extern string outputFolder;
class simulation
{
public:
    simulation(parameters &p); // Reset list of cells and Initializes Simulation
    ~simulation();
    // parameters* currentParameterSet; Philippe: possible improvement to copy the parameter set. needs aq copy function, tricky
    lattice* currentLattice;
    vector <Stromal_cell*> ListSC;
    vector <FDC*> ListFDC;
    vector <B_cell*> ListB_cell;
    vector <T_cell*> ListT_cell;
    vector <Memory_cell*> ListM_cell;
    vector <Plasma_cell*> ListP_cell;
    void BCinflux(double time,parameters &p, lattice &l);
    void Calc_TC(parameters &p, lattice &l, vector<vector3D> &redo_list);
    void Calc_BC(double t,parameters &p,lattice &l, vector<vector3D>&redo_list,vector<int>&going_to_delet);
    void Visualise(double t, parameters &p);
    void Calc_Out(double t,parameters &p,lattice &l, vector<vector3D>&redo_list);
    void clean_dead_cells(lattice &l);
    //This function is to transfer B cells which differentiate to plasma cell from B cell list to plasma cell list
    void transfer_plasma_from_Bcell_list(double t,parameters &p,lattice &l, vector<vector3D>&redo_list);
    vector<int> going_to_delet;
//    short Simulation_ID;
    void InitialCells(lattice& l, parameters& p);
    void simulate(lattice& l, parameters& p); // Executes procedures in GC simulation every time step
    output* currentOutput;
    events* EventOutput;
    stringstream sim_output;
};



#endif // MAIN_H
