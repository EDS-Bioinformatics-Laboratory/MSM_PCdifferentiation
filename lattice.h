#ifndef LATTICE_H
#define LATTICE_H
#include "chemokines3d.h"
#include "cell.h"
#include "random.h"
#include <vector>
using namespace std;

typedef vector3D position;

class lattice
{
public:
    lattice(double _X, double _Y, double _Z);
    lattice(parameters& p, string outputFolder, string CXCL12fname, string  CXCL13fname); // Function that takes parameters and builds the grid
    string print();
    bool insideBorders(vector3D pos);
    bool insideLZ(vector3D &v, parameters& p);//return true if position inside LZ
    int X,Y,Z;
public:
    celltype celltypeat(position& v); // Sometimes we might be interested in celltype, others in the cell pointer
    //Danial: Last Update 21-10-18 01:45, checked --> clear
    int radius_x; //Radius X
    int cx ;  //X coordinate of Center
    int cy ;  //y coordinate of Center
    int cz;  //z coordinate of Center
    
    double GCborder ; //radius^2
    cell* cellat(position &v);
    celltype celltypeat(int x, int y, int z);// sometimes we might be interested in position, others in coordinates
    cell* cellat(int x, int y, int z);
    void putcellat(cell *c);
    vector3D getFreePosition(double minZperc = 0, double maxZperc = 1); // Find random free position in lattice
    pair<bool,vector3D>getFreeNeighbourPosition(cell* cell, vector<vector3D>& neighbours); //Calculate the next position of cell
    gauss_randomize thetas;
    vector3D get_position_mitosis(vector3D &pos);
    vector3D get_random_direction() ;
    vector <vector3D> getNeighbour_nn(vector3D &pos);
    vector <vector3D> getNeighbour_diag(vector3D &pos);
    vector3D get_nn_directed2( cell *tc);
    vector3D getfreeNeighbour_nn(vector3D &pos);
    vector3D getfreeNeighbour_diag(vector3D &pos);
    int is_at_border(vector3D pos);
    void removecellat(position& v); //Function to remove cell pointer from lattice  
public:
    vector < vector < vector < cell* > > > grid;

    double chemoat(molecules whichCXCL, int x, int y, int z); // Need cell size to calculate points in chemgrid
    double chemoat(molecules whichCXCL, position &p); // Need cell size to calculate points in chemgrid
private:
    chemokines3D chemo_grid;
public:
    void putAgFDCat(position& v, FDC *_fdc, double _AgAmount); // note: if ask out of bounds, don't put the Ag
    FDC* getFDCat(int x, int y, int z); //get amount of FDC at position xyz
    FDC* getFDCat(position& p); //get amount of FDC at position p
    double getAgat(int x, int y, int z); //get amount of Ag at position xyz
    double getAgat(position& p); //get amount of Ag at position p
    void removeAgAt(position& p, double RemoveAgAmount = 1); //For Ag consumption
    void AddTotalAmountAginLattice(double Agamount); //
private:
    vector < vector < vector < pair < FDC* , pair < vector < int > , double > > > > > FDC_Ag_grid; // Philippe: later, replace vector<int> by a type or an ID for an antigen
    double TotalAmountAginLattice; //For oserver analyzes
};

#endif // LATTICE_H
