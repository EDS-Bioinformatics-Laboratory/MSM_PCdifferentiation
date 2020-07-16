#ifndef CHEMOKINES3D_H
#define CHEMOKINES3D_H
#include "vector3d.h"
#include <iostream>
#include "parameters.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

enum molecules {CXCL12, CXCL13, N_molecules};
typedef  vector3D position;

class lattice;

class chemokines3D
{
public:
    //chemokines3D();
    chemokines3D(double _X, double _Y, double _Z);
    chemokines3D(parameters& p); // Function that takes parameters and builds the grid
    void loadchemokinesfromHyphasma(string fname, molecules whichCXCL);
    void putchemokines3D(molecules whichCXCL, vector3D &pos, double concentration_ );

    // BE CAREFUL, these positions are not the same as the cell lattice positions !!!
    double concentrationat(molecules whichCXCL, position& p); // Sometimes we might be interested in celltype, others in the cell pointer
    double concentrationat(molecules whichCXCL, int x, int y, int z);// sometimes we might be interested in position, others in coordinates
    // Now talking between 0 and 1 for each dimension. Constraint: both cell and chemokine lattice should cover the same real space in each diomension
    double concentrationatnormalizedposition(molecules whichCXCL, position &p); // Sometimes we might be interested in celltype, others in the cell pointer

    vector < vector < vector < vector < double> > > > chemo_grid;
    int X,Y,Z;
    static const int N_PRE_WORDS = 20;
    string printchemokines3D();
    void  writechemokines(string fname);
};

#endif // CHEMOKINES3D_H
