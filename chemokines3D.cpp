#include "chemokines3d.h"
#include <random>
//#Recheck, danial: this file can be improved significantley
chemokines3D::chemokines3D(double _X, double _Y, double _Z)
{
    X=_X;
    Y=_Y;
    Z=_Z;
}
chemokines3D::chemokines3D(parameters& p )
{
    X = (2 * p.par[radius] / p.par[chemo_dx]) +1;
    Y = (2 * p.par[radius] / p.par[chemo_dx]) +1;
    Z = (2 * p.par[radius] / p.par[chemo_dx]) +1;
    if(( X < 1 ) || (X > 100000000)) cerr << "Error! Diameter chemokine lattice (X):" << X << endl;
    chemo_grid.resize(X);
    for(int x = 0; x < X; ++x){
        if(( Y < 1 ) || (Y > 100000000)) cerr << "Error! Diameter chemokine lattice (Y):" << Y << endl;
        chemo_grid[x].resize(Y);
        for(int y = 0; y < Y; ++y){
            if(( Z < 1 ) || (Z > 100000000)) cerr << "Error! Diameter chemokine lattice (Z):" << Z << endl;
            chemo_grid[x][y].resize(Z);
            for(int z = 0; z < Z; ++z){
                chemo_grid[x][y][z].resize(N_molecules, NAN);
            }
        }
    }
}

void chemokines3D::loadchemokinesfromHyphasma(string fname, molecules whichCXCL)
{
    cerr << "      Load " << fname << " from file ... "<< endl;
    ifstream ff(fname);
    if(!ff) {
        cerr << "ERROR: couldn't find the signal file: " << fname << endl;
        exit(-1);
    }
    string tmp = "";
    for (short a = 0; a < N_PRE_WORDS ; a++) {
       ff >> tmp; //string
    }
    long LN_Index,x,y,z;
    double concentration,d5,d6,d7,d8,d9,d9a;
    int count = 0; //Counter for while loop
    while(!ff.eof()) {
       ff >> LN_Index >> x >> y >> z>>concentration >> d5 >> d6 >> d7 >> d8 >> d9>> d9a; //LNindex x y z conc gradient chemoweight %chemoweight  grad directeion(x,y,z)
       if (ff.bad())
       { cerr<<"Loading chemokines "<< whichCXCL<<" failed!"<<endl;
         exit(1); }
       else {
           vector3D chemo_pos = vector3D((int)x,(int)y,(int)z);
           putchemokines3D(whichCXCL,chemo_pos,concentration);
       }

       if(count >= 1000000) cerr << "infinite loop loading chemokines" << endl;
    }
    ff.close();
    cerr << " done.\n";
    }

void  chemokines3D::putchemokines3D(molecules whichCXCL,vector3D& pos,double concentration_ )
{

    if(( pos.X < 0 ) || ( pos.X > X )||(pos.Y < 0) ||(pos.Y > Y)||( pos.Z < 0 )||( pos.Z > Z ))
        cerr<< "ERROR: Putting chemokine out of lattice!"<<endl;
    if (isnan(concentration_)) cerr<<"Cant put chemokine concentration = NAN"<<endl;
    
    chemo_grid.at(pos.X).at(pos.Y).at(pos.Z).at(whichCXCL) = concentration_;

}

double chemokines3D::concentrationat(molecules whichCXCL, int x, int y, int z)
{

    return chemo_grid[x][y][z][whichCXCL]; //#Rechcek danial
}

double chemokines3D::concentrationat(molecules whichCXCL, position &p)
{
    return concentrationat(whichCXCL,p.X, p.Y, p.Z);
}

string chemokines3D::printchemokines3D(){

    stringstream res;
    res << "Chem_Grid of size: X: " << X << ", Y: " << Y << ", Z: " << Z << endl;
    res << "Chemokine X,Y,Z position"<< ' ' << "CXCL12 Concentration"<< ' ' << "CXCL13 Concentration"<<endl;
    for(int x = 0; x < X; ++x)
    {
        for(int y = 0; y < Y; ++y)
        {
            for(int z = 0; z < Z; ++z)
            {
                 res << x <<" "<< y<< " " << z <<" "<< chemo_grid.at(x).at(y).at(z).at(CXCL12)<< " " << chemo_grid.at(x).at(y).at(z).at(CXCL13)<<endl;
            }
        }
    }
    return res.str();
}

void chemokines3D::writechemokines(string fname)
{
    ofstream myfile;
    myfile.open(fname);
    if (!myfile) cerr << "ERROR! My parameter file empty";
    myfile << printchemokines3D();
    myfile.close();
}
