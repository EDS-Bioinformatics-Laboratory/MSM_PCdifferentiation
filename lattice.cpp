#include "lattice.h"
#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include "vector3d.h"
using namespace std;
lattice::lattice(double _X, double _Y, double _Z) : chemo_grid(_X, _Y, _Z) {
  X = int(_X);
  Y = int(_Y);
  Z = int(_Z);
}
lattice::lattice(parameters& p, string outputFolder, string CXCL12fname, string CXCL13fname)
    : thetas(2), chemo_grid(p) {
  chemo_grid.loadchemokinesfromHyphasma(outputFolder + CXCL12fname, CXCL12);
  chemo_grid.loadchemokinesfromHyphasma(outputFolder + CXCL13fname, CXCL13);
  cerr << "Loaded chemokines from hyphasma to chemokinelatice" << endl;

  X = int(2 * p.par[radius] / p.par[dx] + 1);
  Y = int(2 * p.par[radius] / p.par[dx] + 1);
  Z = int(2 * p.par[radius] / p.par[dx] + 1);
  if ((X < 1) || (X > 1000))
    cerr << "Error! Diameter lattice (X):" << X << endl;
  if ((Y < 1) || (Y > 1000))
    cerr << "Error! Diameter lattice (Y):" << Y << endl;
  if ((Z < 1) || (Z > 1000))
    cerr << "Error! Diameter lattice (Z):" << Z << endl;
  grid.resize(X);
  for (int x = 0; x < X; ++x) {
    grid[x].resize(Y);
    for (int y = 0; y < Y; ++y) {
      grid[x][y].resize(Z, NULL);  // set initial pointers to NULL
    }
  }
  // Danial: Last Update 21-10-18 01:45, checked --> clear
  radius_x = int((X - 1) / 2);  // Radius X
  cx = int((X - 1) / 2);        // X coordinate of Center
  cy = int((Y - 1) / 2);        // y coordinate of Center
  cz = int((Z - 1) / 2);        // z coordinate of Center

  GCborder = double(radius * radius);  // radius^2
  pair<vector<int>, double> AgSequenceAmount({int(p.par[BCR_Length]), 0}, 0);
  pair<FDC*, pair<vector<int>, double> > FDC_Ag(
      0, AgSequenceAmount);  // Could this be a way to add it?
  FDC_Ag_grid.resize(X);
  for (int x = 0; x < X; ++x) {
    FDC_Ag_grid[x].resize(Y);
    for (int y = 0; y < Y; ++y) {
      FDC_Ag_grid[x][y].resize(
          Z, FDC_Ag);
    }
  }
}

// If cell position tries to get outside GC it wont be able to move as border
// nodes will be ocupied by border type.
celltype lattice::celltypeat(int x, int y, int z) {
  double cellcenteredpos =
      double(((x - cx) * (x - cx)) + ((y - cy) * (y - cy)) +
             ((z - cz) * (z - cz)) - (radius_x * radius_x));
  if (isnan(GCborder)) {
    cerr << "lattice:: Error in GCborder: " << GCborder << endl;
  }
  if (isnan(cellcenteredpos)) {
    cerr << "lattice:: Error in cellcenteredpos: " << cellcenteredpos << endl;
  }
  if ((x < 0) || (x >= X) || (y < 0) || (y >= Y) || (z < 0) || (z >= Z)) {
    { return border; }
  } else {
    if (cellcenteredpos > 0) {
      { return border; }
    } else {
      cell* cellthere = grid.at(x).at(y).at(z);
      if (cellthere) {
        return cellthere->cell_type;  // equivanent grid[][][]
      } else {
        return empty;
      }
    }
  }
}

celltype lattice::celltypeat(position& p) {
  return celltypeat(int(p.X), int(p.Y), int(p.Z));
}

cell* lattice::cellat(int x, int y, int z) {
  // If cell position tries to get outside GC it wont be able to move as border
  // nodes will be ocupied by border type.(cells only move there when empty is
  // returned!)
  if (celltypeat(x, y, z) == border) {
    cerr << " ERR: lattice::cellat: Accesses border node: " << x << " , " << y
         << " , " << z << endl;
    return NULL;
  }
  if (celltypeat(x, y, z) == empty) {
    cerr << " ERR: lattice::cellat: Accesses empty node: " << x << " , " << y
         << " , " << z << endl;
    return NULL;
  }
  return grid.at(x).at(y).at(z);
}

cell* lattice::cellat(position& p) {
  return cellat(int(p.X), int(p.Y), int(p.Z));
}

void lattice::putcellat(cell* c) {
  if ((c->position.X < 0) || (c->position.X >= X) || (c->position.Y < 0) ||
      (c->position.Y >= Y) || (c->position.Z < 0) || (c->position.Z >= Z)) {
    cerr << "ERROR lattice::putcellat: Cell out of borders" << c->printcell()
         << endl;
  }

  if (c == NULL) cerr << "Empty cell" << endl;

  if (celltypeat(c->position) != empty)  // Danial: not checked
  {
    cerr << "ERR: Occupied position= " << celltypeat(c->position) << endl;
  } else {
    grid.at(c->position.X).at(c->position.Y).at(c->position.Z) = c;
  }
}

void lattice::removecellat(position& v) {
  // Danial: changed this line for sphere case
  // Elena: using function instead of the same line (just a detail...doesnt
  // change anything!)
  if (celltypeat(v) == border) {
    cerr << "lattice::removing cell out of borders= " << v.X << " " << v.Y
         << " " << v.Z << endl;
  } else {
    grid.at(v.X).at(v.Y).at(v.Z) = NULL;
  }
}


// Convert cell position in lattice (int) to relative position (double) and call
// chemokine grid with it
double lattice::chemoat(molecules whichCXCL, int x, int y, int z) {
  vector3D p = {double(x), double(y), double(z)};
  return chemoat(whichCXCL, p);
}

double lattice::chemoat(molecules whichCXCL, position& p) {
  if (not(insideBorders(p))) {
    return 0;
  }
  return chemo_grid.concentrationat(whichCXCL, p);
}

// Used for cell initial random position Specify precentaje of Z axes you want
// to get random position from (DZ,LZ)
vector3D lattice::getFreePosition(double minZperc, double maxZperc) {
  if ((minZperc < 0.) || (minZperc > 1.) || (maxZperc < 0.) ||
      (maxZperc > 1.) || (minZperc >= maxZperc)) {
    cerr << "ERROR: lattice::getFreePosition(" << minZperc << "," << maxZperc
         << "), wrong fractions !" << endl;
  }

  vector3D freePosition = vector3D(-1, -1, -1);
  bool foundPlace = false;
  int ocupiedspaces = 0;

  while ((!foundPlace) && (ocupiedspaces < 10000000)) {
    ocupiedspaces = ocupiedspaces + 1;

    freePosition.X = random::randomInteger(X);  // extreme values are reached
    freePosition.Y = random::randomInteger(Y);
    // Danial: Dark zone or light zone?.
    freePosition.Z = random::randomInteger(1 + int(double(Z) * minZperc),
                                           int((double(Z)) * maxZperc));
    if (celltypeat(freePosition) == empty)  // DAnial: not chekced
    {
      foundPlace = true;
      ocupiedspaces = 0;
      return freePosition;
    }
  }

  if (ocupiedspaces >= 10000000) {
    cerr << "Can not find random empty position in 10,000,000 try." << endl;
    exit(1);
  }
  return vector3D(-1, -1, -1);
}

// Compute all neighbor Vector3D that are not border!

char shortCellType(celltype ct) {
  switch (ct) {
    case empty: {
      return 'E';
    }
    case FDCell: {
      return 'F';
    }
    case Stromalcell: {
      return 'S';
    }
    case TFHC: {
      return 'T';
    }
    case Centroblast: {
      return 'B';
    }
    case Centrocyte: {
      return 'C';
    }
    case Plasmacell: {
      return 'P';
    }
    case Memorycell: {
      return 'M';
    }
    case border: {
      return '|';
    }
    case cell_type_counter: {
      return '?';
    }
  }
  return '?';
}

string lattice::print() {
  stringstream res;
  res << "Grid of size: X: " << X << ", Y: " << Y << ", Z: " << Z << endl;
  for (int x = 0; x < X; ++x) {
    res << "x=" << x << endl;
    for (int y = 0; y < Y; ++y) {
      res << "y" << (y % 10) << " ";
      for (int z = 0; z < Z; ++z) {
        if (grid.at(x).at(y).at(z) == NULL) {
          res << "_";
        } else
          res << shortCellType(celltypeat(
              x, y,
              z));  // ???Does it make sense to print also cellat(x,y,z) here?..
      }
      res << endl;
    }
    res << endl;
  }
  return res.str();
}

bool lattice::insideLZ(vector3D& v, parameters& p) {
  bool inLZ = true;
  if ((v.X < 0) || (v.X >= X) || (v.Y < 0) || (v.Y >= Y) ||
      (v.Z < (p.par[zoneRatioGC] * Z)) || (v.Z >= Z)) {
    inLZ = false;
  }
  return inLZ;
}

void lattice::putAgFDCat(position& v, FDC* _fdc, double _AgAmount) {
  if (_fdc == NULL) cerr << "Error: Empty FDC" << endl;
  if (_AgAmount < 0 || _AgAmount > 100000)
    cerr << "Error: Ag out of range" << endl;
  FDC_Ag_grid.at(v.X).at(v.Y).at(v.Z).second.second += _AgAmount;
  FDC_Ag_grid.at(v.X).at(v.Y).at(v.Z).first = _fdc;
}

FDC* lattice::getFDCat(int x, int y, int z) {
  // If cell position tries to get outside GC it wont be able to move as border
  // nodes will be ocupied by border type.(cells only move there when empty is
  // returned!)
  if ((x < 0) || (x >= X) || (y < 0) || (y >= Y) || (z < 0) || (z >= Z)) {
    cerr << "lattice::getFDCat: ERROR: Cell out of borders" << x << " , " << y
         << " , " << z << endl;
    return NULL;
  }
  return FDC_Ag_grid.at(x).at(y).at(z).first;
}
FDC* lattice::getFDCat(position& v) {
  return getFDCat(int(v.X), int(v.Y), int(v.Z));
}

double lattice::getAgat(int x, int y, int z) {
  // If cell position tries to get outside GC it wont be able to move as border
  // nodes will be ocupied by border type.(cells only move there when empty is
  // returned!)
  if (not(insideBorders(vector3D(x, y, z)))) {
    cerr << "Error, get Ag out of borders." << x << " , " << y << " , " << z
         << endl;
    exit(1);
  }
  return FDC_Ag_grid.at(x).at(y).at(z).second.second;
}

double lattice::getAgat(position& v) { return getAgat(v.X, v.Y, v.Z); }

void lattice::removeAgAt(position& p, double RemoveAgAmount) {
  if (not(insideBorders(p))) {
    cerr << "Error, remove Ag out of borders" << p.print() << endl;
    exit(1);
  }
  FDC_Ag_grid.at(p.X).at(p.Y).at(p.Z).second.second -=
      RemoveAgAmount;  // Ag amount in lattice
}

void lattice::AddTotalAmountAginLattice(double Agamount) {
  TotalAmountAginLattice += Agamount;
}


bool lattice::insideBorders(vector3D pos) {
  return (celltypeat(pos) != border);
}
//#Recheck, danial:improvment
// This finds 6 near neighbours, 2 in each plan
vector<vector3D> lattice::getNeighbour_nn(vector3D& pos) {
  vector<vector3D> neighbours_nn;
  neighbours_nn.reserve(6);
  double tmp[3] = {pos.X, pos.Y, pos.Z};
  neighbours_nn.push_back(vector3D(tmp[0] - 1, tmp[1], tmp[2]));
  neighbours_nn.push_back(vector3D(tmp[0], tmp[1] - 1, tmp[2]));
  neighbours_nn.push_back(vector3D(tmp[0], tmp[1], tmp[2] - 1));
  neighbours_nn.push_back(vector3D(tmp[0] + 1, tmp[1], tmp[2]));
  neighbours_nn.push_back(vector3D(tmp[0], tmp[1] + 1, tmp[2]));
  neighbours_nn.push_back(vector3D(tmp[0], tmp[1], tmp[2] + 1));
  return neighbours_nn;
}

//#Recheck, danial:improvment
//This finds 12 diag neighbours
vector<vector3D> lattice::getNeighbour_diag(vector3D& pos) {
  vector<vector3D> neighbours_diag;
  neighbours_diag.reserve(12);
  double tmp[3] = {pos.X, pos.Y, pos.Z};
  neighbours_diag.push_back(vector3D(tmp[0] - 1, tmp[1], tmp[2] - 1));
  neighbours_diag.push_back(vector3D(tmp[0] - 1, tmp[1], tmp[2] + 1));
  neighbours_diag.push_back(vector3D(tmp[0] + 1, tmp[1], tmp[2] - 1));
  neighbours_diag.push_back(vector3D(tmp[0] + 1, tmp[1], tmp[2] + 1));

  neighbours_diag.push_back(vector3D(tmp[0] - 1, tmp[1] - 1, tmp[2]));
  neighbours_diag.push_back(vector3D(tmp[0] - 1, tmp[1] + 1, tmp[2]));
  neighbours_diag.push_back(vector3D(tmp[0] + 1, tmp[1] - 1, tmp[2]));
  neighbours_diag.push_back(vector3D(tmp[0] + 1, tmp[1] + 1, tmp[2]));

  neighbours_diag.push_back(vector3D(tmp[0], tmp[1] - 1, tmp[2] - 1));
  neighbours_diag.push_back(vector3D(tmp[0], tmp[1] - 1, tmp[2] + 1));
  neighbours_diag.push_back(vector3D(tmp[0], tmp[1] + 1, tmp[2] - 1));
  neighbours_diag.push_back(vector3D(tmp[0], tmp[1] + 1, tmp[2] + 1));
  return neighbours_diag;
}
//#Recheck, danial:improvment
vector3D lattice::getfreeNeighbour_nn(vector3D& pos) {
  vector<vector3D> neighbours_nn;
  neighbours_nn.reserve(6);

  vector<vector3D> freeneighbours_nn;
  vector3D SelectdNeighbour(-1, -1, -1);

  neighbours_nn = getNeighbour_nn(pos);

  for (unsigned int i = 0; i < neighbours_nn.size(); i++) {
    if (insideBorders(neighbours_nn[i])) {
      if (celltypeat(neighbours_nn[i]) == empty) {
        freeneighbours_nn.push_back(neighbours_nn[i]);
      }
    }
  }
  if (freeneighbours_nn.size() > 0) {
    int x = random::randomInteger(0, int(freeneighbours_nn.size()));
    SelectdNeighbour.X = freeneighbours_nn[x].X;
    SelectdNeighbour.Y = freeneighbours_nn[x].Y;
    SelectdNeighbour.Z = freeneighbours_nn[x].Z;
  }
  return SelectdNeighbour;
}
//#Recheck, danial:improvment
vector3D lattice::getfreeNeighbour_diag(vector3D& pos) {
  vector<vector3D> neighbours_diag;
  neighbours_diag.reserve(12);
  vector<vector3D> freeneighbours_diag;
  vector3D SelectdNeighbour(-1, -1, -1);
  neighbours_diag = getNeighbour_diag(pos);

  for (unsigned int i = 0; i < neighbours_diag.size(); i++) {
    if (insideBorders(neighbours_diag[i])) {
      if (celltypeat(neighbours_diag[i]) == empty) {
        freeneighbours_diag.push_back(neighbours_diag[i]);
      }
    }
  }
  if (freeneighbours_diag.size() > 0) {
    int x = random::randomInteger(0, int(freeneighbours_diag.size()));
    SelectdNeighbour.X = freeneighbours_diag[x].X;
    SelectdNeighbour.Y = freeneighbours_diag[x].Y;
    SelectdNeighbour.Z = freeneighbours_diag[x].Z;
  }
  return SelectdNeighbour;
}
//#Recheck, danial:improvment
vector3D lattice::get_position_mitosis(vector3D& pos) {
  vector3D place(-1, -1, -1);
  place = getfreeNeighbour_nn(pos);
  if (place.X == -1 || place.Y == -1 || place.Z == -1) {
    place = getfreeNeighbour_diag(pos);
    if (place.X == -1 || place.Y == -1 || place.Z == -1) {
      return vector3D(-1, -1, -1);
    } else {
      return vector3D(place);
    }
  }
  return vector3D(place);
}

vector3D lattice::get_random_direction() {
  vector3D random_direction(-1.0, -1.0, -1.0);
  double phi = random::randomDouble(2.0 * 3.141592654);
  random_direction.X = double(cos(phi));
  random_direction.Y = double(sin(phi));
  double theta = random::randomDouble(3.141592654);
  random_direction.X *= sin(theta);
  random_direction.Y *= sin(theta);
  random_direction.Z = cos(theta);
  return (random_direction);
}
//#Recheck, danial:improvment
vector3D lattice::get_nn_directed2(cell* c1) {
  vector3D selected(-1, -1, -1);
  vector<int> index;
  vector3D normalized_polarity = c1->polarity.getNormalizedVector();

  double tmp = max(abs(normalized_polarity.X), abs(normalized_polarity.Y));
  tmp = max(abs(tmp), abs(normalized_polarity.Z));

  if (fabs(tmp - normalized_polarity.X) < 1e-6) {
    index.push_back(3);
  }
  if (fabs(tmp - normalized_polarity.Y) < 1e-6) {
    index.push_back(4);
  }
  if (fabs(tmp - normalized_polarity.Z) < 1e-6) {
    index.push_back(5);
  }
  if (fabs(tmp + normalized_polarity.X) < 1e-6) {
    index.push_back(0);
  }
  if (fabs(tmp + normalized_polarity.Y) < 1e-6) {
    index.push_back(1);
  }
  if (fabs(tmp + normalized_polarity.Z) < 1e-6) {
    index.push_back(2);
  }

  vector<vector3D> neighbours;
  neighbours.reserve(6);
  neighbours = getNeighbour_nn(c1->position);
  vector<vector3D> accepted_neighbour;
  for (unsigned int i = 0; i < index.size(); i++) {
    vector3D tmp_neighbour = {neighbours[index[i]].X, neighbours[index[i]].Y,
                              neighbours[index[i]].Z};
    if (insideBorders(tmp_neighbour)) {
      accepted_neighbour.push_back(neighbours[index[i]]);
    }
  }
  if (accepted_neighbour.size() > 0) {
    int x = fabs(accepted_neighbour.size());
    selected = accepted_neighbour[random::randomInteger(x)];
  }
  return selected;
}

int lattice::is_at_border(vector3D pos) {
  // In this function we do not consider the polarity of cell when it is at the
  // border.
  vector<vector3D> neighbours = getNeighbour_nn(pos);
  for (int i = 0; i < neighbours.size(); i++) {
    if (celltypeat(neighbours.at(i)) == border) {
      return 1;
    }
  }
  return 0;
}

