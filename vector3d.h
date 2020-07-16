#ifndef VECTOR3D_H
#define VECTOR3D_H
#include <string>
#include <vector>

using namespace std;

class vector3D {
 public:
  vector3D(double _X, double _Y, double _Z);
  double X, Y, Z;
  vector3D operator+(vector3D v);
  vector3D operator-(vector3D v);
  vector3D operator*(double v);
  double operator[](int i);
  bool operator==(vector3D& v);
  bool compare(double v);

  // Calculate norm of vector
  double getNorm();
  // Normalize vector that changes XYZ
  void NormalizeYourself();
  // Calculate normalized vector
  vector3D getNormalizedVector();
  string print();
};

double getdist(vector3D& u, vector3D& v);
double getScalarproduct(
    vector3D& u,
    vector3D& v);  // Calculate scalar product of 2 vectors (u and v)
double getProjection(vector3D& u,
                     vector3D& v);  // Calculate projection of 2 vectors

#endif  // VECTOR3D_H
