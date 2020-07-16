#include "vector3d.h"
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iostream>
#include <numeric>
#include <iostream>
#include <vector>
#include <functional>

vector3D::vector3D(double _X, double _Y, double _Z) {
  X = _X;
  Y = _Y;
  Z = _Z;
}

// Redefenition of operator +
vector3D vector3D::operator+(vector3D Addend_vector) {
  vector3D Sum = vector3D(0., 0., 0.);
  Sum.X = X + Addend_vector.X;
  Sum.Y = Y + Addend_vector.Y;
  Sum.Z = Z + Addend_vector.Z;
  return Sum;
}
// Redefenition of operator -
vector3D vector3D::operator-(vector3D Subtrahend_vector) {
  vector3D Difference = vector3D(0., 0., 0.);
  Difference.X = X - Subtrahend_vector.X;
  Difference.Y = Y - Subtrahend_vector.Y;
  Difference.Z = Z - Subtrahend_vector.Z;
  return Difference;
}
// Redefenition of operator *
vector3D vector3D::operator*(double Scalar) {
  vector3D Product = vector3D(1., 1., 1.);
  Product.X = X * Scalar;
  Product.Y = Y * Scalar;
  Product.Z = Z * Scalar;
  return Product;
}

// Redefinition of operator == (Ture/False) if X,Y,Z coordinates are equal.
bool vector3D::operator==(vector3D &v) {
  double error = 1e-10;
  return ((abs(X - v.X) <= error) && (abs(Y - v.Y) <= error) &&
          (abs(Z - v.Z) <= error));
}

//#Recheck danial: why?
bool vector3D::compare(double a) {
  double error = 1e-10;
  return ((abs(X - a) <= error) && (abs(Y - a) <= error) &&
          (abs(Z - a) <= error));
}
// Redefinition of operator [] (access by index).
double vector3D::operator[](int i) {
  if (i == 0) {
    return X;
  } else if (i == 1) {
    return Y;
  } else if (i == 2) {
    return Z;
  } else {
    cout << "Incorrect index for vector3D" << endl;
    exit(1);
  }
}

// This function returns the norm of a vector
double vector3D::getNorm() {
  double norm = sqrt((X * X) + (Y * Y) + (Z * Z));
  return norm;
}

// This function normalizes the vector
void vector3D::NormalizeYourself() {
  vector3D normalize = vector3D(X / getNorm(), Y / getNorm(), Z / getNorm());
  X = normalize.X;
  Y = normalize.Y;
  Z = normalize.Z;
}

// Create a new normalized vector
vector3D vector3D::getNormalizedVector() {
  double Norm_of_vector = getNorm();
  if (fabs(Norm_of_vector) > 1e-9) {
    X /= Norm_of_vector;
    Y /= Norm_of_vector;
    Z /= Norm_of_vector;
  } else {
    X = 0;
    Y = 0.;
    Z = 0.;
  }
  vector3D normalized_vector = vector3D(X, Y, Z);
  return normalized_vector;
}

// Calculate scalar product of 2 vectors
double getScalarproduct(vector3D &u, vector3D &v) {
  double scalar_product = (u.X * v.X) + (u.Y * v.Y) + (u.Z * v.Z);
  return scalar_product;
  // It is possible to change this with built-in inner_product function
  //    vector <double> a(3);
  //    vector <double> b(3);
  //    for (int i=0;i<3;i++)
  //    {
  //        a[i]=u[i];
  //        b[i]=v[i];
  //    }
  //
  //    return double (inner_product(a.begin(), a.end(), b.begin(), 0.0));
}

// This function calculates distance between 2 vectors
double getdist(vector3D &u, vector3D &v) {
  double dist =
      sqrt(pow((u.X - v.X), 2) + pow((u.Y - v.Y), 2) + pow((u.Z - v.Z), 2));
  return dist;
}

// This function computes projections of neighbors
//#Recheck danial: What is this for ?!?
double getProjection(vector3D &u, vector3D &v) {
  double projection =
      getScalarproduct(u, v) /
      (u.getNorm() * v.getNorm());  //  getScalarproduct() is outside of
                                    //  vector3D class so it is called as a
                                    //  function; getNorm() is a property of
                                    //  vector (to access .getNorm())
  return projection;
}

string vector3D::print() {
  stringstream res;
  res << " X: " << X << ", Y: " << Y << ", Z: " << Z << " ";
  return res.str();
}

