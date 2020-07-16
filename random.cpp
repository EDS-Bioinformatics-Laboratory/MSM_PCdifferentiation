#ifdef _WIN32
#include <windows.h>
#include <time.h>
#elif __APPLE__
#elif __linux__
#endif
#include "random.h"
#include <cmath>
#include <random>
#include <iostream>
#include <algorithm>
#include "output.h"

//Danial: There is a good space for improvement in this file. All can be summerized in one single class and there are parts that can be replaced by built-in or inline functions,... .

using namespace std;
const double gauss_randomize::Ny = 10;
int dtphase_frequency[4][21]={{0}};

//#Recheck
double random::cell_cycle_time (double i, int current_cycle)
{
    double phaselength;
    double arg;
    double thisphase = i;
    double width = 1. * double(thisphase); //"fraction_of_phase_width"=1 later change with param, from hyphasma
    phaselength = 3.0 * thisphase;
    while (phaselength <= 0. || phaselength >= 2. * thisphase) {
        // This samples from a normal distribution with mean thisphase and width width
        arg = (2. * random::randomDouble(1.0) - 1.);
        phaselength = thisphase + sqrt(2.) * width * inverse_erf(arg);
    }
    int dt_index = int (phaselength / (2.0 * thisphase / double (21. - 1.)) + 0.5);
    dtphase_frequency[current_cycle][dt_index]++;
    return phaselength;
    
}

int random::randomInteger(int a){
    return int (double (rand()) * a / (double (RAND_MAX) + double (1)));
}


int random::randomInteger(int a,int b){
    return  (a + (rand() / ((RAND_MAX /(b - a))+1)));
}



double random::randomDouble(double bis){
    return double (rand()) * bis / (double (RAND_MAX) + double (1)); //Danial: Added to be sure we use same distribution.
}

gauss_randomize::gauss_randomize() {
    arraydim = 100;
    field.reserve(100);
}

gauss_randomize::gauss_randomize(short dataset) {
    
    if (dataset == 1) { gauss_initialize(1.05,1.05,181,0,3.141592654); } else if (dataset == 2) {
        cyster_initialize();
    }
}

double gauss_randomize::gauss(double &x, double &x0, double &width) {
    return exp(-1. * (x - x0) * (x - x0) / (width * width));
}
void gauss_randomize::gauss_initialize(double x0, double width, int Nx, double xmin, double xmax) {
    /* Values between xmin and xmax at precision (xmax-xmin)/(Nx-1) are saved in an array.
     * [N-1 in the denominator ensures that there are really Nx possible values including
     * both interval limits [xmin,xmax]]
     * The frequency of occurrence corresponds to a Gaussian distribution
     * centered at x0 with width "width".
     * The dimension of array is a result of the minimum ymin of the Gaussian in the
     * considered interval [xmin,xmax]. The vertical resolution is then given by
     * ymin/Ny, where Ny=10 in standard situation.
     */
    // Calculate the x-interval
    double dx = (xmax - xmin) / (double (Nx - 1));
    // Find the minimum of the Gaussian (ymin)
    double ymin = 2.;
    for (int i = 0; i < Nx; i++) {
        double x = xmin + double (i) * dx;
        double y = gauss(x,x0,width);
        if (y < ymin) { ymin = y; }
    }
    if (ymin < 0.001) { ymin = 0.001; }
    double dy = ymin / Ny;
    // Calculate the total number of necessary entries in the array
    arraydim = 0;
    for (int i = 0; i < Nx; i++) {
        double x = xmin + double (i) * dx;
        arraydim += int (gauss(x,x0,width) / dy + 0.5);
    }
    int ind = 0;
    for (int i = 0; i < Nx; i++) {
        double x = xmin + double (i) * dx;
        int n = int (gauss(x,x0,width) / dy + 0.5);
        for (int j = 0; j < n; j++) {
            field.push_back(x);
            ++ind;
        }
    }
}
double gauss_randomize::cyster07angle_wt(int angle) {
    switch (angle) {
        case 5:
            return 23;
            break;
        case 15:
            return 61;
            break;
            
        case 25:
            return 85;
            break;
            
        case 35:
            return 97;
            break;
            
        case 45:
            return 100;
            break;
            
        case 55:
            return 82;
            break;
            
        case 65:
            return 79;
            break;
            
        case 75:
            return 62;
            break;
            
        case 85:
            return 53;
            break;
            
        case 95:
            return 44;
            break;
            
        case 105:
            return 40;
            break;
            
        case 115:
            return 35;
            break;
            
        case 125:
            return 30;
            break;
            
        case 135:
            return 25;
            break;
            
        case 145:
            return 20;
            break;
            
        case 155:
            return 17;
            break;
            
        case 165:
            return 9;
            break;
            
        case 175:
            return 4;
            break;
    }
    return -1;
}
void gauss_randomize::cyster_initialize() {
    /* Does the same as the Gauss-constructor but using a fixed dataset.
     * Different datasets may be stored here and may be called by the variable dataset.
     */
    // Use Allen 2007 Figure S5 C (green) for the turning angle distribution
    // define the number of x-values
    int Nx = 18;
    // define the minimum x-value
    int xmin = 5;
    // Calculate the x-interval
    int dx = 10;
    // Find the minimum occuring value
    double ymin = 4.;
    // Define the vertical resolution (Ny is a local constant double)
    double dy = ymin / Ny;
    // Calculate the total number of necessary entries in the array
    arraydim = 0;
    for (int i = 0; i < Nx; i++) {
        int x = xmin + i * dx;
        arraydim += int (cyster07angle_wt(x) / dy + 0.5);
    }
    // Define the array
    field.reserve(arraydim) ;
    // Fill in the entries
    int ind = 0;
    for (int i = 0; i < Nx; i++) {
        int x = xmin + i * dx;
        int n = int (cyster07angle_wt(x) / dy + 0.5);
        double xrad = double (x) * 3.141592654 / 180.;
        for (int j = 0; j < n; j++) {
            field.push_back(xrad);
            ++ind;
        }
    }
    cout << "Allen et al. 2007 turning-angle-distribution for WT initialised.\n";
}
gauss_randomize::~gauss_randomize() {
    field.clear();
}
double gauss_randomize::get_distribution_value() {

    return field[random::randomInteger(arraydim)];
}

double inverse_erf(double x) {
    const double pi = 3.141592654;
    const double a = 8. * (pi - 3) / (3. * pi * (4. - pi));
    double z = 0;
    double kl = (2. / (pi * a) + 0.5 * log(1. - x * x));
    //  cout<<x<<"-->"<<kl<<">"<<sqrt(kl-log(1.-x*x)/a)<<"\n";
    z = sqrt(sqrt(kl * kl - log(1. - x * x) / a) - kl);
    if (x < 0) { z *= -1.; }
    return z;
}
