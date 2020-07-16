#ifndef i_random
#define i_random
#include <random>
#include <vector>
using namespace std;

struct random
{
//    static std::mt19937* generator;
//    static bool initialized;
//    static void initializeRandom(unsigned int init_seed=0);
    static int randomInteger (int); // careful, borders are reached
    static int randomInteger(int a,int b);
    static double randomDouble (double);
    static double cell_cycle_time (double,int);
//    static double randomGaussian(double, double);
//    static double boundedRandomGaussian(double mean , double SD, double vmin, double vmax);
//    template<typename T>
//
////    static void shuffle(vector<T> v){
////        if(!initialized) initializeRandom();
////        std::shuffle(v.begin(), v.end(), *generator);
////    }
//
//    void shuffle(vector<double>);
//    void shuffle(vector<int>);
};
double inverse_erf(double x);

class gauss_randomize {
public:
    gauss_randomize();
    gauss_randomize(short dataset);
    ~gauss_randomize();
    double get_distribution_value();
private:
    void gauss_initialize(double, double, int, double, double);
    double gauss(double&,double&,double&);
    void cyster_initialize();
    double cyster07angle_wt(int angle);
    vector <double> field;
    int arraydim;
    
    ////§§§ Philippe 21-03-2017
    static const double Ny;
};
/*class permutations {
  public:
   permutations(int dim);
   vector<int> get_rand_permutation();
  private:
   int n_permut;
   int dimension;
   vector< vector <int> > savedPermutations;
};*/

#endif
