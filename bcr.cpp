#include "bcr.h"
#include "random.h"
#include<cmath>
#include <sstream>
using namespace std;
vector <vector <int>> Shape_Space_Seeder;

// Constructor of BCR for the shape-space method
BCR::BCR(parameters& p)
{   //BCR constructor, sets the number of mutations from germline to 0, takes the size of BCR from parameter file (default: 4) and sets the BCR receptor using a random seed from a pool of randomly produced seeds (default: 100).
    nMutFromGermline = 0.;
    BCReceptor.resize(p.par[BCR_Length], 0);
    BCReceptor=Shape_Space_Seeder[random::randomInteger(p.par[BCR_pool])];
}

void initialize_Seeds (parameters &p)
{
    // This function creates a pool of random seeds (default: 100).
    Shape_Space_Seeder.resize(p.par[BCR_pool],vector<int>(p.par[BCR_Length]));

    for (int j=0;j<p.par[BCR_pool];j++)
    {
        for (int i =0; i< p.par[BCR_Length];i++)
        {
            Shape_Space_Seeder[j][i]=random::randomInteger(p.par[Nmax]+1); // default: Nmax is 9 --> each digit is between 0-9 so we use random(10) to include 9 itself.
        }
//        for ( int i = 0, inc = 1; i < p.par[BCR_Length]; ++inc, i += inc )
//        {
//        cout << i << endl;
//            Shape_Space_Seeder[j][i]=5;
//        }
//        for ( int i = 1, inc = 1; i < p.par[BCR_Length]; ++inc, i += inc )
//        {
//        cout << i << endl;
//            Shape_Space_Seeder[j][i]=4;
//        }
    }
}

//#Recheck danial: This can get improved
long double BCR::mutateBCR(parameters& p)
{
    long double aff1= getMyAffinity4Ag(p);
    if(BCReceptor.size() != 4) cerr << "Error:MutateBCR: BCR size is: "<< BCReceptor.size()<< endl;
    if(random::randomDouble(1) < pMut)
    {
        // This "step" is the size of affinity change (delta_affinity). setting it as it is now will improve the average affinity of output cells by making change in affinity smaller for cells with high affinity.
        
        double step=0;
        double max_step=2.5;
        double aff = double(int(getMyAffinity4Ag(p) * 10))/ 10;
        step= (-1*max_step*aff)+max_step;     //#temporary Danial: To make the probability of mutation independent of the affinity and keeping the average affinity of output cells still high, it is possible to change the mutation step size
        nMutFromGermline += 1; //Increase numer of mutations from the beginning founder
        short int randomL; //random location
        short int randomC; //random change
        randomL = random::randomInteger(4);
        randomC = random::randomInteger(2);
        while ((BCReceptor[randomL]==9 && randomC==1)|| (BCReceptor[randomL]== 0 && randomC==0))
        {
            randomL = random::randomInteger(4);
            randomC = random::randomInteger(2);
        }
        if (randomC==1)
        {BCReceptor[randomL] += step;}
        if (randomC==0)
        {BCReceptor[randomL] -= step;}
    }
    long double aff2= getMyAffinity4Ag(p);
    return double(aff2-aff1); // Returns the change in affinity due to the mutation
}
//This function calculates affinity using shape-space concept by computing the distance between Ag and bcr
double BCR::getMyAffinity4Ag (parameters &p)
{
    double dist = 0.;
    for(unsigned int i = 0; i < BCReceptor.size(); i++)
    {
       dist += fabs(double (3. - BCReceptor[i])); // Philippe: now the target antigen is 3333, modify later
    }
    double dist2=dist*dist;
    return exp(-1. * dist2 / (2.8 * 2.8)); //Danial: #Recheck
}

string BCR::print_BCR(){
    stringstream res;
    int L = (int) BCReceptor.size();
    for(int i = 0; i < L; ++i)
    {
        res<<BCReceptor[i];
    }
    return res.str();
}




