#ifndef BCR_H
#define BCR_H
#include <vector>
#include "parameters.h"
#include <string>
void initialize_Seeds (parameters &p);
class BCR
{
    
    //This class includes the implementation of B-cell receptor
public:
    
    //Shape-space method
    BCR(parameters &p);
    vector<int> BCReceptor;
    double pMut; // Probability of mutation
    int    nMutFromGermline; // Number of mutations from the beginning founder cell
    long double mutateBCR(parameters &p);
    double getMyAffinity4Ag(parameters& p);
    string print_BCR();

};

#endif // BCR_H
