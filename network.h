#ifndef NETWORK_H
#define NETWORK_H

#include "odesolver.h"
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace std;

struct network:
        public odesolver {
public:
            network();
            enum{init_p, init_b, init_r,bcr,cd40,Mu_p, Mu_b, Mu_r, sigma_p, sigma_b, sigma_r, l_p,l_b,l_r,k_p, k_b, k_r, NbParameters};
            enum{b, r, p, NbVariables};

            virtual void derivatives(const vector<double> &x, vector<double> &dxdt, const double t);
            void initialise();
            void setBaseParameters();
            void setDinamicParameters(double p, double b, double r, double cd40_,double bcr_);
            void action(string name, double parameter){
                if(!name.compare("wash")){
                    if((parameter > 1.0) || (parameter < 0)) {cerr << "ERR: odesolverMinLatent::action(" << name << ", " << parameter << "), wrong parameter value\n"; return;}
                    return;
                }
            }

            double hillA(double V, double K){
                if(V < 0) return 0;
                return V*V / (V*V + K*K);
            }
            double hillI(double V, double K){
                if(V < 0) return 1;
                return K*K / (V*V + K*K);
            }
            double sq(double v1){
                return v1*v1;
            }
};

#endif // NETWORK_H
