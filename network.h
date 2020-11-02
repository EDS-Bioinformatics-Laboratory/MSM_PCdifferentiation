#ifndef NETWORK_H
#define NETWORK_H

#include "odesolver.h"
#include "parameters.h"
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace std;
enum DLBCL_pars {Healthy_pars,M2,M8,M2M8,NbModels_pars};
enum DLBCL_dervs {Healthy_dervs,M1,M3A,M3B,M4,NbModels_dervs};

struct network:
        public odesolver {
public:
            network();
            enum{init_p, init_b, init_r,bcr,cd40,Mu_p, Mu_b, Mu_r, sigma_p, sigma_b, sigma_r, l_p,l_b,l_r,k_p, k_b, k_r, NbParameters};
            enum{b, r, p, NbVariables};
            DLBCL_pars dlbcl_pars;
            DLBCL_dervs dlbcl_dervs;
            virtual void derivatives(const vector<double> &x, vector<double> &dxdt, const double t);
            void initialise();
            void setBaseParameters(parameters& p);
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
