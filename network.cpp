#include "network.h"

network::network(): odesolver(NbVariables, NbParameters)
{
        name = string("ODEmodel_PCdifferentiation");
        dt = 0.002; // initial time step -> then it is adaptive
        print_all_secs = 1200; //every how many seconds it is plotting

        // Name of variables

        names[b] = string("BCL6");
        names[r] = string("IRF4");
        names[p] = string("PRDM1");

        // Name of parameters

        paramNames[init_p] = "Initial levels PRDM1";
        paramNames[init_b] = "Initial levels BCL6";
        paramNames[init_r] = "Initial levels IRF4";
        paramNames[Mu_p] = "Basal transcription rate PRDM1";
        paramNames[Mu_b] = "Basal transcription rate BCL6";
        paramNames[Mu_r] = "Basal transcription rate IRF4";
        paramNames[sigma_p] = "Maximum induced transcription rate PRDM1";
        paramNames[sigma_b] = "Maximum induced transcription rate BCL6";
        paramNames[sigma_r] = "Maximum induced transcription rate IRF4";
        paramNames[l_p] = "Degradation transcription rate PRDM1";
        paramNames[l_b] = "Degradation transcription rate BCL6";
        paramNames[l_r] = "Degradation transcription rate IRF4";
        paramNames[k_p] = "Dissociation constant PRDM1";
        paramNames[k_b] = "Dissociation constant BCL6";
        paramNames[k_r] = "Dissociation constant IRF4";
        paramNames[bcr] = " CC-Ag Interaction";
        paramNames[cd40] = "CC-TFHC Interaction";


        backSimulated.clear();

    }

//Elena: set enum depending on parameter file.

//Elena: set network parameters depending on DLBCL enum
    void network::setBaseParameters(parameters& p){
        params.clear();     // to make sure they are all put to zero
        params.resize(NbParameters, 0.0);
        // NOTE: Parameters in units of time in martinez are per 4 hrs. I convert them to per hr. (value / 4)
        // Elena: Note this is only set once when a Bcell is created.
        
        DLBCL_pars dlbcl_pars_type = (DLBCL_pars) p.par[dlbcl_parameters]; //Elena: take DLBCL_type from parameter file for parameters
        dlbcl_pars = dlbcl_pars_type;
        
        DLBCL_dervs dlbcl_dervs_type = (DLBCL_dervs) p.par[dlbcl_derivatives]; //Elena: take DLBCL_type from parameter file for derivatives
        dlbcl_dervs = dlbcl_dervs_type;
        
        if(dlbcl_pars_type  ==NAN ||dlbcl_pars_type < 0 || dlbcl_pars_type > (DLBCL_pars) 10 )
            cerr<<"Error: network::setBaseParameters dlbcl_pars_types are: "<<dlbcl_pars_type <<"; Check parameters"<<endl;
        if( dlbcl_dervs_type==NAN || dlbcl_dervs_type < 0 ||  dlbcl_dervs_type > (DLBCL_dervs) 10)
            cerr<<"Error: network::setBaseParameters dlbcl_dervs_types are: "<<dlbcl_dervs_type <<"; Check parameters"<<endl;
        
        switch (dlbcl_pars)
        {
            case Healthy_pars:
            {
                params[Mu_b] = 2 ;//M/t //Elena: DLBCL: Uncomment to run reference model
                params[Mu_r] = 0.1 ;//M/t //Elena: DLBCL: Uncomment to run reference model
                params[sigma_r] = 2.6 ;//M/t //Elena: DLBCL: Uncomment to run reference model
                break;
            }
            case M2:
            {
                params[Mu_b] = 2 * 10 ;//M/t  //Elena: DLBCL: Uncomment to run M2 model
                params[Mu_r] = 0.1 ;//M/t //Elena: DLBCL: Uncomment to run reference model
                params[sigma_r] = 2.6 ;//M/t //Elena: DLBCL: Uncomment to run reference model
                break;
            }
            case M8:
            {
                params[Mu_b] = 2 ;//M/t //Elena: DLBCL: Uncomment to run reference model
                params[Mu_r] = 0.1 *2;//M/t //Elena: DLBCL: Uncomment to run M8 model (2-fold)
                params[sigma_r] = 2.6 *2;//M/t //no units//Elena: DLBCL: Uncomment to run M8 model (2- fold) ;
                break;
            }
            case M2M8:
            {
                params[Mu_b] = 2 * 10 ;//M/t  //Elena: DLBCL: Uncomment to run M2 model
                params[Mu_r] = 0.1 *2;//M/t //Elena: DLBCL: Uncomment to run M8 model (2-fold)
                params[sigma_r] = 2.6 *2;//M/t //no units//Elena: DLBCL: Uncomment to run M8 model (2- fold) ;
                break;
            }
            case NbModels_pars:
            {
                cerr<<"Error: Wrong parameter DLBCL identifyer"<<endl;
                break;
            }
        }

        params[sigma_p] = 9 ; //no units //Elena: DLBCL: Uncomment to run reference model
        params[sigma_b] = 100;//no units //Elena: DLBCL: Uncomment to run reference model
        params[Mu_p] = 1E-6 ; //M/t //Elena: DLBCL: Uncomment to run reference model
        params[l_p] = 1  ; //1/t //Elena: DLBCL: Uncomment to run reference model
        params[l_b] = 1  ; //1/t //Elena: DLBCL: Uncomment to run reference model
        params[l_r] = 1  ; //1/t //Elena: DLBCL: Uncomment to run reference model
        params[k_p] = 1 ;//M //Elena: DLBCL: Uncomment to run reference model
        params[k_b] = 1 ;//M //Elena: DLBCL: Uncomment to run reference model
        params[k_r] = 1 ;//M //Elena: DLBCL: Uncomment to run reference model
        params[init_p] = 0.1; //M //Elena: DLBCL: Uncomment to run reference model
        params[init_b] = 11.258310; //M //Elena: DLBCL: Uncomment to run reference model
        params[init_r] = 0.1; //M //Elena: DLBCL: Uncomment to run reference model
        params[bcr] = 0; //binary 0/1 //Elena: DLBCL: Uncomment to run reference model
        params[cd40] = 0; //binary 0/1 //Elena: DLBCL: Uncomment to run reference model

        setBaseParametersDone();
    }


    //These parameters depend on the state of each cell and are set in main loop during the cell interaction.
    void network::setDinamicParameters(double p, double b, double r, double bcr_,double cd40_){
        params[init_p] = p;
        params[init_b] = b;
        params[init_r] = r;
        params[bcr] = bcr_;
        params[cd40] = cd40_;


    }


    void network::initialise(){ // don't touch to parameters !
        val.clear();
        val.resize(NbVariables, 0.0);
        init.clear();
        init.resize(NbVariables, 0.0);

        init[p]   = params[init_p]; //M
        init[b]   = params[init_b]; //M
        init[r]   = params[init_r]; //M

        val[p] = init[p];
        val[b] = init[b];
        val[r] = init[r];

        t = 0;
        initialiseDone();
    }

//Elena: define network derivatives depending on DLBCL enum
    void network::derivatives(const vector<double> &x, vector<double> &dxdt, const double t){
        
        switch(dlbcl_dervs){
            case Healthy_dervs:{
                
                dxdt[b] = params[Mu_b] + params[sigma_b] * hillI(x[p], params[k_p]) * hillI(x[b], params[k_b]) * hillI(x[r], params[k_r])  - (params[l_b] + (params[bcr] * hillI(x[b], params[k_b]))) * x[b]  ; //Elena: DLBCL: Uncomment to run reference model.
                
                dxdt[p] = params[Mu_p] + params[sigma_p] * hillI(x[b], params[k_b]) + params[sigma_p] * hillA(x[r], params[k_r]) - params[l_p] * x[p]   ; //Elena: DLBCL: Uncomment to run reference model.
                
                break;
            }
            case M1:{
                 
                dxdt[b] = params[Mu_b] + params[sigma_b] * hillI(x[p], params[k_p]) * hillI(x[r], params[k_r])  - (params[l_b] + (params[bcr] * hillI(x[b], params[k_b]))) * x[b]  ; //Elena: DLBCL: Uncomment to run M1 model where BCL6 autoinhibitory loop has been removed.
                dxdt[p] = params[Mu_p] + params[sigma_p] * hillI(x[b], params[k_b]) + params[sigma_p] * hillA(x[r], params[k_r]) - params[l_p] * x[p]   ; //Elena: DLBCL: Uncomment to run reference model.
                
                break;
            }
            case M3A:{
                               
                 dxdt[b]= params[Mu_b] + params[sigma_b] * hillI(x[p], params[k_p]) * hillI(x[b], params[k_b])  - (params[l_b] + (params[bcr] * hillI(x[b], params[k_b]))) * x[b]  ; //Elena: DLBCL: Uncomment to run M3A model where IRF4 inhibitory effect over BCL6 has been removed.
                dxdt[p] = params[Mu_p] + params[sigma_p] * hillI(x[b], params[k_b]) + params[sigma_p] * hillA(x[r], params[k_r]) - params[l_p] * x[p]   ; //Elena: DLBCL: Uncomment to run reference model.
                break;
            }
            case M3B:{

                dxdt[b]= params[Mu_b] + params[sigma_b] * hillI(x[r], params[k_r]) * hillI(x[b], params[k_b])  - (params[l_b] + (params[bcr] * hillI(x[b], params[k_b]))) * x[b]  ; //Elena: DLBCL: Uncomment to run M3B model where BLIMP1 inhibitory effect over BCL6 has been removed.
                dxdt[p] = params[Mu_p] + params[sigma_p] * hillI(x[b], params[k_b]) + params[sigma_p] * hillA(x[r], params[k_r]) - params[l_p] * x[p]   ; //Elena: DLBCL: Uncomment to run reference model.
                break;
            }
            case M4:{
                
                dxdt[b]= params[Mu_b] + params[sigma_b] * hillI(x[r], params[k_r]) * hillI(x[b], params[k_b])  - (params[l_b] + (params[bcr] * hillI(x[b], params[k_b]))) * x[b]  ; //Elena: DLBCL: Uncomment to run M3B model where BLIMP1 inhibitory effect over BCL6 has been removed.
                dxdt[p] = 0  ; //Elena: DLBCL: Uncomment to run M4 model.
                break;
            }
            case NbModels_dervs:{
                cerr<<"Error:parameter DLBCL identifyer"<<endl;
                break;
            }
        }
        
        dxdt[r] = params[Mu_r] + params[sigma_r] * (hillA(x[r], params[k_r])) + (params[cd40] * hillI(x[b], params[k_b])) - params[l_r] * x[r]   ; //Elena: DLBCL: Uncomment to run reference model.
    }
