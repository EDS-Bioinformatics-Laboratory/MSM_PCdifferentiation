#ifndef i_setparam
#define i_setparam

///§§§ Philippe 2017-10-12
#include <fstream>

#include <string>
#include <vector>
#include "dynarray.h"
using namespace std;

/** @brief How to add new parameters :
 * Add them in the parameter file, two lines each
 *      syntax :
 *          kennwerte(keyword):
 *          value
 * New parameters should be entered in the class :
 *      Werte,
 *          as new fields.
 * The following functions have to be updated with this new parameters :
 *      Werte()
 *      Werte.fPut(...);
 *      Werte.fGet(...);
 *      Werte.show();
 * External use is done using :
 *      parameter par;
 *      par.wahl(...)
 *      par.Value.xxx (field for this parameter)
 */

/** @brief Type definition for suffixes with increasing ID each time (0000, 0001 etc ...)*/
typedef char suffix[5];
/** @brief to increase the suffix by 1. */
void addchar(suffix &tmp);

/** @brief Maximum dimension of tables (fields of Werte), ex : max number of initial sequences */
const int MAXDIM = 2000;
const int MAXDIMSMALL = 30;
const int MAXKENNTEXTSIZE = 2000;           // Philippe 01.02.2016

enum representation {
   GC,tumour,islet
};
enum Ig_classes {
   IgM,IgG,IgE,IgA,nIg_classes
};

/** @brief: Public class to carry parameter values (as fields): */
class Werte {
  public:
   /** @brief constructor, calls ini2d by default : */
   Werte();
   bool show_missing_pars;

   Werte(const Werte &) { cerr << "Err : Werte copy constructor called" << endl; }
   /** @brief fills the fields with default values for running a 2D simulation. Note : static
    * default
    * values are overrided */
   void ini2d();
   /** @brief fills the fields with default values for running a 3D simulation. Note : static
    * default
    * values are overrided */
   void ini3d();

   /** @brief : function that associates the ID (integer) of a parameter to its name/line (kenntext)
    * to be read in the parameter file */
   const char* kenntext(int nummer);

   /** @brief : Output all the parameter values in a new parameter file, formatted with
    * 'kenntext:\nvalue\n' */
   ofstream &fPut(ofstream &s);

   /** @brief : parse a parameter file and fills all the fields of Werte with it. When a parameter
    * is
    * not found (its kennwerte), an
    *  error is raised. */
   void fGet(char * parname, bool transform2rate);

   /** @brief (tool function :) Looks for a specific parameter in a parameter file :
    *  reopens a parameter file using a specified ifstream, and just stops at the level of the
    * description of parameter n (kenntext(n)),
    *  meaning the next time s >> is used, the value for this particular parameter should follow.
    *  @param parname   (file name)   parameter file to parse
    *  @param s         (the ifstream that will be used afterwards to get the value of the
    * parameter)
    *  @return 1 if found, 0 if not found (short)   */
   short fFind(char * parname, ifstream &s, int n);



   static void compareParameterSets( Werte &w1, Werte &w2);

   /** ------ List of parameters fields ----- */

   // General properties:
   // ===================
   // Betriebssystem: 0=Unix; 1=Windows
   short system,outputfiles,timevalues,show_Ki67,safety_checks;
   representation show_mode;
   long ini_random,late_ini_random;
   // Output-Restriction
   vector<int> file_output;
   // Border conditions: cyclic=1, fixed=0
   // short int CyclicBorder;
   /* Calculate as it is (=0)
    * Use int-Population for calculation (=1)
    * Use Rounded values for calculation (=2)
    * Use old values for calc. and refresh them after a passed cycle (=3)
    */
   // short int UseIntPart;

   ////Philippe : to include later in main
   string prefixSigFiles;


   // Dimension of arrays:
   long CB_Narray,CC_Narray,TC_Narray,FDC_Narray,OUT_Narray,STROMA_Narray,BETA_Narray;

   // Shape space:
   // ============
   // Dimension des Shapespace
   int DimShapeSpace;
   // Metrik
   short int metrik;
   // Number of B-Cell States, Number per Dimension
   long int SSStates, SSRangePerDim;
   // Total Number of presented Antigen Epitops (int-type):
   // Number of Antigen Peaks in its Shapespace (int-type):
   int totalA,APeakNumber;
   // Selection term weighted with relative number of presented epitops =1
   // not weighted =0
   // Random oder vorgegeben
   dynarray<long int> takeA;
   // Antigen fractions
   vector<double> ag_fraction;
   // Width and amplitude of gaussian affinity-weight function
   double GammaGauss,amplitudeGauss;

   // Sequence space:        // Philippe 01-02-2016
   // ============
   // Type of affinity function : 0= standard (saham's)  1= normalized to max_affinity_cluster 2=
   // sliding windows of size max_affinity_cluster
   int type_affinity_function;

   int use_logarithmic_seq_affinity;
   // Use sequence space (1/0): (if not, use shape space)
   int use_sequence_space;
   // Length of sequences:
   int size_sequences;
   // Probability of mutation, per base.
   double sequence_mut_per_base;
   // Number of Initial Antigen sequences (int-type):
   int init_antigen_sequences;
   // Fix Antigen Sequence presentation (max 1000 values):
   vector<string> initAntigenSeqs;
   // Maximum Hamming distance between antigens:
   int max_hamming_antigens;
   // Minimum Hamming distance between antigens:
   int min_hamming_antigens;

   // note : the initial number of BCRs is determined by totalBSS : total B seeder cells, already
   // defined in the centroblasts part.
   // Fix initial Repertoire distribution (max 1000 values):
   vector<string> initBCRSeqs;
   // Maximum Initial Hamming distance between BCRs:
   int max_hamming_BCRs;
   // Minimum affinity of initial BCRs to Antigens:
   double min_initial_affinity_BCRs;
   // Maximum affinity of initial BCRs to Antigens:
   double max_initial_affinity_BCRs;

   // Fix initial Repertoire distribution for T cells:
   vector<string> initTCRSeqs;
   // Maximum Initial Hamming distance between TCRs:
   int max_hamming_TCRs;
   // Minimum affinity of initial TCRs to Antigens:
   double min_initial_affinity_TCRs;
   // Maximum affinity of initial TCRs to Antigens:
   double max_initial_affinity_TCRs;




   ////§§§ Philippe 2017-10-12
   bool use_predefined_tc_dynamics;
   vector<double> TCdynT;
   vector<double> TCdynNb;




   // Specifity of sequences affinity (double R):
   double R_affinity;
   // Optimum affinity cluster size (affinity doesn't increase beyond it):
   int max_affinity_cluster;

   // Arup space:        // Philippe 20-03-2016
   // ============
   int use_arup_space;
   int arup_length_sequences;
   int arup_N_conserved;
   int arup_N_mutates;
   int arup_N_shielded;
   int arup_nb_ini_antigens;

   vector<string> arup_ini_antigens;
   vector<double> arup_ag_fraction;

   int arup_nb_mutations_gen_strains;
   double arup_threshold_activation;
   double arup_h_min;
   double arup_h_max;

   vector<string> arup_ini_bcrs;
   double arup_mutation;
   double arup_proba_lethal_mut;
   double arup_proba_affecting_mut;
   double arup_proba_silent_mut;

   vector<double> arup_law_mut_Xs;
   vector<double> arup_law_mut_Densities;

   double arup_alpha;
   double arup_hprime_min;
   double arup_hprime_max;
   double arup_hmut_min;
   double arup_hmut_max;

   // Signals:
   // ========
   // mode of signal treatment
   short signal_mode;
   short objects_transparent;
   // Dirichlet boundary for signal differ2CC
   double bound_differ2CC,bound_CXCL12,bound_CXCL13,bound_ab,bound_ag,bound_SEMA4D;
   // critical concentration for desensitisation
   double CXCL12crit,CXCL13crit;
   double CXCL12recrit,CXCL13recrit;
   // Diffusion constant
   double D_differ2CC,D_CXCL12,D_CXCL13,D_antibody,D_antigen,D_SEMA4D;
   // take some signals from file
   dynarray<bool> fix_signals;
   // define properties of a glucose field
   bool fix_glucose_gradient;
   double const_dynamic_glucose_field,fix_glucose_gradient_min,fix_glucose_gradient_max;

   // Space properties:
   // =================
   // Dimension of lattice
   short int DimSpace;
   // Radius of GC in microm
   double GC_radius;
   // Grid size in every direction
   double gridsize[3];
   // lattice constant in microm
   double dx,dx_signal;
   // Shape of reaction volume
   short vol_shape;
   // obstacles:
   short obstacles;
   double wall_level,collagen_density,collagen_cluster;
   int wall_width,slit_number,slit_width;

   // Time and phases:
   // ================
   // Starting time for mutation
   double Start_Mutation,Start_Differentiation;
   // Starting time for Output
   double StartOutput;
   // Stop new BC influx
   double newBCinflux_stop;
   // Zeitschritte, Anfang und Ende
   double deltat, tmin, tmax;
   // Zahl der Zeitschritte zwischen 2 Outputs
   // ### Hier lieber die Zeit zwischen 2 Outputs einlesen und umrechnen
   long int ToFileStep;

   // Cells in general
   // ================
   // adhesion
   double adhesion_time;    // time needed to establish adhesion of a cell fragment with others
   // chemotaxis
   double chemo_max,        // maximum weight relative to diffusion,
          chemo_steep,  // steepness of reduction for smaller chemokine gradients
          chemo_half;   // chemokine-gradient of half weight (of chemo_max)

   /// Philippe 2017-10-22
   double chemo_max_tc,        // maximum weight relative to diffusion,
          chemo_steep_tc,  // steepness of reduction for smaller chemokine gradients
          chemo_half_tc;   // chemokine-gradient of half weight (of chemo_max)
   bool use_specific_tc_chemotaxis;
    /// Philippe 2017-10-26
   double proba_TC_CC_interaction;


   // Nutrient consumption in general
   double use_glucose,use_oxygen,use_glucose_pro,use_oxygen_pro,critical_nutrient;
   double bound_glucose,bound_oxygen;
   double D_glucose,D_glucose_H2O,D_oxygen,D_oxygen_H2O;
   // Macrophages
   double p_macrophage;

   // Motility
   bool allow_exchange;    // Makes exchange of contact inhibited cells possible >=v7.05.1

   // random polarity
   short use_specific_turning_angles;

   // Centroblasts or blast1:
   // =======================
   // Total Number of initial B-Cells:
   long int totalB,totalBss;
   // Rate of new BC influx to GC:
   double newBCinflux_rate;
   // width of smooth BC influx stop:
   double smooth_stopBCinflux;
   // radius of centroblasts
   double CB_radius;
   // position in shape space
   dynarray<long int> takeB;
   // restrictions for the position of seeder cells in shape space:
   double min_seeder_dist,max_seeder_dist;
   // position in space
   dynarray<long int> posCB;

   // Proliferation, growth, and mutation
   double
      proliferate, // Rate per hr
      dx_CB,       // maximal distance for CB-proliferation from dividing cell
      grow,
      tolight,
      CB_maxvolume4differ2CC, // fraction of total volume up to which CB may differentiate
      mutation,mutation_after_tc,mutation_after_dec_tc,mutation_affinity_exponent,
      CB_fixed_times_of_divisions,CB_fixed_times_of_divisions_in_expansion,CB2OUT_prob;
   int reset_antigen_after_collection;
   short present_specific_ag2TC;
   short fixed_time_of_divisions_mode;
   bool smooth_differentiation,smooth_dif2out,exit2tz;
   double smooth_differentiation_time,smooth_dif2out_time;
   double CB_dt_G0,CB_dt_G1,CB_dt_G2,CB_dt_S,CB_dt_M;
   double CB_dtphase_width;
   double t_inject_BrdU, deltat_inject_BrdU, BrdU_detection_threshold;
   int n_inject_BrdU;
   bool transmit_CC_delay_to_CB_cycle;
   bool retain_ag,ag_loaded_CB_diff2output,ag_loaded_CC_directly2TFH,
        ag_loaded_CB_stop_mutation,ag_deleted_in_fresh_CC;
   double divide_ag_asymmetric,asymmetric_polarity_index,smooth_PI,BC_ag_preloaded;

   // receptors
   short CBreceptor_use;
   double
      CBreceptor_dissociation,
      CBreceptor_binding,
      CBreceptor_total,
      CBreceptor_activation;

   // adhesion
   double CB_max_adhesion;    // maximum adhesion force in % of full stickness

   // Motility:
   double D_CB;             // Diffusion (alternative to cell velocity -- if v is used take -1)
   double v_CB;             // CB-velocity (alternative to diffusion -- if D is used take -1)
   double v_CB_width;       // defines a width of Gauss distributed v_CB values (-1 for fixed)
   double CB_smoothmove;    // Distributes a barycenter movement thought to overcome one
                            // lattice constant on several time steps. In each time step
                            // only a subset of fragments is moved for values >1.
   double CB_persistence;   // average time gap in minutes between changes of direction of the
                            // cell polarity.
   short CB_v_modi;         // modus of velocity state treatment
   short CB_n_v_states;     // # of velocity states
   double v_CB_switch_deltat;    // Mean duration in a v-state in minutes
   double v_CB_factor;      // for 2 velocities: the factor by which the velocity is reduced

   // Shape parameter
   double distance_tolerance,half_tolerance_deformation;
   /* For cell-fragment-movement:
    * tolerance for the distance to the barycenter
    * to chose the target lattice point. The real
    * tolerance parameter is calculated according
    * to apparent deformation. Therefore, the deformation
    * at which the tolerance is half maximal is given
    * separately. */
   double CB_D_cytosol;     // Diffusion constant for fragments in the cytosol.
   // =D_CB means equal diffusion for fragments and cell
   // in the non-deforming limit --> surface tension of the cell.
   double v_CB_cytosol;     // Alternative to CB_D_cytosol (one of both has to be set to -1).
   double CB_elongation;    // Cell elongation by active movement
                            // 1 means the shift of the barycenter per time step is
                            // one lattice constant; larger values lead to cell elongation.
   double CB_K_elongation;    // Elongation in units of spherical cell radius, at which
                              // the reshaping force is half maximal (in Hill equation).

   // blast2:
   // =======================
   // Total Number of initial B-Cells:
   long total_blast2;
   // radius of centroblasts
   double blast2_radius;
   // position in shape space
   // dynarray<long int> takeB;
   // position in space
   dynarray<long int> pos_blast2;
   // maximal distance for CB-proliferation from dividing cell
   double dx_blast2;
   // Raten pro Zeit (per h)
   double blast2_proliferate,blast2_grow;
   // For cell-movement: tolerance for the distance to the barycenter
   // to chose the target lattice point
   double blast2_distance_tolerance,blast2_half_tolerance_deformation;
   // Diffusion
   double D_blast2;

   // Centrocytes:
   // ============
   // Dauern
    //Elena: Network parameters
   double cd40;
   double bcr;
   double BLIMP1th;
   double pBLIMP1;
   double pIRF4;
   double pBCL6;
   double type_CD40_signal;
    double dlbcl_parameters;
    double dlbcl_derivatives;
    
   double CC_test_delay,CC_ICAM_delay;
   // Test with FDC necessary?
   unsigned short CC_FDC_selection;
   bool collectFDCsignals;
   double collectFDCperiod;
   double prob2kill_noFDCcontactBCs;
   // Wahrscheinlichkeiten
   double TCell,output,output_DEC,FDCsignalling;
   double shrink,apoptosis,apoptosis4FDCselected,macrophage,ignore_affinity,
          selection,ccdiff,ccdiff_delay,ccdiff_delay_DEC,final_differentiation_rate,
          tc_search_duration_per_FDCcontact, // add this to the search time for each collected
                                             // Ag-portion
          tc_search_duration_fixed, // BC search duration for Tfh for mode==1
          TC_time,          // Duration of TC-CC interaction in hours
          TC_time_width,    // width of this duration in hours
          TC_rescue_time;   // Minimum duration for selection
   short mode_of_setting_TC_time, tc_search_duration_mode;
   bool negativeTCselection, BCstaysonTCbyTCtime, multipleTFHcontacts;
   // Motility:
   double D_CC;             // Diffusion (alternative to cell velocity -- if v is used take -1)
   double CXCR5down;        // rate of CXCR5 downregulation (-1 for none)
   double CXCR4down;        // rate of CXCR4 downregulation (-1 for none)
   double v_CC;             // CB-velocity (alternative to diffusion -- if D is used take -1)
   double v_CC_width;       // defines a width of Gauss distributed v_CC values (-1 for fixed)
   double CC_persistence;   // average time gap in minutes between changes of direction of the
                            // cell polarity.
   short CC_v_modi;         // modus of velocity state treatment
   short CC_n_v_states;     // # of velocity states
   double v_CC_switch_deltat;    // Mean duration in a v-state in minutes
   double v_CC_factor;      // for 2 velocities: the factor by which the velocity is reduced
   short use_ab_dynamics;    // =0 for old affinity model;
                             // =1 feedback of output average affinity on binding threshold
                             // =2 same using max affinity of output
                             // =3 using average affinity of produced and injected antibodies as
                             // threshold
   double initial_ab_affinity;    // -1 for take seeder cell average
   bool ignore_apoptotic_CC;    // 0: include them; 1: ignore them for cell number analysis
   short CC_apoptotic_motility_mode;
   double p_apo_randomwalk;

   // class switch
   short do_switch_classes;
   static const int switch_dimension = nIg_classes * nIg_classes;
   double switch_matrix[switch_dimension];
   double IgE_BCRlevel,IgE_factor_cellcycle,IgE_factor_divisions,CC_IgE_prob_CXCR5down;

   ///§§§ Philippe 23/03/2017
   double IgG_BCRlevel;             // 373
   double IgG_factor_cellcycle;     // 374
   double IgG_factor_divisions;     // 375
   double Founder_IgX[nIg_classes]; // 376
   double decay_proba_switch;       // 377 for all probabilities of switching
   double IgG_factor_leaving;       // 378
   double Affinity_threshold_IgG;   // 379
   double Int_Antigen_threshold_IgG;// 380
   double tc_help_IgG;              // 381
   double stddev_initial_divisions; // 386  if != 0, then apply it
   double stddev_DND;               // 387  if != 0, then apply it

   double time_tc_selection_block;      // 382
   double factor_tc_selection_block;    // 383  if == 1, then no blocking
   double time_DND_block;               // 384
   double factor_DND_block;             // 385  if == 1, then no blocking
   double factor_founder_div_block;     // 388
   ///§§§
   int    mode_tc_selection_block;      // 389

   //bool authorize_tc_selection;

   //double switching_clock_mean, switching_clock_stddev;

   bool pMHC_dependent_division, signal_dependent_number_of_divisions;
   double pMHC_dependent_P_max, pMHC_dependent_K, pMHC_dependent_nHill, pMHC_dependent_P_min,
          pMHC_dependent_P_standard, pMHC_dependent_pMHC_of_2divisions;
   double TFHsignal_dependent_K, TFHsignal_of_P0divisions;

   bool ICOSL_dependent_Tfh_signals, ICOSL_memory;
   short ICOSL_upregulation_mode;
   double ICOSL_upregulation_time;

   // Philippe 2018-01-05 New version from git
   double dT_FoxO, dT_mTORC1, nFoxO, KFoxO;

   // TC:
   // ===============
   // Total Number of initial T-Cells:
   long int totalTC;
   // radius of TC
   double TC_radius;
   // Motility:
   double v_TC;             // TC-velocity
   double v_TC_width;       // defines a width of Gauss distributed v_TC values (-1 for fixed)
   double v_TC_CC;          // TC-velocity if encountering a CC
   double TC_persistence;   // average time gap in minutes between changes of direction of the
                            // cell polarity.
   double north_weight;     // tendency "a" to walk north: p=(1-a)r+an with r: random, n: north
   // Selection of CCs
   // From v16.07.15 and before and controls whether CCs need to interact with TFH:
   short TC_CC_selection;    // 1 if TC rescue CC from apoptosis
   /* Tfh-specific helper function, introduced in August 2016.
    *  Requires a Tfh-specific property Tfh_quality.
    *  Tfh_quality is used to scale the helper function up and down and is different from the
    *  <TCell> which reduces the selection probability globally. Tfh_quality is used to scale
    *  the number of divisions induced in selected BCs. */
   /* Mode of Tfh-CC interaction :
    *  0: Constant selection probability as in v16.07.15 and before (no Tfh-specific variation)
    *     --> Tfh_quality=1 for all Tfh.
    *  1-2: Each Tfh is a helper with a specific Tfh_quality \in [0,1]
    *  1: --> Tfh_quality is sampled from a Gaussian with maximum 1 and width given below
    *  2: --> Tfh_quality is derived from a position in AffinitySpace given below
    */
   short TFH_CC_selection_mode;
   double TFH_CC_selection_gauss_width;
   long int TFH_ASpos;
   // Division:
   bool do_TC_division;
   double TC_doubling, TC_meancycle, TC_cyclewidth, dx_TC;
   int TC_Ndivisions;

   // OUT:
   // ===============
   double mk_ab,pm_differentiation_time;
   double v_OUT;             // OUT-velocity
   double v_OUT_width;       // defines a width of Gauss distributed v_OUT values (-1 for fixed)
   double OUT_persistence;   // average time gap in minutes between changes of the polarity

   // Antibodies
   // =============
   // Antigen-Antibody binding
   double ic_k_on,ic_k_off,ag_threshold;
   int antibodies_resolution;
   double antibodies_production,antibodies_degradation;
   double k_ic_exp_min,k_ic_exp_max;
   double N_GC,V_blood;
   double inject_antibody,injected_antibody_affinity,inject_antibody_time;
   long inject_antibody_ASindex;

   // Photoactivation
   // ================
   bool photoactivation;
   double photoactivation_t0,
          photoactivation_x0,photoactivation_y0,photoactivation_z0,
          photoactivation_delta_x,photoactivation_delta_y,photoactivation_delta_z;
   bool def_DEC205,inject_antiDEC205OVA,DEC205_induce_CBdifferentiation,retain_DEC205_ag;
   double def_DEC205_t0,inject_antiDEC205OVA_t0,antiDEC205OVA_tend,p_DEC205,
          TC_dec205ova_time,TC_factor_dec205ova,DEC205_p_factor,DEC205_forces_output;

   // FDC:
   // ===============
   // total number
   int FDCnumber;
   // Percentage of FDC network in the total GC volume
   double FDCnetwork;
   // position of FDCs
   dynarray<long int> posFDC;
   // FDC: length of arms
   int FDClength;
   // FDC: dendrites treated transparent
   short int FDCtransparent;
   // signal production of FDC by vesicle exocytosis or continuous
   short FDCvesicle;
   // Signal production rate
   double mksignal,mkCXCL12,mkCXCL13,mk_SEMA4D;
   // Antigen comsumption
   double ag_per_FDC,ag_saturation_FDC;
   // How to distribute Ag on FDC fragments
   short ag_distribution_mode,ag_detection_mode;

   // BETA-cells:
   // =======================
   long BETA_Nini;            // Total Number of initial betacells:
   dynarray<long int> BETA_pos;         // position in space

   // Proliferation, growth, motility
   short
      BETA_v_modi,        // modus of velocity state treatment
      BETA_n_v_states;    // # of velocity states
   double
      BETA_radius,
      BETA_proliferate,   // Rate per hr
      BETA_max_pro,       // maximal distance for CB-proliferation from dividing cell
      BETA_grow,
      BETA_shrink,
      BETA_max_adhesion,  // maximum adhesion force in % of full stickness
   // motility
      BETA_persistence,   // average time gap in minutes between changes of direction of the
   // cell polarity.
      BETA_v,             // velocity
      BETA_v_factor,      // for 2 velocities: the factor by which the velocity is reduced
      BETA_v_switch_deltat, // Mean duration in a v-state in minutes
      BETA_v_cytosol,     // Strength of reshaping forces (cytosolic elements speed)
      BETA_elongation,    // Cell elongation by active movement
                          // 1 means the shift of the barycenter per time step is
   // one lattice constant; larger values lead to cell elongation.
      BETA_K_elongation,  // Elongation in units of spherical cell radius, at which
   // the reshaping force is half maximal (in Hill equation).
      BETA_distance_tolerance,
      BETA_half_tolerance_deformation,
   /* For cell-fragment-movement:
    * tolerance for the distance to the barycenter
    * to chose the target lattice point. The real
    * tolerance parameter is calculated according
    * to apparent deformation. Therefore, the deformation
    * at which the tolerance is half maximal is given
    * separately. */
      BETA_smoothmove;   // Distributes a barycenter movement thought to overcome one
   // lattice constant on several time steps. In each time step
   // only a subset of fragments is moved for values >1.

   // Cell tracking:
   // ===============
   // Zahl der Punkte auf der v-Achse fuer die v-Verteilung der diffusiven Bewegung
   int v_resolution,s_resolution,alpha_resolution;
   // intervall of velocities in microns/min (used for time intervall averaged velocities)
   double delta_v,delta_s,delta_alpha,trackfrom,trackuntil,track_delta_t;
   // Number of cells of different celltypes to be tracked
   int tALL,tCB,tCC,tOUT,tTC,tBETA;

  private:
   /** @brief a common initializer for fields values in 2D and 3D (called bt ini2d and ini3d) */
   void inialld();
};

class betaWerte {
  public:
   betaWerte();
   void ini();
   const char* kenntext(int nummer);
   ofstream &fPut(ofstream &s);
   short fFind(char * parname, ifstream &s, int n);
   void fGet(char * parname);
   void show();

   // The parameters:
   // ===================

   short use_Nernst,set_leakage_zero,use_inactivation,
         use_dynamic_tau_K_V,use_dynamic_tau_Na_V,use_dynamic_tau_fNa_V,
         use_dynamic_H_K_Ca,use_voltage_gating_K_Ca,
         use_dynamic_IP3,use_dynamic_tau_IP3;

   double dt;    // step size dt = 1 min (internally this is half of the value given here)
   double dy;      // ### maximum tolerance of double t-step deviation (in cell fractions)
   double t_0;                // inital t
   double t_max;              // final t
   double dt_output;     // step size of writing in output file (every hour)

   enum beta_proteins {
      NaK,K_ATP,K_V,K_Ca,sK_Ca,Na_V,fNa_V,NCX,PMCA,
      Ca_L,Ca_T,SERCA,IP3,gap,N_beta_proteins
   };
   double * rho;

   double V_0,
          R_bc,Sur_ER,Vol_ER,C_m,
          K_ext,Na_ext,Ca_ext,Ca_ER_0,
          Vbar_Ca_delta,
          cal_0,K_cal,buf_0,K_buf,buf_ER_0,K_buf_ER,
          glu_0,IP3_0,K_0,Na_0,Ca_0,
          T,
          rho_NaK,
          Ihat_NaK,   H_NaK,     n_NaK,     H2_NaK,    n2_NaK,     alpha_NaK,
   // rho_K_ATP,
          gbar_K_ATP, tau_K_ATP, s_h_K_ATP,  kappa_K_ATP,
   // rho_K_V,
          gbar_K_V,   tau_K_V,   V_h_K_V,    kappa_K_V, theta_K_V, W_h_K_V,    lambda_K_V,
   // rho_K_Ca,
          gbar_K_Ca,  H_K_Ca,    n_K_Ca,    V_h_K_Ca,  kappa_K_Ca, tau_K_Ca,
   // rho_sK_Ca,
          gbar_sK_Ca, C_sK_Ca,   kappa_sK_Ca,tau_sK_Ca,
   // rho_Na_V,
          gbar_Na_V,  tau_Na_V,  V_h_Na_V,   kappa_Na_V,theta_Na_V,W_h_Na_V,   lambda_Na_V,
   // rho_fNa_V,
          gbar_fNa_V, tau_fNa_V, V_h_fNa_V,  kappa_fNa_V,
   // rho_NCX,
          Ihat_NCX,   H_NCX,     n_NCX,      alpha_NCX,
   // rho_PMCA,
          Ihat_PMCA,  H_PMCA,    n_PMCA,     alpha_PMCA,
   // rho_Ca_L,
          gbar_Ca_L,  tau_Ca_L,  V_h_Ca_L,   kappa_Ca_L,theta_Ca_L,W_h_Ca_L,
          lambda_Ca_L,C_Ca_L,    n_Ca_L,
   // rho_Ca_T,
          gbar_Ca_T,  tau_Ca_T,  V_h_Ca_T,   kappa_Ca_T,  theta_Ca_T,W_h_Ca_T,   lambda_Ca_T,
   // rho_SERCA,
          Ihat_SERCA, H_SERCA, n_SERCA,
   // rho_IP3,
          g_IP3_max,  gbar_IP3,  C_IP3_act,    n_IP3_act, tau_IP3,
          Cbar_IP3_inh, n_IP3_inh, theta_IP3,
          P_IP3,     kappa_IP3,
          k_IP3_plus,k_IP3_minus,C_P,n_P,
   // rho_gap,
          gbar_gap,   tau_gap;
   bool gap_dynamic;

   enum random_laws {
      equal,poisson,gauss,N_random_laws
   };
   bool randomise_beta_proteins;
   random_laws randomisation_type;
   double randomisation_range;
};

// Klasse zur Verwaltung von Werte:
class hyphasmaParameter {
  private:
   static const short namelength = 500;

  public:
   ////§§§ Philippe 21-03-2017
   static const double N_A; // mol^-1
   static void goon();

   hyphasmaParameter() { }
   hyphasmaParameter(const hyphasmaParameter &) { cerr << "Parameter Copy Constructor called" << endl; }
   char name[namelength];
   char logfile[namelength];

   // Parameterwerte:
   Werte Value;
   // Parameter values for betacells:
   betaWerte betaValue;

   // reset the random-number generator in a parameter file to a value determined by time
//   void reset_random_ini(const char*, int);

   // Externer Aufruf:
   int wahl(const char*, bool transform2rate, bool show_missing);

  private:
//   char first;

   // Dateiname mit Suffix:
   char namesuffix[namelength];

   // temporaere Variablen zur Pufferung:
   char ntmp[namelength];
   Werte vtmp;

   // Prozeduren:
   void in_datei(suffix suff, suffix log);
   void save();
   void recover();
   void read(bool transform2rate);

 ///§§§ Philippe : put it public to save a parameter as different files
  public:
   void write(string filename = string("")); // if not given, will be saved with name nameSuffix
};

#endif
