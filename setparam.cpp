/* Routinen fuer die Erstellung und Veraenderung der
 *      Paramterdatei
 */
#include "setparam.h"
#include <string.h>
#include <math.h>
#include <fstream>
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Parameter-Klassen-Konstruktoren
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// The following values are set in a way that parameter
// files from previous versions remain compatible and
// can still be used. Default parameter values have to
// be set such that the program can still run as before
// the introduction of the new parameters if they are not
// specified.
// ########################################################

   ////§§§ Philippe 21-03-2017
const double hyphasmaParameter::N_A = 6.02205e+23; // mol^-1


// increases the numeric text tmp by +1
void addchar(suffix &tmp) {
   // suffix declared in setparam.h
   int i = 3;
   char weiter = 1;
   while (i >= 0 && weiter == 1) {
      if (tmp[i] != '9') {
         ++tmp[i];
         weiter = 0;
      } else {
         tmp[i] = '0';
         --i;
      }
   }
   if (weiter == 1) {
      cout << "Too large number in addchar(..) !\n";
   }
}
// initializes the default parameter values that are common between 2D and 3D simulations
void Werte::inialld() {
   // Set System to Unix
   system = 0;
   // Safety-checks (0=less, 1=more, 2=a lot)
   safety_checks = 1;
   // Keine Vorarbeit des Zufallsgenerators
   ini_random = 0;
   late_ini_random = 0;
   // Schreibe alles raus
   outputfiles = 0;
   // Modus of space representation
   show_mode = GC;
   // Hebe Ki67 hervor?
   show_Ki67 = 0;
   // Lese Zeiten nicht Raten
   timevalues = 1;
   // Initial dimensions of arrays
   CB_Narray = 20000;
   CC_Narray = 20000;
   TC_Narray = 2000;
   FDC_Narray = 300;
   OUT_Narray = 100;
   STROMA_Narray = 400;
   BETA_Narray = 0;

   prefixSigFiles = "";

   // Signals:
   D_differ2CC = 200;   // microm^2/h
   D_CXCL12 = 1000;     // microm^2/h
   D_CXCL13 = 1000;     // microm^2/h
   D_antibody = 2000;   // microm^2/h
   D_antigen = 2000;    // micron^2/h
   D_SEMA4D = 1000;
   signal_mode = 2;   // QUANTA
   bound_differ2CC = 0.0;
   bound_CXCL12 = 0.0;
   bound_CXCL13 = 0.0;
   bound_ab = 0.0;
   bound_ag = 0.0;
   bound_SEMA4D = 0.0;
   CXCL12crit = -1;
   CXCL13crit = -1;
   CXCL12recrit = -1;
   CXCL13recrit = -1;
   objects_transparent = 1;

   // Dimension of Shapespace
   DimShapeSpace = 4;
   // Metrik ist nicht euklidisch
   metrik = 1;
   // Number of B-Cell states in SS
   SSStates = 10000;
   // Resulting Range per dimension
   SSRangePerDim = long (pow(SSStates, (1. / DimShapeSpace)));
   // A shapespace initialization
   totalA = 1000;
   APeakNumber = 1;
   ag_fraction.reserve(100);
   ag_fraction.clear();
   for (int i = 0; i < 100; i++) {
      ag_fraction[i] = -1;
   }
   // Width of gaussian affinity weight function:
   GammaGauss = 2.8;
   amplitudeGauss = 1.0;

   use_arup_space = 0;
   arup_length_sequences = 46;
   arup_N_conserved = 18;
   arup_N_mutates = 22;
   arup_N_shielded = 6;
   arup_nb_ini_antigens = 1;

   // arup_ini_antigens.push_back(string("1111111111111111111111111111111111111111111111"));
   // arup_ag_fraction.push_back(1.0);

   arup_nb_mutations_gen_strains = 11;
   arup_threshold_activation = 10.8;
   arup_h_min = -0.18;
   arup_h_max = 0.9;

   // arup_ini_bcrs.push_back();
   arup_mutation = 0.003;
   arup_proba_lethal_mut = 0.3;
   arup_proba_affecting_mut = 0.2;
   arup_proba_silent_mut = 0.5;

   /*arup_law_mut_Xs.push_back(-6.4);
    *  arup_law_mut_Xs.push_back(-6.2);
    *  arup_law_mut_Densities.push_back(0.019047619);
    *  arup_law_mut_Densities.push_back(0.003809524);*/

   arup_alpha = 2;
   arup_hprime_min = -1.5;
   arup_hprime_max = 1.5;
   arup_hmut_min = -1.5;
   arup_hmut_max = 1e6;

   // No Restriction:
   for (int i = 0; i < MAXDIM; i++) {
      file_output[i] = -1;
      takeA[i] = -1;
      takeB[i] = -1;
      posCB[i] = -1;
      pos_blast2[i] = -1;
      BETA_pos[i] = -1;
      initAntigenSeqs[i] = string("-1");
      initBCRSeqs[i] = string("-1");
      initTCRSeqs[i] = string("-1");
      // the vectors for arup space are started empty, and are filled depending on the output -> not
      // allocated before, but can access their size by .size()
   }
   for (int i = 0; i < MAXDIMSMALL; i++) {
      fix_signals[i] = false;
   }



   //// Philippe 2017-10-12
   use_predefined_tc_dynamics = false;
   TCdynT.clear();
   TCdynNb.clear();
   proba_TC_CC_interaction = 1.0;





   // Sequence space:
   type_affinity_function = 0;   // 0: saham's model, 1: sahams renormalized by a cluster size. 2:
                                 // sliding window
   use_sequence_space = 0;
   size_sequences = 50;
   sequence_mut_per_base = 0.003; // from chakraborty's paper
   init_antigen_sequences = 10;
   max_hamming_antigens = 25;
   min_hamming_antigens = 1;
   max_hamming_BCRs = 25;
   min_initial_affinity_BCRs = 0.1;
   max_initial_affinity_BCRs = 0.9;
   max_hamming_TCRs = 25;
   min_initial_affinity_TCRs = 0.1;
   max_initial_affinity_TCRs = 0.8;
   R_affinity = 2.0;
   max_affinity_cluster
      = 10;   // depending of type_affinity_function : 1/ useless 2/ norm size 3/ size of the
              // sliding window

   // For lattice GCx
   vol_shape = 0;
   obstacles = 0;
   wall_level = 0.3;
   wall_width = 2;
   slit_number = 5;
   slit_width = 1;
   collagen_density = 0.3;
   collagen_cluster = 0.00333;

   // Timescales in h
   tmin = 0.;
   tmax = 504.;
   // Start Mutation
   Start_Mutation = 72.;
   Start_Differentiation = 72.;
   newBCinflux_stop = 72.;
   // Start Output
   StartOutput = 120.;

   // Adhesion
   adhesion_time = 0.166667;   // 10 seconds
   CB_max_adhesion = 0.0;
   // Chemotaxis
   chemo_max = 10.;        // Faktor
   chemo_steep = 1.e+10;   // l/mol
   chemo_half = 2.e-10;    // mol/l

   /// Philippe 2017-10-22
   use_specific_tc_chemotaxis = false;
   chemo_max_tc = 10.;        // Faktor
   chemo_steep_tc = 1.e+10;   // l/mol
   chemo_half_tc = 2.e-10;    // mol/l


   // Motility
   allow_exchange = false;
   use_specific_turning_angles = 0;

   // Photoactivation
   photoactivation = 0;
   photoactivation_t0 = 96.;
   photoactivation_x0 = 110.;
   photoactivation_y0 = 150.;
   photoactivation_z0 = 150.;
   photoactivation_delta_x = 100.;
   photoactivation_delta_y = 20.;
   photoactivation_delta_z = 20.;

   def_DEC205 = false;
   def_DEC205_t0 = 96.;
   p_DEC205 = 0.15;
   inject_antiDEC205OVA = false;
   inject_antiDEC205OVA_t0 = 96.;
   antiDEC205OVA_tend = -1.;    // default is that antiDEC205OVA acts for one round of selection
   TC_dec205ova_time = 0.0;     // unchanged TC-BC contact times
   TC_factor_dec205ova = 1.0;   // factor of TC number increase upon dec205ova stimulation
   DEC205_p_factor = 1.0;       // factor of prolongation of the CB state
   DEC205_induce_CBdifferentiation = false;
   DEC205_forces_output = 0.0;
   retain_DEC205_ag = false;

   // CB
   // totall number of initial B-Cells
   totalB = 3;
   totalBss = totalB;
   newBCinflux_rate = 0.0;     // per hour
   smooth_stopBCinflux = -1;   // time [hr] of
   min_seeder_dist = -1;
   max_seeder_dist = -1;
   // Diffusionkonstanten berechnet aus Stokes und Blut-Viskositaet
   D_CB = 5.;                   // microm^2/h
   v_CB = 2.;                   // i.e. standard is to use diffusion for definition of movement
   v_CB_width = -1;             // fixed v_CB value
   CB_v_modi = 1;               // i.e. random velocity state changes
   CB_n_v_states = 1;           // i.e. only one velocity state
   v_CB_switch_deltat = -1.0;   // only one v_state
   v_CB_factor = 1.0;           // one velocity state only, i.e. velocity factor = 1.0
   v_CB_cytosol = -1.0;         // use CB_D_cytosol for surface tension
   distance_tolerance = 0.3;
   half_tolerance_deformation = 0.01;   // always minimum tolerance
   mutation = 0.5;
   mutation_after_tc = -1.0;
   mutation_after_dec_tc = -1.0;
   mutation_affinity_exponent = 0.0;
   CBreceptor_use = 0;
   CBreceptor_dissociation = 1.0;
   CBreceptor_binding = 1.0;
   CBreceptor_activation = 0.5;
   CBreceptor_total = 1;
   CB_D_cytosol = D_CB;
   CB_elongation = 1.;
   CB_K_elongation = 2.;
   CB_smoothmove = 1.;
   CB_persistence = 2.;                 // no persistence
   CB_maxvolume4differ2CC = 1.0;        // no restriction for proliferation from volume
   CB_fixed_times_of_divisions = -1.;   // # of cell cycle times before susceptible to
                                        // differentiation
   fixed_time_of_divisions_mode = 0;    // calculate according duration of monoclonal expansion
   CB_fixed_times_of_divisions_in_expansion = 12.0;
   CB_dt_G0 = 0.0;
   CB_dt_G1 = 0.0;
   CB_dt_G2 = 0.0;
   CB_dt_S = 0.0;
   CB_dt_M = 0.0;
   CB_dtphase_width = 0.0;
   transmit_CC_delay_to_CB_cycle = false;
   t_inject_BrdU = -1;
   deltat_inject_BrdU = -1;
   n_inject_BrdU = 0;
   BrdU_detection_threshold = 0.05;
   CB2OUT_prob = -1.;
   exit2tz = 0;   // random walk of output cells
   retain_ag = false;
   divide_ag_asymmetric = 0.;
   asymmetric_polarity_index = 1.0;
   smooth_PI = 0.0;
   ag_deleted_in_fresh_CC = true;
   ag_loaded_CB_diff2output = false;
   ag_loaded_CC_directly2TFH = false;
   ag_loaded_CB_stop_mutation = false;
   BC_ag_preloaded = -1;

   // blast2:
   total_blast2 = 0;
   D_blast2 = 30;   // microns^2/h
   blast2_radius = 5;   // microns
   dx_blast2 = 1;             // microns (1 means only on next neighbors)
   blast2_proliferate = 10;   // hours
   blast2_grow = 1;           // hours
   blast2_distance_tolerance = 0.3;
   blast2_half_tolerance_deformation = 0.01;

   // CCs
   D_CC = -15.;            // microm^2/h
   CXCR4down = -1.;        // no time based down regulation of CXCL13 sensitivity
   CXCR5down = -1.;        // no time based down regulation of CXCL13 sensitivity
   CC_FDC_selection = 1;   // CC have to see FDC!
   collectFDCsignals = false;
   collectFDCperiod = 0.;
   prob2kill_noFDCcontactBCs = -1.0;
   present_specific_ag2TC = 0; // provide 0=total, 1=max-value, 2=TC-specific-value to TC
   v_CC = 5.;
   v_CC_width = -1;
   v_OUT = 5.;
   v_OUT_width = -1;
   CC_v_modi = 1;
   CC_n_v_states = 1;
   v_CC_factor = 1.;
   v_CC_switch_deltat = -1.;
   CC_persistence = 2.;         // no persistence
   OUT_persistence = 2.;        // no persistence
   use_ab_dynamics = 0;         // compatible to versions before hyphasma7.06.4
   initial_ab_affinity = -1.;   // use average seeder cell affinity
   CC_ICAM_delay = -1.;
   multipleTFHcontacts = false;
   negativeTCselection = true;
   ignore_apoptotic_CC = false;
   CC_apoptotic_motility_mode = 1;
   p_apo_randomwalk = 0.;
   reset_antigen_after_collection = -1;
   ignore_affinity = -1.;

   do_TC_division = 0;     // boolean
   TC_doubling = 24.0;     // time [hours] to be converted in a rate
   TC_meancycle = 10.0;    // hours
   TC_cyclewidth = 0.5;    // fraction [%]
   TC_Ndivisions = 1;      // number
   dx_TC = 0.;             // distance in microns
   tc_search_duration_mode = 0;
   tc_search_duration_fixed = 3.0;
   tc_search_duration_per_FDCcontact = 0.75;
   mode_of_setting_TC_time = 0; // fixed TC_time
   TC_time = 2.1;          // hours
   TC_time_width = 0.5;    // hours
   TC_rescue_time = 2.0;   // hours
   BCstaysonTCbyTCtime = 0;

   ///§§§ Philippe
   time_tc_selection_block = 24;       // (hours) -> no blocking    382
   factor_tc_selection_block = 1.0;    //  -> no blocking           383
   time_DND_block = 24;                //                           384
   factor_DND_block = 1.0;             //  -> no blocking           385
   factor_founder_div_block = 1.0;     //                           388
   ///§§§
   mode_tc_selection_block = 0;        //  -> no blocking           389

   pMHC_dependent_division = false;
   signal_dependent_number_of_divisions = false;
   pMHC_dependent_P_standard = 2.0;
   pMHC_dependent_P_max = 6.0;
   pMHC_dependent_P_min = 1.0;
   pMHC_dependent_K = 8.0;
   TFHsignal_dependent_K = 1.5;
   pMHC_dependent_nHill = 1.0;
   pMHC_dependent_pMHC_of_2divisions = 2.0;
   TFHsignal_of_P0divisions = 0.75;

   ICOSL_dependent_Tfh_signals = false;
   ICOSL_memory = false;
   ICOSL_upregulation_mode = 0;
   ICOSL_upregulation_time = 0.;

   // 2017-01-05 New version from git
   dT_FoxO = 1.0; // FoxO down at zero is getting back to 1=100% in this time
   nFoxO = 1.0;
   KFoxO = 3.0;
   dT_mTORC1 = 12.0; // mTORC1 reaches 1.0 in this time


   // class switch
   do_switch_classes = 0;
   // the following generates a switch matrix with 0 everywhere and 1 in the diagonal
   int switch_counter = 0;
   for (int i = 0; i < switch_dimension; i++) {
      switch_matrix[i] = 0;
      if (switch_counter == 0) {
         switch_matrix[i] = 1;
         switch_counter = int (nIg_classes) + 1;
      }
      --switch_counter;
   }
   IgE_BCRlevel = 0.3;

   ///§§§ Philippe 26/03/2017
   IgG_BCRlevel = 1.0;
   IgG_factor_cellcycle = 1.0;     // 374
   IgG_factor_divisions = 1.0;     // 375
   for(int i = 1; i < nIg_classes; ++i){
      Founder_IgX[i] = 0;
   }
   Founder_IgX[IgM] = 1.0;         // 376
   decay_proba_switch = 0;       // 377 for all probabilities of switching
   IgG_factor_leaving = 1.0;       // 378
   Affinity_threshold_IgG = 0.0;   // 379
   Int_Antigen_threshold_IgG = 0.0;// 380
   tc_help_IgG = 1.0;              // 381
   stddev_initial_divisions = 0; // 386  if != 0, then apply it
   stddev_DND = 0;               // 387  if != 0, then apply it



   CC_IgE_prob_CXCR5down = 0.;
   IgE_factor_cellcycle = 1.0;
   IgE_factor_divisions = 1.0;

   // output
   output = 0.2;
   output_DEC = output;

   TCell = 1.;
   FDCsignalling = 1.;

   // T cells
   TC_radius = 3.;        // microns
   v_TC = 8.;             // microns/min
   v_TC_width = -1;
   v_TC_CC = 0.;          // microns/min
   TC_persistence = 2.;   // min
   TC_CC_selection = 0;   // no active selection by TC
   north_weight = 0.;     // TC do not tend to go north in the reaction volume

   // FDCs
   for (int i = 0; i < 100; i++) {
      posFDC[i] = -1;
   }
   FDClength = 10;   // microm
   ag_per_FDC = -260.;
   ag_saturation_FDC = 20.;
   ag_distribution_mode = 0;
   ag_detection_mode = 0;

   // antibodies
   ag_threshold = 1.e-08;   // Mol
   ic_k_on = 1.e06;         // /(Mol s)
   ic_k_off = 1.e-03;       // /s

   antibodies_resolution = 0;        // 209;
   antibodies_production = 1.e-17;   // 210; in mol per hour and cell
   antibodies_degradation = 30;   // 214; days
   k_ic_exp_min = 5.5;   // 211; exponent of value in 1/Mol
   k_ic_exp_max = 10.5;   // 212; exponent of value in 1/Mol
   pm_differentiation_time = 24.;   // 213; hours
   N_GC = 1000;   // 215;
   V_blood = 0.01;   // 216; liter
   inject_antibody = 0.;
   injected_antibody_affinity = 8.5;
   inject_antibody_ASindex = -1;
   inject_antibody_time = 96.;

   mk_SEMA4D = -1.0e-08;   // Mol/h FDC

   // Nutrients
   use_glucose = -95.e-18;   // mol/sec
   use_oxygen = -20.e-18;    // mol/sec
   use_glucose_pro = -95.e-18;
   use_oxygen_pro = -20.e-18;
   bound_glucose = 0.;
   bound_oxygen = 0.;
   D_glucose_H2O = 4.146e+04;     // microns^2/min
   D_oxygen_H2O = 1.464e+05;      // microns^2/min
   D_glucose = 6.3e+03;           // microns^2/min
   D_oxygen = 1.05e+05;           // microns^2/min
   critical_nutrient = 2.5e-08;   // Mol^2
   fix_glucose_gradient = false;
   fix_glucose_gradient_min = 2.0;       // mM
   fix_glucose_gradient_max = 12.0;      // mM
   const_dynamic_glucose_field = 10.0;   // in glucose resting value

   // Cell tracking:
   tALL = 0;
   tCB = 0;
   tCC = 0;
   tOUT = 0;
   tTC = 0;
   tBETA = 0;
   trackfrom = 0.0;
   trackuntil = 1.0;           // hours
   track_delta_t = 1. / 60.;   // hours (1 minute)
   v_resolution = 100;
   delta_v = 2.;   // micron/min
   alpha_resolution = 36;
   delta_alpha = 5.;   // grad
   s_resolution = 100;
   delta_s = 0.1;   // micron/min

   // BETA-cells
   BETA_Nini = 0;
   BETA_radius = 5;          // microns
   BETA_proliferate = -1;    // hours (-1: no proliferation)
   BETA_max_pro = 1;         // microns (1=next neighbours)
   BETA_grow = 1;            // hours
   BETA_shrink = 1;          // hours
   BETA_max_adhesion = 0;    // maximum adhesion force in % of full stickness
   BETA_persistence = 0.1;   // min
   BETA_v = 0.1;             // micron/min
   BETA_v_modi = 1;          // random velocity mode
   BETA_n_v_states = 1;      // one state only
   BETA_v_factor = 1;
   BETA_v_switch_deltat = -1;   // min
   BETA_v_cytosol = 0.1;   // micron/min
   BETA_elongation = 1;
   BETA_K_elongation = 2;
   BETA_distance_tolerance = 0.2;
   BETA_half_tolerance_deformation = 0.2;
   BETA_smoothmove = 1;
}
void Werte::ini2d() {
   // Achtung!!! Die Standardwerte sind von res0x29a!
   // Mit "### statt xxx" sind Verbesserungsvorschlaege gemacht.
   // Standardwerte unabhaengig von der Dimension festlegen:
   inialld();

   // time steps
   deltat = 0.01;
   ToFileStep = 2400;

   // Antigene:
   takeA[0] = 3333;

   // Biological parameters
   proliferate = 9.;   // h
   grow = 3.;          // h;   (2 0.9 proliferate dx^2) / pi r_CB^2 = 3h
   shrink = 2.;        // h
   tolight = 2.8;      // ln(2)/6 h
   smooth_differentiation = false;
   smooth_differentiation_time = 3.;   // hours
   smooth_dif2out = false;
   smooth_dif2out_time = 3.;      // hours
   apoptosis = 10.;               // h
   apoptosis4FDCselected = -1.;   // infinite
   p_macrophage = 96.3;   // for necrotic CB (in h) (value from Gernot's model 138, *ln(2)=96 ???
                          // mit ihm checken #####)
   macrophage = 0.01;     // h
   selection = 2.;        // h
   ccdiff = 7.;           // h
   ccdiff_delay = -1.;    // h
   ccdiff_delay_DEC = -1.;             // h
   final_differentiation_rate = -1.;   // hr
   CC_test_delay = 0.1;                // h
   mksignal = 8.5;                     // /h FDC

   totalTC = 10;   // 50

   // FDC:
   FDCnumber = 20;
   FDCvesicle = 0;
   // Take 2/3 of GC volume for FDC network
   FDCnetwork = 0.6667;
   FDCtransparent = 1;
   // Position von FDCs in den oberen 70% der Flaeche (2D)
   posFDC[0] = 150;
   posFDC[1] = 200;
   posFDC[2] = 298;
   posFDC[3] = 324;
   posFDC[4] = 394;
   posFDC[5] = 511;
   posFDC[6] = 650;
   posFDC[7] = 692;
   posFDC[8] = 703;
   posFDC[9] = 751;
   posFDC[10] = 806;
   posFDC[11] = 842;
   posFDC[12] = 950;
   posFDC[13] = 969;
   posFDC[14] = 1030;
   posFDC[15] = 1050;
   posFDC[16] = 1093;
   posFDC[17] = 1150;
   posFDC[18] = 1219;
   posFDC[19] = 1250;
   mkCXCL12 = 0.0;   // Mol/(h FDC)
   mkCXCL13 = 0.0;   // Mol/(h FDC)

   // CBs:
   CB_radius = 4.5;   // microm
                      // maximal CB-proliferation distance in microm
   // dx_CB=8.*CB_radius;
   dx_CB = 0.;
   // Position der seeder CBs passend zu obigen FDCs
   // posCB[0]=385;
   // posCB[1]=892;
   // posCB[2]=1190;

   // OUT
   mk_ab = -3.0e-08;   // number per hour

   // For lattice GCx
   DimSpace = 2;
   // Radius of GC in microm
   GC_radius = 220.;
   // Separately in each dimension:
   gridsize[0] = 220.;
   gridsize[1] = 220.;
   gridsize[2] = 0.;
   // lattice constant in microm
   dx = 10.;
   // lattice constant for signal grid in microm (<0. : =dx)
   dx_signal = -1.;
}
void Werte::ini3d() {
   inialld();
   // Timescales in h
   deltat = 0.05;
   ToFileStep = 480;

   // Biological parameters
   proliferate = 6.;    // h
   grow = 0.666667;     // h  // #### insert: 2 pi r_CB^3 p(incl ln(2)) / (3 0.9 dx^3)
   shrink = 0.666667;   // h
   tolight = 3.;        // h
   apoptosis = 6.;      // h
   apoptosis4FDCselected = -1.;
   macrophage = 0.01;       // h
   selection = 3.;          // h
   ccdiff = 3.5;            // h
   ccdiff_delay = -1;       // h
   ccdiff_delay_DEC = -1;   // h
   CC_test_delay = 2.;      // h
   mksignal = 27.4;         // /h FDC

   totalTC = 0;   // 500

   // FDC:
   FDCnumber = 172;
   FDCvesicle = 1;
   // Position von FDCs in den oberen 50% der Flaeche (3D)
   FDCnetwork = 0.5;
   FDCtransparent = 1;
   // Antigene:
   takeA[0] = 1170;

   // CBs
   CB_radius = 7.5;   // microm
   // maximal CB-proliferation distance in microm
   dx_CB = 8. * CB_radius;

   // OUT
   mk_ab = 0.;   // Mol per hour

   // For lattice GCx
   DimSpace = 3;
   // Radius of GC in microm
   GC_radius = 160.;
   // Separately in each dimension:
   gridsize[0] = 160.;
   gridsize[1] = 160.;
   gridsize[2] = 160.;
   // lattice constant in microm
   dx = 10.;
   // lattice constant for signal grid in microm (<0. : =dx)
   dx_signal = -1.;
}
Werte::Werte()
   : file_output(MAXDIM, 0),
   takeA(MAXDIM, 0, 1),
   initAntigenSeqs(MAXDIM, string("-1")),
   initBCRSeqs(MAXDIM, string("-1")),
   initTCRSeqs(MAXDIM, string("-1")),
   fix_signals(MAXDIMSMALL, 0, 1),
   takeB(MAXDIM, 0, 1),
   posCB(MAXDIM, 0, 1),
   pos_blast2(MAXDIM, 0, 1),
   posFDC(100, 0, 1),
   BETA_pos(MAXDIM, 0, 1) {
   // Initialisiere standardmaessig fuer 2d Rechnungen
   ini2d();
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Parameter File-Formate:
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const char* Werte::kenntext(int nummer) {
   switch (nummer) {
      // free: 402+

   /// Philippe
      // General
   case 390: {
       return "Use signal files with following prefix ('none' if none)";
   }
      case 1: {
         return "System: 0=Unix; 1=Windows";
      }
      break;

      case 2: {
         return "Zufallsgenerator initialisieren";
      }
      break;

      case 3: {
         return "Extra Zufallsgenerator am Anfang";
      }
      break;

      case 4: {
         return "File-Output-Restrictions";
      }
      break;

      case 5: {
         return "Zufallsgenerator post-proliferation initialisieren";
      }
      break;

      case 6: {
         return "Angabe von Zeiten =1 (nicht Raten =0)";
      }
      break;

      case 7: {
         return "Was wird rausgeschrieben (0=alles; 1=wenig)";
      }
      break;

      case 8: {
         return "Hebe Ki67 hervor? (0=nein; 1=ja)";
      }
      break;

      case 9: {
         return "Do checks? (0=no; 1=few; 2=yes)";
      }
      break;

      case 124: {
         return "Mode of spatial output (0=GC 1=tumour ...)";
      }
      break;

      case 172: {
         return "Initial array dimension of CB";
      }
      break;

      case 173: {
         return "Initial array dimension of CC";
      }
      break;

      case 174: {
         return "Initial array dimension of TC";
      }
      break;

      case 175: {
         return "Initial array dimension of OUT";
      }
      break;

      case 176: {
         return "Initial array dimension of FDC";
      }
      break;

      case 177: {
         return "Initial array dimension of STROMA";
      }
      break;

      case 178: {
         return "Initial array dimension of BETA";
      }
      break;

      // shape space and antigen
      case 10: {
         return "Use metric (1:N_mutation; 2:euclidian)";
      }
      break;

      case 11: {
         return "Dimension of Shapespace (<11)";
      }
      break;

      case 12: {
         return "Number of antibody types";
      }
      break;

      case 13: {
         return "Number per Dimension (Fixed,int-type)";
      }
      break;

      case 14: {
         return "Total initial Number of presented Antigen Epitops (int-type)";
      }
      break;

      case 15: {
         return "Number of initial Anitgen Peaks in its Shapespace (int-type)";
      }
      break;

      case 16: {
         return "Fix Epitop presentation (max 10 values)";
      }
      break;

      case 323: {
         return "Fraction of Ags (non-fixed Ag enter with same fraction)";
      }
      break;

      case 18: {
         return "Width of gaussian affinity weight function";
      }
      break;

      case 121: {
         return "Amplitude of Gauss affinity weight function (0<a<=1)";
      }
      break;

      // Sequence shape
      case 349: {
         return "Use logarithmic affinity";
      }
      break;

      case 307: {
         return "Number of Initial Antigen sequences (int-type)";
      }
      break;

      case 308: {
         return "Maximum Hamming distance between antigens";
      }
      break;

      case 309: {
         return "Minimum Hamming distance between antigens";
      }
      break;

      case 310: {
         return "Fix Antigen Sequence presentation (max 1000 values)";
      }
      break;

      case 311: {
         return "Fix initial Repertoire distribution (max 1000 values)";
      }
      break;

      case 312: {
         return "Maximum Initial Hamming distance between BCRs";
      }
      break;

      case 313: {
         return "Minimum affinity of initial BCRs to Antigens";
      }
      break;

      case 314: {
         return "Maximum affinity of initial BCRs to Antigens";
      }
      break;

      case 315: {
         return "Length of sequences";
      }
      break;

      case 348: {
         return "Mutation proba per base per division (sequence space only)";
      }
      break;

      case 316: {
         return "Specifity of sequences affinity (double R)";
      }
      break;

      case 317: {
         return "Fix initial Repertoire distribution for T cells";
      }
      break;

      case 318: {
         return "Maximum Initial Hamming distance between TCRs";
      }
      break;

      case 319: {
         return "Minimum affinity of initial TCRs to Antigens";
      }
      break;

      case 320: {
         return "Maximum affinity of initial TCRs to Antigens";
      }
      break;

      case 321: {
         return "Use sequence space (1/0)";
      }
      break;

      case 326: {
         return "Optimum affinity cluster size (affinity doesn't increase beyond it)";
      }
      break;

      case 327: {
         return
            "Type of affinity function : 0= standard (saham's)  1= normalized to max_affinity_cluster 2= "
            "sliding windows of size max_affinity_cluster";
      }
      break;

      // space lattice in general
      case 20: {
         return "Dimension of lattice";
      }
      break;

      case 21: {
         return "Lattice constant of space grid";
      }
      break;

      case 123: {
         return "Lattice constant of signal grid (-1. for =space grid) (in microm)";
      }
      break;

      case 22: {
         return "Radius of GC (microm)";
      }
      break;

      case 23: {
         return "Shape of reaction volume (0=sphere; 1=cube)";
      }
      break;

      case 110: {
         return "Obstacles: 0=no; 1=wall; 2=wall+slits; 3=random";
      }
      break;

      case 111: {
         return "Position of wall (% of volume height)";
      }
      break;

      case 112: {
         return "Width of wall (in points)";
      }
      break;

      case 113: {
         return "Number of slits";
      }
      break;

      case 114: {
         return "Width of slits (in points)";
      }
      break;

      case 115: {
         return "Density of random obstacles (% of reaction volume)";
      }
      break;

      case 116: {
         return "Clustering of obstacles (% of random obstacles per obstacle)";
      }
      break;

      case 200: {
         return "length of x-axis for vol-shape>1 (micron) (2D: hor; 3D: depth)";
      }
      break;

      case 201: {
         return "length of y-axis for vol-shape>1 (micron) (2D: vert; 3D: hor)";
      }
      break;

      case 202: {
         return "length of z-axis for vol-shape>1 (micron) (3D: vert)";
      }
      break;

      // Cells in general
      case 117: {
         return "Maximum weight of chemokine to random polarity";
      }
      break;

      case 118: {
         return "Steepness of weight reduction with chemokine gradient";
      }
      break;

      case 119: {
         return "Chemokine gradient of half weight";
      }
      break;


       /// Philippe 2017-10-22 T cell specific chemotaxis
       // Cells in general
       case 393: {
          return "T cell specific Maximum weight of chemokine to random polarity";
       }
       break;

       case 394: {
          return "T cell specific Steepness of weight reduction with chemokine gradient";
       }
       break;

       case 395: {
          return "T cell specific Chemokine gradient of half weight";
       }
       break;
       case 396: {
           return "use_specific_tc_chemotaxis";
       }
       break;

   /// Philippe 26-10-2017
   case 397: {
       return "proba_TC_CC_interaction";
   }
   break;



      case 120: {
         return "Duration for establishing adhesion binding (min)";
      }
      break;

      case 122: {
         return "Macrophagocytosis of necrotic cells";
      }
      break;

      case 125: {
         return "Allow exchange of cells during movement (0,1)";
      }
      break;

      case 128: {
         return "Set polarity with specific distributions (0:no,1:gauss,2:cyster)";
      }
      break;

      // Centroblasts in shape space and space
      case 30: {
         return "Total Number of initial B-Cells";
      }
      break;

      case 99: {
         return "Total Number of initial B-Cells types";
      }
      break;

      case 297: {
         return "Rate of new BC flux into the GC (0=none) [# of cells/hr]";
      }
      break;

      case 298: {
         return "Smoothness width of switch-off of new BC influx (-1=no) [hr]";
      }
      break;

      case 305: {
         return "Seeder cells have a SS distance to antigen of more than (-1: no limit) [metric]";
      }
      break;

      case 306: {
         return "Seeder cells have a SS distance to antigen of less than (-1: no limit) [metric]";
      }
      break;

      case 31: {
         return "Centroblast radius (microm)";
      }
      break;

      case 32: {
         return "Maximal distance for CB proliferation (microm)";
      }
      break;

      case 33: {
         return "Fix initial Centroblast distribution (max 10 values)";
      }
      break;

      case 34: {
         return "Fix initial Centroblast position (max 10 values)";
      }
      break;

      case 35: {
         return "CB use receptors for signal-induced differentiation (0=no 1=ss 2=dyn)";
      }
      break;

      case 36: {
         return "Bound-receptor-threshold for CB-differentiation (fraction)";
      }
      break;

      case 37: {
         return "Dissociation constant for CB-differentiation signal (quanta or h)";
      }
      break;

      case 38: {
         return "Binding constant for CB-differentiation signal (h)";
      }
      break;

      case 39: {
         return "Total number of receptors on CB";
      }
      break;

      case 168: {
         return "Duration of CXCR4 expression in hours (-1 for no downregulation)";
      }
      break;

      case 100: {
         return "Surface tension as cytosol diffusion in microns^2/h";
      }
      break;

      case 101: {
         return "CB elongation during active move (% of CB radius)";
      }
      break;

      case 130: {
         return "CB elongation for half maximum reshaping force (in CB radii)";
      }
      break;

      case 102: {
         return "Smoothness of CB movement (# of timesteps for dx=a move, default 1)";
      }
      break;

      case 103: {
         return "Persistence of CB polarity (time gap in min)";
      }
      break;

      case 104: {
         return "Mean cell velocity (in microns/min)";
      }
      break;

      case 356: {
         return "Width of target CB velocity (percentage of mean; -1 for fixed)";
      }
      break;

      case 105: {
         return "Modus of velocity (1=random; 2=1+weights; 3=polarity-coupled; 4=adhesion)";
      }
      break;

      case 131: {
         return "Number of velocity states";
      }
      break;

      case 106: {
         return "Factor of velocity reduction from v_state max to min";
      }
      break;

      case 107: {
         return "Maximum adhesion force of CB (% of full stickness)";
      }
      break;

      case 108: {
         return "Surface tension as fragment velocity in microns/min";
      }
      break;

      case 109: {
         return "Persistence of v-states in min";
      }
      break;

      case 50: {
         return "Diffusion constant of CB in microm^2/h";
      }
      break;

      case 53: {
         return "Tolerance for distance to barycenter for movement to target point in %";
      }
      break;

      case 57: {
         return "Half tolerance at deformation of (0.01 - infty)";
      }
      break;

      case 19: {
         return "Fraction of volume up to which proliferation is possible (1=no restriction)";
      }
      break;

      case 70: {
         return "Proliferation rate";
      }
      break;

      case 227: {
         return "Number of required cell cycles before differentiation (cell cycle times)";
      }
      break;

      case 275: {
         return
            "Mode of division number in expansion phase (0: calculate; 1: as later; 2: code specified)";
      }
      break;

      case 293: {
         return
            "Fixed number of divisions in expansion phase (chose 3 in the mode of division above)";
      }
      break;

      case 71: {
         return "Mutation probability:";
      }
      break;

      case 235: {
         return "Mutation probability after TC selection (-1: for unchanged; or [0..1])";
      }
      break;

      case 245: {
         return "Mutation probability after DEC205-competent TC selection (-1: for same; [0..1])";
      }
      break;

      case 236: {
         return "Affinity-dependent mutation upon TC contact (<=0: none; affinity-exponent)";
      }
      break;

      case 72: {
         return "Rate for differentiation of centroblasts to centrocytes";
      }
      break;

      case 203: {
         return "Smooth onset of CB differentiation (1=yes, 0=no)";
      }
      break;

      case 250: {
         return "Width of smooth onset of CB differentiation (hours)";
      }
      break;

      case 249: {
         return "Probability of CB to OUT differentiation (<=0 for none)";
      }
      break;

      case 264: {
         return "Exit of OUT cells towards T zone (1: yes; 0: random walk)";
      }
      break;

      case 256: {
         return "Duration of CB cell cycle phase G1 (hours)";
      }
      break;

      case 257: {
         return "Duration of CB cell cycle phase S (hours)";
      }
      break;

      case 258: {
         return "Duration of CB cell cycle phase G2 (hours)";
      }
      break;

      case 259: {
         return "Duration of CB cell cycle phase M (hours)";
      }
      break;

      case 260: {
         return "Duration of CB cell cycle phase G0 (hours)";
      }
      break;

      case 261: {
         return "Width of Gaussian variation of phases (fraction of average duration)";
      }
      break;

      case 263: {
         return "Transfer CC delay of differentiation to CB cycle time (1: yes; 0: no)";
      }
      break;

      case 352: {
         return "Inject BrdU at time (hours) [-1 for none]";
      }
      break;

      case 353: {
         return "Repeat BrdU injection every (time interval in hours) [-1 for none]";
      }
      break;

      case 354: {
         return "Number of BrdU injections (set Repeat BrdU ... >0)";
      }
      break;

      case 355: {
         return "Detection threshold of BrdU for read-out (% of max)";
      }
      break;

      case 267: {
         return "Keep loaded antigen after BC selection for recycling (1:yes; 0: no)";
      }
      break;

      case 268: {
         return "Divide antigen asymmetrically on daughter BC (1: 100%; 0: no)";
      }
      break;

      case 273: {
         return "Distribute antigen asymmetric (0.5 < polarity index < 1.0)";
      }
      break;

      case 276: {
         return "Width of smooth distribution around the polarity index [% of PI]";
      }
      break;

      case 269: {
         return "Antigen-retaining daughter BCs differentiate to output (1: yes; 0: no)";
      }
      break;

      case 274: {
         return "Retained antigen is deleted in fresh CC (1:yes; 0:no)";
      }
      break;

      case 270: {
         return
            "CC differentiated from antigen-loaded CB directly interact with TFH (1: yes; 0: no)";
      }
      break;

      case 271: {
         return "Antigen pre-loaded BCs suppress mutation (1:yes; 0: as others)";
      }
      break;

      case 277: {
         return "BC are pre-loaded with antigen at start (-1:no; # of antigen portions)";
      }
      break;

      // TC
      case 24: {
         return "Number of initial TC";
      }
      break;

      case 25: {
         return "Radius of TC in microns";
      }
      break;

      case 26: {
         return "Velocity of TC in microns/min";
      }
      break;

      case 359: {
         return "Width of target TC velocity (percentage of mean; -1 for fixed)";
      }
      break;

      case 27: {
         return "Velocity of TC in interaction with CC in microns/min";
      }
      break;

      case 28: {
         return "Persistence of polarity of TC in min";
      }
      break;

      case 242: {
         return "Tendency for TC to stay in the LZ (weight in 0..1)";
      }
      break;

      case 299: {
         return "TC divide (0=no; 1=yes)";
      }
      break;

      case 300: {
         return "TC doubling time (hours, <0 for mechanistic triggering of division)";
      }
      break;

      case 301: {
         return "Duration of TC cell cycle (hours)";
      }
      break;

      case 302: {
         return "Width of Gaussian variation of cell cycle (fraction of average duration)";
      }
      break;

      case 303: {
         return "Number of required cell cycles before return to quiescence (cell cycles)";
      }
      break;

      case 304: {
         return "Maximal distance for TC proliferation (microm)";
      }
      break;

      // FDCs
      case 40: {
         return "Total number of FDCs";
      }
      break;

      case 41: {
         return "Fix FDC position (max 100 values)";
      }
      break;

      case 42: {
         return "Length of FDC arms (microm)";
      }
      break;

      case 43: {
         return "Percentage of FDCnetwork in GC volume";
      }
      break;

      case 44: {
         return "FDC dendrites treated transparent (yes=1; no=0)";
      }
      break;

      case 68: {
         return "Threshold Ag-concentration for binding CC (in Mol)";
      }
      break;

      case 29: {
         return "Presented Ag per FDC in units of threshold (-1: no consumption)";
      }
      break;

      case 46: {
         return
            "Ag saturation per FDC-fragment in units of threshold (1: constant binding probability)";
      }
      break;

      case 324: {
         return
            "Ag distribution per FDC fragment (0: as ag_fraction; 1: 1 random Ag with prob=ag_fraction)";
      }
      break;

      case 325: {
         return "Ag detection on FDC fragment (0: highest affinity Ag; 1: highest Ag amount)";
      }
      break;

      // antibodies and immune complexes
      case 83: {
         return "Diffusion constant of soluble antibodies in microm^2/h";
      }
      break;

      case 47: {
         return "Antibody production rate in Mol per hour and output cell";
      }
      break;

      case 48: {
         return "Boundary for antibodies";
      }
      break;

      case 126: {
         return "Ag presenting Ab adapt to PC-produced Ab-quality (1=yes, 0=no)";
      }
      break;

      case 208: {
         return "Initial soluble Ab-affinity for antigen (0..1, -1 for seeder-affinity)";
      }
      break;

      case 67: {
         return "k_off for dissociation of immune complex (in /s)";
      }
      break;

      case 66: {
         return "k_on for building immune complex (in /(Mol s); -1 for no binding)";
      }
      break;

      case 209: {
         return "Resolution of systemic antibodies in affinity bins (0=no systemic antibodies)";
      }
      break;        // antibodies_resolution=0

      case 210: {
         return
            "Production rate of systemic antibodies by output population (mol per hour and cell)";
      }
      break;        // antibodies_production=1.e-17

      case 214: {
         return "Antibody degradation half time (days)";
      }
      break;        // antibodies_degradation=30 days

      case 211: {
         return "Exponent of lowest affinity dissociation constant of immune complex (1/Mol)";
      }
      break;        // k_ic_exp_min=5.5

      case 212: {
         return "Exponent of highest affinity dissociation constant of immune complex (1/Mol)";
      }
      break;        // k_ic_exp_max=10.5

      case 213: {
         return "Half time of plasma cell differentiation to antibody producing cell (hours)";
      }
      break;        // pm_differentiation_time=log(2.)/24.;

      case 215: {
         return "Number of GC generating antibody producing plasma cells";
      }
      break;        // N_GC=1000

      case 216: {
         return "Blood volume (l)";
      }
      break;        // V_blood=0.01 l

      case 217: {
         return "Inject antibodies: concentration (mol/l; <=0 for none)";
      }
      break;

      case 218: {
         return "Inject antibodies: exponent of affinity in l/mol";
      }
      break;

      case 322: {
         return "Inject antibodies: index in affinityspace (-1 to use exponent of affinity):";
      }
      break;

      case 234: {
         return "Inject antibodies: time of administration (hours)";
      }
      break;

      // Signals
      case 169: {
         return "Take signals from file and fix dynamics (0 or 1)";
      }
      break;

      case 54: {
         return "Signal diffusion mode (0:QUANTA; 1: EULER; 2: ADI)";
      }
      break;

      case 56: {
         return "Objects are transparent for signals (1:yes; 0: no)";
      }
      break;

      case 45: {
         return "FDC signal production mode (1: Vesikel; 0: continuous)";
      }
      break;

      case 52: {
         return "Diffusion constant of differentiation signal in microm^2/h";
      }
      break;

      case 79: {
         return "Rate of diff2CC signal production";
      }
      break;

      case 55: {
         return "Boundary for Differ2CC-Signal";
      }
      break;

      case 69: {
         return "Diffusion constant of chemotaxis signal in microm^2/h";
      }
      break;

      case 59: {
         return "Rate of CXCL13 production";
      }
      break;

      case 58: {
         return "Boundary for CXCL13";
      }
      break;

      case 171: {
         return "Critical CXCL13 concentration for desensitisation [Mol]";
      }
      break;

      case 238: {
         return "Critical CXCL13 concentration for resensitisation [Mol, -1 for none]";
      }
      break;

      case 165: {
         return "Diffusion constant of CXCL12 (microm^2/hr)";
      }
      break;

      case 166: {
         return "Rate of CXCL12 production (Mol/(hr cell))";
      }
      break;

      case 167: {
         return "Boundary for CXCL12 (Mol)";
      }
      break;

      case 170: {
         return "Critical CXCL12 concentration for desensitisation [Mol]";
      }
      break;

      case 237: {
         return "Critical CXCL12 concentration for resensitisation [Mol, -1 for none]";
      }
      break;

      case 86: {
         return "Diffusion constant of soluble antigen in microm^2/h";
      }
      break;

      case 84: {
         return "Boundary for antigen";
      }
      break;

      case 87: {
         return "Semaphorin SEMA4D production in Mol per hour and FDC";
      }
      break;

      case 88: {
         return "Boundary of Semaphorin SEMA4D";
      }
      break;

      case 89: {
         return "Diffusion constant of Semapohrin SEMA4D in microm^2/h";
      }
      break;

      // Nutrients
      case 140: {
         return "Glucose consumption by viable cells (mol/sec)";
      }
      break;

      case 141: {
         return "Oxygen consumption by viable cells (mol/sec)";
      }
      break;

      case 142: {
         return "Glucose consumption by proliferating cells (mol/sec)";
      }
      break;

      case 143: {
         return "Oxygen consumption by proliferating cells (mol/sec)";
      }
      break;

      case 144: {
         return "Boundary for glucose (Mol)";
      }
      break;

      case 145: {
         return "Boundary for oxygen (Mol)";
      }
      break;

      case 146: {
         return "Diffusion constant of soluble glucose (microns^2/min)";
      }
      break;

      case 147: {
         return "Diffusion constant of soluble oxygen (microns^2/min)";
      }
      break;

      case 148: {
         return "Diffusion constant of glucose in tissue (microns^2/min)";
      }
      break;

      case 149: {
         return "Diffusion constant of oxygen in tissue (microns^2/min)";
      }
      break;

      case 150: {
         return "Critical nutrient product concentration for necrosis (Mol^2)";
      }
      break;

      // Centrocytes
      case 51: {
         return "Diffusion constant of CC/out in microm^2/h";
      }
      break;

      case 132: {
         return "Persistence of CC polarity (time gap in min)";
      }
      break;

      case 133: {
         return "Mean CC velocity (in microns/min)";
      }
      break;

      case 357: {
         return "Width of target CC velocity (percentage of mean; -1 for fixed)";
      }
      break;

      case 134: {
         return "Modus of CC velocity (1=random; 2=1+weights; 3=polarity-coupled; 4=adhesion)";
      }
      break;

      case 135: {
         return "Number of CC velocity states";
      }
      break;

      case 136: {
         return "Factor of CC velocity reduction from v_state max to min";
      }
      break;

      case 137: {
         return "Persistence of v-states in min";
      }
      break;

      case 129: {
         return "Duration of CXCR5 expression in hours (-1 for until FDC-contact)";
      }
      break;

      case 283: {
         return "Motility mode of apoptotic CC (0=no; 1=CXCL13; 2=random; 3=CXCL12)";
      }
      break;

      case 284: {
         return
            "Half life of CXCL12 sensitivity of apoptotic CC before random walk (hours, -1:none)";
      }
      break;

      case 85: {
         return "CC need to see FDC before selection (1=yes, 0=no)?";
      }
      break;

      case 253: {
         return "CC collect antigen by serial FDC encounters (1=yes, 0=no)";
      }
      break;

      case 254: {
         return "Duration of CC collection of antigen by serial FDC encounters (h)";
      }
      break;

      case 351: {
         // prob2kill_noFDCcontactBCs
         return "Probability to kill BCs that failed to collect antigen (-1: for always)";
      }
      break;

      case 17: {
         return "Present collected Ag to TC [0=all; 1=max-specific; 2=TC specificity]";
      }
      break;

      case 49: {
         return "Positive selection by TC-CC-interaction (0=no, 1=yes)";
      }
      break;

      case 360: {
         return "Derive CC selection from multiple contacts with Tfh (0=single, 1=multiple)";
      }
      break;

      case 233: {
         return "Negative selection by TC-CC-interaction (0=no, 1=yes)";
      }
      break;

      case 362: {
         return "Mode of setting TC-CC-interaction time (0=fixed; 1=Gauss; 2=affinity, see above)";
      }
      break;

      case 138: {
         return "Duration of TC-CC-interaction in hours";
      }
      break;

      case 361: {
         return "Width of TC-CC-interaction duration in hours (set mode = 1)";
      }
      break;

      case 139: {
         return "Minimum duration of TC-CC-polarisation for CC-rescue in hours";
      }
      break;

      case 350: {
         // BCstaysonTCbyTCtime
         return "BC stay in contact when signalling is above selection threshold (0=no, 1=yes)";
      }
      break;

       ///§§§ Philippe
   case 382: { // time_tc_selection_block
      return "time of blocking T cell selection (CD40) ";
   }
   break;
   case 383: { // factor_tc_selection_block
      return "factor blocking tc selection";
   }
   break;
   case 384: { // time_DND_block
      return "time of blocking DND (CD40)";
   }
   break;
   case 385: { // factor_DND_block
      return "Factor of CD40 blocking on DND (div = coeff * expected div)";
   }
   break;
   case 388: { //   factor_founder_div_block = 1.0;     //                           388
      return "Factor of CD40 blocking on founder division number (div = coeff * original number)";
   }
   break;
   case 389: { //   mode_tc_selection_block
      return "Mode of tc selection block (0=none, 1=tc_rescue_time, 2=tc_search_duration, 3=tc_interact_time)";
   }
   break;




      case 278: {
         return "Do class switch (0=no, 1=at TFH-selection, 2=at division)";
      }
      break;

      case 279: {
         return "Class switch probabilities (sum of each line must be 1)";
      }
      break;

      case 280: {
         return "BCR expression level of IgE BCs (percentage of normal)";
      }
      break;

      case 281: {
         return "Cell cylce time of IgE BCs (factor; 1:as normal)";
      }
      break;

      case 282: {
         return "Number of IgE BC divisions (factor; 1:as normal)";
      }
      break;

      case 285: {
         return "Failure prob of IgE-CC CXCR5 upregulation (0:no; >0:CXCR5 down; <0:CXCR4 up)";
      }
      break;

      case 286: {
         return "Use pMHC-dependent TFH-induced number of divisions (0,1)";
      }
      break;

      case 287: {
         return "pMHC-dependent division number Hill: maximum (P_max)";
      }
      break;

      case 288: {
         return "pMHC-dependent division number Hill: Hill-coefficient (n_P)";
      }
      break;

      case 289: {
         return "pMHC-dependent division number Hill: half (K_P; -1 for calc with A_0)";
      }
      break;

      case 290: {
         return "pMHC-dependent division number Hill: pMHC for 2 divisions (A_0)";
      }
      break;

      case 291: {
         return "pMHC-dependent division number Hill: minimum (P_min)";
      }
      break;

      case 292: {
         return "pMHC-dependent division number Hill: standard (P_0)";
      }
      break;

      case 294: {
         return
            "Reset amount of collected antigen after collection-phase (antigen portions; <=0: no reset)";
      }
      break;

      case 295: {
         return "Ignore affinity when binding antigen on FDCs (-1: no; mean binding probability)";
      }
      break;

      // Output
      case 265: {
         return "Persistence of OUT polarity (time gap in min)";
      }
      break;

      case 266: {
         return "Mean OUT velocity (in microns/min)";
      }
      break;

      case 358: {
         return "Width of target OUT velocity (percentage of mean; -1 for fixed)";
      }
      break;

      // Time and GC-phases
      case 60: {
         return "Time steps (in h)";
      }
      break;

      case 61: {
         return "Beginning time (negative tolight=0) (in h)";
      }
      break;

      case 62: {
         return "End time (in h)";
      }
      break;

      case 63: {
         return "Time steps between output";
      }
      break;

      case 296: {
         return "Stop new BC influx to GC (in h)";
      }
      break;

      case 64: {
         return "Start Mutation and Differentiation (in h)";
      }
      break;

      case 127: {
         return "Start Mutation independent of Differentiation (in h)";
      }
      break;

      case 65: {
         return "Start Output (in h)";
      }
      break;

      // Cell dynamics
      case 73: {
         return "Rate of positive selection at FDCs";
      }
      break;

      case 206: {
         return "Probability of signalling during functional CC-FDC contact";
      }
      break;

      case 74: {
         return "T-Cell-Selection probability";
      }
      break;

      case 75: {
         return "Rate of selected CC-differentiation";
      }
      break;

      case 262: {
         return "Delay of selected CC-differentiation (hours; -1: none)";
      }
      break;

      case 246: {
         return "Delay of selected CC-differentiation upon DEC-binding (hours; -1: none)";
      }
      break;

      case 76: {
         return "Output probability toward Antibody production";
      }
      break;

      case 244: {
         return "Output probability of DEC205-OVA activated BC toward Antibody production";
      }
      break;

      case 204: {
         return "Smooth onset of differentiation to output (1=yes, 0=no)";
      }
      break;

      case 251: {
         return "Width of smooth onset of differentiation to output (hours)";
      }
      break;

      case 255: {
         return "Delay of differentiation to output after selection (hours; -1: none)";
      }
      break;

      case 77: {
         return "Rate of apoptosis";
      }
      break;

      case 207: {
         return "Life time of FDC selected CC (hours, -1=infinite)";
      }
      break;

      case 252: {
         return "Ignore apoptotic CC in cell number analysis (0=no; 1=yes)";
      }
      break;

      case 78: {
         return "Rate of macrophage transport of dead cells";
      }
      break;

      case 80: {
         return "Rate of cell growth";
      }
      break;

      case 81: {
         return "Rate of cell shrinking";
      }
      break;

      case 82: {
         return "Time gap between affinity tests";
      }
      break;

      case 205: {
         return "Duration of motility suppression for each affinity test (-1=none, hr)";
      }
      break;

      // blast2 in space
      case 90: {
         return "Total Number of initial blast2-Cells";
      }
      break;

      case 91: {
         return "Blast2 radius (microm)";
      }
      break;

      case 92: {
         return "Maximal distance for blast2 proliferation (microm)";
      }
      break;

      case 93: {
         return "Diffusion constant of blast2 in microm^2/h";
      }
      break;

      case 94: {
         return "Fix initial blast2 position (max 10 values)";
      }
      break;

      case 95: {
         return "Blast2 proliferation";
      }
      break;

      case 96: {
         return "Blast2 cell growth";
      }
      break;

      case 97: {
         return "Blast2 tolerance for distance to barycenter for fragment diffusion %";
      }
      break;

      case 98: {
         return "Blast2 half tolerance at deformation of (0.01 - infty)";
      }
      break;

      // Cell tracking
      case 151: {
         return "Total number of tracked cells (0 = use numbers for each cell type)";
      }
      break;

      case 152: {
         return "Number of tracked CB";
      }
      break;

      case 153: {
         return "Number of tracked CC";
      }
      break;

      case 154: {
         return "Number of tracked output cells";
      }
      break;

      case 155: {
         return "Number of tracked TC";
      }
      break;

      case 156: {
         return "Tracking start time (hours)";
      }
      break;

      case 157: {
         return "Tracking end time (hours)";
      }
      break;

      case 158: {
         return "Tracking time interval (hours)";
      }
      break;

      case 159: {
         return "Number of intervals for the speed-distribution";
      }
      break;

      case 160: {
         return "Width of one interval for speed-distribution (micron/min)";
      }
      break;

      case 161: {
         return "Number of intervals for the elongation-distribution";
      }
      break;

      case 162: {
         return "Width of one interval for elongation-distribution (axis-ratio)";
      }
      break;

      case 163: {
         return "Number of intervals for the turning-angle-distribution";
      }
      break;

      case 164: {
         return "Width of one interval for turning-angle-distribution (degree)";
      }
      break;

      // Photoactivation
      case 219: {
         return "Do photoactivation of cells (0=no, 1=yes)";
      }
      break;

      case 220: {
         return "Start photoactivation at time (hours)";
      }
      break;

      case 221: {
         return "Coordinates of the photoactivation area (x,y,z in microns)";
      }
      break;

      case 222: {
         return "Size of the photoactivation area (x,y,z in microns)";
      }
      break;

      case 228: {
         return "Use DEC205 receptors on B cells (0=no, 1=yes)";
      }
      break;

      case 229: {
         return "Time of setting DEC205 receptors on B cells (hours)";
      }
      break;

      case 230: {
         return "Fraction of DEC205+/+ B cells";
      }
      break;

      case 231: {
         return "Inject anti-DEC205-OVA for BCR independent activation (0=no, 1=yes)";
      }
      break;

      case 232: {
         return "Time of anti-DEC205-OVA injection (hours)";
      }
      break;

      case 239: {
         return "Duration of anti-DEC205-OVA activity (hours, -1 for one selection)";
      }
      break;

      case 240: {
         return "Duration of TC-BC-contact upon DEC205OVA activation (hours, 0 for unchanged)";
      }
      break;

      case 241: {
         return "Increase number of TCs upon DEC205OVA activation (factor of normal)";
      }
      break;

      case 243: {
         return "Increase number of BC divisions upon DEC205OVA activation (factor of normal)";
      }
      break;

      case 247: {
         return "Induce CB differentiation upon DEC205OVA activation (0: no; 1: yes)";
      }
      break;

      case 248: {
         return "Probability of CB to OUT differentiation if DEC205OVA is bound (<=0 for none)";
      }
      break;

      case 272: {
         return "Antigen acquired via DEC205 is retained (0: no; 1: yes)";
      }
      break;

      // BETA-CELL
      case 179: {
         return "Total Number of initial betacells";
      }
      break;

      case 180: {
         return "betacell radius (microm)";
      }
      break;

      case 181: {
         return "betacell cell cycle (hours) (-1 for none)";
      }
      break;

      case 182: {
         return "Maximal distance for betacell proliferation (microm)";
      }
      break;

      case 183: {
         return "betacell inverse growth rate (hours)";
      }
      break;

      case 184: {
         return "betacell inverse shrinkage rate (hours)";
      }
      break;

      case 185: {
         return "Maximum adhesion force of betacells (% of full stickness)";
      }
      break;

      case 186: {
         return "Persistence of betacell polarity (time gap in min)";
      }
      break;

      case 187: {
         return "Mean betacell velocity (in microns/min)";
      }
      break;

      case 188: {
         return "Modus of betacell velocity (1=random; 2=1+weights; 3=polarity-coupled; 4=adhesion)";
      }
      break;

      case 189: {
         return "Number of betacell velocity states";
      }
      break;

      case 190: {
         return "Factor of betacell velocity reduction from v_state max to min";
      }
      break;

      case 191: {
         return "Persistence of betacell v-states in min";
      }
      break;

      case 192: {
         return "Surface tension as betacell fragment velocity in microns/min";
      }
      break;

      case 193: {
         return "betacell elongation during active move (% of radius)";
      }
      break;

      case 194: {
         return "betacell elongation for half maximum reshaping force (in radii)";
      }
      break;

      case 195: {
         return "betacell tolerance for distance to barycenter for movement to target point in %";
      }
      break;

      case 196: {
         return "Half betacell tolerance at deformation of (0.01 - infty)";
      }
      break;

      case 197: {
         return "Smoothness of betacell movement (# of timesteps for dx=a move, default 1)";
      }
      break;

      case 198: {
         return "Fix initial betacell position";
      }
      break;

      case 199: {
         return "Number of tracked betacells";
      }
      break;

      case 223: {
         return "Define a fixed glucose gradient field (0:no; 1:yes)";
      }
      break;

      case 224: {
         return "Define a constant glucose field (<0:no; >0:times resting value)";
      }
      break;

      case 225: {
         return "Define fixed glucose gradient minimum value (mM)";
      }
      break;

      case 226: {
         return "Define fixed glucose gradient minimum value (mM)";
      }
      break;

      case 328: {
         return "Use arup space (1/0) [use_arup_space]";
      }
      break;

      case 329: {
         return "Length of sequences [arup_length_sequences]";
      }
      break;

      case 330: {
         return "Nb conserved residues [arup_N_conserved]";
      }
      break;

      case 331: {
         return "Nb mutated residues [arup_N_mutates]";
      }
      break;

      case 332: {
         return "Nb shielded residues (the rest) [arup_N_shielded]";
      }
      break;

      case 333: {
         return "Number of Initial Antigen sequences (int-type) [arup_nb_ini_antigens]";
      }
      break;

      case 334: {
         return
            "Fix Antigen Sequence presentation. Order is Conserved(should be 1), Variable, shielded (should be 1). max 1000 values [arup_ini_antigens[...]]";
      }
      break;

      case 335: {
         return
            "Fraction of Arup Ags (non-fixed Ag enter with same fraction) [arup_ag_fraction[...]]";
      }
      break;

      case 336: {
         return
            "Number of mutations for mutated strains (if more initial Ag have to be generated) [arup_nb_mutations_gen_strains]";
      }
      break;

      case 337: {
         return "Activation threshold (kcal/mol) [arup_threshold_activation]";
      }
      break;

      case 338: {
         return "Initial interval for BCR sequences [arum_h_min, arup_h_max]";
      }
      break;

      case 339: {
         return "Fix initial Repertoire distribution (max 1000 values) [arup_ini_bcrs[...]]";
      }
      break;

      case 340: {
         return "Mutation rate per sequence per division per residue [arup_mutation]";
      }
      break;

      case 341: {
         return "probability of a mutation being lethal [arup_proba_lethal_mut]";
      }
      break;

      case 342: {
         return "probability of a mutation being affecting [arup_proba_affecting_mut]";
      }
      break;

      case 343: {
         return "probability silent [arup_proba_silent_mut]";
      }
      break;

      case 344: {
         return
            "distribution of affinity changes by mutation [nbLinesToRead, arup_law_mut_Xs[...], arup_law_mut_Densities[...]]";
      }
      break;

      case 345: {
         return "Coefficient of unshielding [arup_alpha]";
      }
      break;

      case 346: {
         return "Bounding of the shielding h' [arup_hprime_min, arup_hprime_max]";
      }
      break;

      case 347: {
         return "Bounding of the residues facing mutated h' [arup_hmut_min, arup_hmut_max]";
      }
      break;

      case 363: {
         return
            "Mode of BC search for Tfh help (0: death-rate; 1: fixed time; 2: pMHC-determined time)";
      }
      break;

      case 364: {
         return "Fixed time of BC search for Tfh (hours, use mode=1)";
      }
      break;

      case 365: {
         return "Added duration of BC search for Tfh per successful FDC contact (hours, use mode=2)";
      }
      break;

      case 366: {
         return "Use Tfh-signal-dependent number of divisions (0, 1-> use parameters below)";
      }
      break;

      case 367: {
         return "Tfh-signal-dependent division number Hill: half (K_P; -1 for calc with A_0)";
      }
      break;

      case 368: {
         return "Tfh-signal-dependent division number Hill: signal for P_0 divisions (A_0)";
      }
      break;

      case 369: {
         return "Make Tfh-signals to BCs dependent on BC-ICOSL expression (0=no; 1=yes):";
      }
      break;

      case 370: {
         return "Mode of ICOSL upregulation on BCs (0=fixed to 1; 1=Hill with time below):";
      }
      break;

      case 371: {
         return "Delay of ICOSL upregulation (hours, set previous par to 1):";
      }
      break;

      case 372: {
         return "BC keep memory of ICOSL-upregulation in the previous round (0=no; 1=yes):";
      }
      break;

       ///§§§ Philippe 25/03/2017
      case 373: {
         return "Increased antigen internalization for IgG cells (factor: IgG_BCRlevel):";
      }
      break;
      case 374: {
         return "Increased cell cycle length? for IgG cells (IgG_factor_cellcycle):";
      }
      break;
      case 375: {
         return "Increased cell divisions for IgG cells (IgG_factor_divisions):";
      }
      break;
      case 376: {
         return "Repartition of already switched founder cells (Founder_IgX[nIg_classes]):";
      }
      break;
      case 377: {
         return "All switching probabilities exponentially decay over time (decay_proba_switch, 1/hr):";
      }
      break;
      case 378: {
         return "Increased probability of leaving the GC for IgG cells (IgG_factor_leaving):";
      }
      break;
      case 379: {
         return "IgG switching happens only on cells above an affinity threshold (Affinity_threshold_IgG):";
      }
      break;
      case 380: {
         return "IgG switching happens only on cells that internalized higher amount of antigen (Int_Antigen_threshold_IgG):";
      }
      break;
      case 381: {
         return "IgG switching happens only on cells with a t cell help threshold (tc_help_IgG):";
      }
      break;
      case 386:{
          return "standard deviation number of initial divisions (stddev_initial_divisions)"; // 386  if != 0, then apply it
      break;
      }
      case 387:{
          return "standard deviation, dynamic number of division (stddev_DND)";               // 387  if != 0, then apply it
      break;
      }









       //// Philippe 2017-10-12
   case 391: {
       return "Use predefined tc dynamics";
       break;
   }
   case 392: {
       return "Predefined tc dunamics, series of time-tcNumbers (Nb_lines t1 N1 t2 N2 ...)";
       break;
   }

    // Philippe 2018-01-05 New version from git
   case 398: {
     return "FoxO reconstitution time from 0 to 1 [hours] (dT_FoxO):";
   }
     break;
   case 401: {
     return "mTORC1 production time from 0 to 1 [hours] (dT_mTORC1):";
   }
     break;
   case 399: {
     return "pMHC of half max FoxO deflection by B-T-interaction [Ag portions] (KFoxO):";
   }
     break;
   case 400: {
     return "Hill-coefficient of pMHC-dependent FoxO deflection [#] (nFoxO):";
   }
     break;







   }
   return "error";
}
ofstream&Werte::fPut(ofstream &s) {
   const char * u;
   if (timevalues == 1) {
      u = " (h)";
   } else {
      u = " (/h)";
   }

   int i;
   s << "=========================================\n";
   s << "General properties:\n";
   s << "=========================================\n";
   s << kenntext(1) << ":\n";
   s << system << "\n";
   s << kenntext(9) << ":\n";
   s << safety_checks << "\n";
   s << kenntext(2) << ":\n";
   s << ini_random << "\n";
   s << kenntext(5) << ":\n";
   s << late_ini_random << "\n";
   s << kenntext(6) << ":\n";
   s << timevalues << "\n";
   s << kenntext(7) << ":\n";
   s << outputfiles << "\n";
   s << kenntext(8) << ":\n";
   s << show_Ki67 << "\n";
   s << kenntext(124) << ":\n";
   s << show_mode << "\n";
   s << kenntext(172) << ":\n";
   s << CB_Narray << "\n";
   s << kenntext(173) << ":\n";
   s << CC_Narray << "\n";
   s << kenntext(174) << ":\n";
   s << TC_Narray << "\n";
   s << kenntext(175) << ":\n";
   s << OUT_Narray << "\n";
   s << kenntext(176) << ":\n";
   s << FDC_Narray << "\n";
   s << kenntext(177) << ":\n";
   s << STROMA_Narray << "\n";
   s << kenntext(178) << ":\n";
   s << BETA_Narray << "\n";

   // s << "Use cyclic border conditions (=1):\n" << CyclicBorder << "\n";
   // s << "Use Int-Populations (0-3):\n" << UseIntPart << "\n";
   s << "\n";

   s << "=========================================\n";
   s << "Shape space:\n";
   s << "=========================================\n";
   s << kenntext(11) << ":\n" << DimShapeSpace << "\n";
   s << kenntext(10) << ":\n" << metrik << "\n";
   // s << "Mutate always (1) for t>0 only (0):\n" << MutateAlways << "\n";
   s << kenntext(12) << ":\n" << SSStates << "\n";
   s << kenntext(13) << ":\n" << SSRangePerDim << "\n";
   s << kenntext(14) << ":\n" << totalA << "\n";
   s << kenntext(15) << ":\n" << APeakNumber << "\n";
   s << kenntext(16) << ":\n";
   i = 0;
   while (takeA[i] != -1 && i < 10) {
      s << takeA[i] << "\n";
      ++i;
   }
   s << "-1\n";
   s << kenntext(323) << ":\n";
   i = 0;
   while (ag_fraction[i] != -1 && i < 100) {
      s << ag_fraction[i] << "\n";
      ++i;
   }
   s << "-1\n";
   s << kenntext(18) << ":\n" << GammaGauss << "\n";
   s << kenntext(121) << ":\n" << amplitudeGauss << "\n";
   s << "\n";

   s << "=========================================\n";
   s << "Sequence Space:\n";
   s << "=========================================\n";
   s << kenntext(321) << ":\n" << use_sequence_space << "\n";
   s << kenntext(327) << ":\n" << type_affinity_function << "\n";
   s << kenntext(349) << ":\n" << use_logarithmic_seq_affinity << "\n";
   s << kenntext(315) << ":\n" << size_sequences << "\n";
   s << kenntext(348) << ":\n" << sequence_mut_per_base << "\n";
   s << kenntext(307) << ":\n" << init_antigen_sequences << "\n";
   s << kenntext(310) << ":\n";
   i = 0;
   while (initAntigenSeqs[i].compare(std::string("-1")) && i < MAXDIM) {
      s << initAntigenSeqs[i] << "\n";
      ++i;
   }

   ///§§§ Philippe 28/03/2017
   s << "-1\n";
   s << kenntext(308) << ":\n" << max_hamming_antigens << "\n";
   s << kenntext(309) << ":\n" << min_hamming_antigens << "\n";
   s << kenntext(311) << ":\n";
   i = 0;
   while (initBCRSeqs[i].compare(std::string("-1")) && i < MAXDIM) {
      s << initBCRSeqs[i] << "\n";
      ++i;
   }
   ///§§§ Philippe 28/03/2017
   s << "-1\n";
   s << kenntext(312) << ":\n" << max_hamming_BCRs << "\n";
   s << kenntext(313) << ":\n" << min_initial_affinity_BCRs << "\n";
   s << kenntext(314) << ":\n" << max_initial_affinity_BCRs << "\n";
   s << kenntext(317) << ":\n";
   i = 0;
   while (initTCRSeqs[i].compare(std::string("-1")) && i < MAXDIM) {
      s << initTCRSeqs[i] << "\n";
      ++i;
   }
   ///§§§ Philippe 28/03/2017
   s << "-1\n";
   s << kenntext(318) << ":\n" << max_hamming_TCRs << "\n";
   s << kenntext(319) << ":\n" << min_initial_affinity_TCRs << "\n";
   s << kenntext(320) << ":\n" << max_initial_affinity_TCRs << "\n";
   s << kenntext(316) << ":\n" << R_affinity << "\n";
   s << "\n";
   s << kenntext(326) << ":\n" << max_affinity_cluster << "\n";

   s << "=========================================\n";
   s << "Signals:\n";
   s << "=========================================\n";


   ///// Philippe: changed here
   s << kenntext(390) << ":\n" << ((prefixSigFiles.size() == 0) ? string("none") : prefixSigFiles) << "\n";

   s << kenntext(169) << ":\n";
   for (i = 0; i < MAXDIMSMALL; i++) {
      s << fix_signals[i] << " ";
   }
   s << "\n";
   s << kenntext(52) << ":\n" << D_differ2CC << "\n";
   s << kenntext(55) << ":\n" << bound_differ2CC << "\n";
   s << kenntext(79) << "/(h FDC):\n" << mksignal << "\n";

   s << kenntext(165) << ":\n" << D_CXCL12 << "\n";
   s << kenntext(166) << ":\n" << mkCXCL12 << "\n";
   s << kenntext(167) << ":\n" << bound_CXCL12 << "\n";
   s << kenntext(170) << ":\n" << CXCL12crit << "\n";
   s << kenntext(237) << ":\n" << CXCL12recrit << "\n";

   s << kenntext(69) << ":\n" << D_CXCL13 << "\n";
   s << kenntext(59) << " (in Mol/(h FDC)):\n" << mkCXCL13 << "\n";
   s << kenntext(58) << " (in Mol):\n" << bound_CXCL13 << "\n";
   s << kenntext(171) << ":\n" << CXCL13crit << "\n";
   s << kenntext(238) << ":\n" << CXCL13recrit << "\n";

   s << kenntext(89) << ":\n" << D_SEMA4D << "\n";
   s << kenntext(87) << ":\n" << mk_SEMA4D << "\n";
   s << kenntext(88) << " (in Mol):\n" << bound_SEMA4D << "\n";

   s << kenntext(86) << ":\n" << D_antigen << "\n";
   s << kenntext(84) << " (in Mol):\n" << bound_ag << "\n";

   s << kenntext(54) << ":\n" << signal_mode << "\n";
   s << kenntext(56) << ":\n" << objects_transparent << "\n";
   s << kenntext(123) << ":\n" << dx_signal << "\n";
   s << "\n";

   s << "=========================================\n";
   s << "BrdU staining:\n";
   s << "=========================================\n";
   s << kenntext(352) << ":\n" << t_inject_BrdU << "\n";
   s << kenntext(353) << ":\n" << deltat_inject_BrdU << "\n";
   s << kenntext(354) << ":\n" << n_inject_BrdU << "\n";
   s << kenntext(355) << ":\n" << BrdU_detection_threshold << "\n";

   s << "=========================================\n";
   s << "Nutrients:\n";
   s << "=========================================\n";
   s << kenntext(140) << " (in mol/sec):\n" << use_glucose << "\n";
   s << kenntext(141) << " (in mol/sec):\n" << use_oxygen << "\n";
   s << kenntext(142) << " (in mol/sec):\n" << use_glucose_pro << "\n";
   s << kenntext(143) << " (in mol/sec):\n" << use_oxygen_pro << "\n";
   s << kenntext(144) << " (in Mol):\n" << bound_glucose << "\n";
   s << kenntext(145) << " (in Mol):\n" << bound_oxygen << "\n";
   s << kenntext(146) << " (in micron^2/min):\n" << D_glucose_H2O << "\n";
   s << kenntext(147) << " (in micron^2/min):\n" << D_oxygen_H2O << "\n";
   s << kenntext(148) << " (in micron^2/min):\n" << D_glucose << "\n";
   s << kenntext(149) << " (in micron^2/min):\n" << D_oxygen << "\n";
   s << kenntext(150) << " (in Mol^2):\n" << critical_nutrient << "\n";
   s << kenntext(223) << " :\n" << fix_glucose_gradient << "\n";
   s << kenntext(225) << " :\n" << fix_glucose_gradient_min << " mM\n";
   s << kenntext(226) << " :\n" << fix_glucose_gradient_max << " mM\n";
   s << kenntext(224) << " :\n" << const_dynamic_glucose_field << " times glucose resting value\n";
   s << "\n";

   s << "=========================================\n";
   s << "Space properties:\n";
   s << "=========================================\n";
   s << kenntext(20) << ":\n" << DimSpace << "\n";
   s << kenntext(21) << ":\n" << dx << "\n";
   s << kenntext(22) << ":\n" << GC_radius << "\n";
   s << kenntext(23) << ":\n" << vol_shape << "\n";
   s << kenntext(200) << ":\n" << gridsize[0] << "\n";
   s << kenntext(201) << ":\n" << gridsize[1] << "\n";
   s << kenntext(202) << ":\n" << gridsize[2] << "\n";
   s << kenntext(110) << ":\n" << obstacles << "\n";
   s << kenntext(111) << ":\n" << wall_level << "\n";
   s << kenntext(112) << ":\n" << wall_width << "\n";
   s << kenntext(113) << ":\n" << slit_number << "\n";
   s << kenntext(114) << ":\n" << slit_width << "\n";
   s << kenntext(115) << ":\n" << collagen_density << "\n";
   s << kenntext(116) << ":\n" << collagen_cluster << "\n";
   s << "\n";

   s << "=========================================\n";
   s << "Time and phases:\n";
   s << "=========================================\n";
   s << kenntext(60) << ":\n" << deltat << "\n";
   s << kenntext(61) << ":\n" << tmin << "\n";
   s << kenntext(62) << ":\n" << tmax << "\n";
   s << kenntext(63) << ":\n" << ToFileStep << "\n";
   s << kenntext(296) << ":\n" << newBCinflux_stop << "\n";
   s << kenntext(64) << ":\n" << Start_Differentiation << "\n";
   s << kenntext(127) << ":\n" << Start_Mutation << "\n";
   s << kenntext(65) << ":\n" << StartOutput << "\n";
   s << "\n";

   s << "=========================================\n";
   s << "Cell tracking:\n";
   s << "=========================================\n";
   s << kenntext(151) << ":\n" << tALL << "\n";
   s << kenntext(152) << ":\n" << tCB << "\n";
   s << kenntext(153) << ":\n" << tCC << "\n";
   s << kenntext(154) << ":\n" << tOUT << "\n";
   s << kenntext(155) << ":\n" << tTC << "\n";
   s << kenntext(199) << ":\n" << tBETA << "\n";
   s << kenntext(156) << ":\n" << trackfrom << "\n";
   s << kenntext(157) << ":\n" << trackuntil << "\n";
   s << kenntext(158) << ":\n" << track_delta_t << "\n";
   s << kenntext(159) << ":\n" << v_resolution << "\n";
   s << kenntext(160) << ":\n" << delta_v << "\n";
   s << kenntext(161) << ":\n" << s_resolution << "\n";
   s << kenntext(162) << ":\n" << delta_s << "\n\n";
   s << kenntext(163) << ":\n" << alpha_resolution << "\n";
   s << kenntext(164) << ":\n" << delta_alpha << "\n";
   s << "\n";

   s << "=========================================\n";
   s << "Cell photoactivation:\n";
   s << "=========================================\n";
   s << kenntext(219) << ":\n" << photoactivation << "\n";
   s << kenntext(220) << ":\n" << photoactivation_t0 << "\n";
   s << kenntext(221) << ":\n" << photoactivation_x0 << "\n" << photoactivation_y0 << "\n"
     << photoactivation_z0
     << "\n";
   s << kenntext(222) << ":\n" << photoactivation_delta_x << "\n" << photoactivation_delta_y
     << "\n"
     << photoactivation_delta_z << "\n";
   s << kenntext(228) << ":\n" << def_DEC205 << "\n";
   s << kenntext(229) << ":\n" << def_DEC205_t0 << "\n";
   s << kenntext(230) << ":\n" << p_DEC205 << "\n";
   s << kenntext(231) << ":\n" << inject_antiDEC205OVA << "\n";
   s << kenntext(232) << ":\n" << inject_antiDEC205OVA_t0 << "\n";
   s << kenntext(239) << ":\n" << antiDEC205OVA_tend << "\n";
   s << kenntext(240) << ":\n" << TC_dec205ova_time << "\n";
   s << kenntext(241) << ":\n" << TC_factor_dec205ova << "\n";
   s << kenntext(243) << ":\n" << DEC205_p_factor << "\n";
   s << kenntext(247) << ":\n" << DEC205_induce_CBdifferentiation << "\n";
   s << kenntext(248) << ":\n" << DEC205_forces_output << "\n";
   s << kenntext(272) << ":\n" << retain_DEC205_ag << "\n";
   s << "\n";

   s << "=========================================\n";
   s << "Cells:\n";
   s << "=========================================\n";
   s << kenntext(117) << ":\n" << chemo_max << "\n";
   s << kenntext(118) << " l/mol :\n" << chemo_steep << "\n";
   s << kenntext(119) << " mol/l :\n" << chemo_half << "\n";

   /// Philippe 2017-10-22
   s << kenntext(393) << ":\n" << chemo_max_tc << "\n";
   s << kenntext(394) << " l/mol :\n" << chemo_steep_tc << "\n";
   s << kenntext(395) << " mol/l :\n" << chemo_half_tc << "\n";
   s << kenntext(396) << ":\n" << ((use_specific_tc_chemotaxis) ? 1 : 0) << "\n";
   /// Philippe 2017-10-26
   s << kenntext(397) << ":\n" << ((proba_TC_CC_interaction) ? 1 : 0) << "\n";



   s << kenntext(120) << ":\n" << adhesion_time << "\n\n";
   s << kenntext(122) << ":\n" << p_macrophage << "\n\n";
   s << kenntext(125) << ":\n" << allow_exchange << "\n\n";
   s << kenntext(128) << ":\n" << use_specific_turning_angles << "\n\n";

   s << "=========================================\n";
   s << "Centroblasts:\n";
   s << "=========================================\n";
   s << kenntext(30) << ":\n" << totalB << "\n";
   s << kenntext(297) << ":\n" << newBCinflux_rate << "\n";
   s << kenntext(298) << ":\n" << smooth_stopBCinflux << "\n";
   s << kenntext(31) << ":\n" << CB_radius << "\n";
   s << kenntext(32) << ":\n" << dx_CB << "\n";
   s << kenntext(99) << ":\n" << totalBss << "\n";
   s << kenntext(305) << ":\n" << min_seeder_dist << "\n";
   s << kenntext(306) << ":\n" << max_seeder_dist << "\n";
   s << kenntext(33) << ":\n";
   i = 0;
   while (takeB[i] != -1 && i < 10) {
      s << takeB[i] << "\n";
      ++i;
   }
   s << "-1\n";
   s << kenntext(34) << ":\n";
   i = 0;
   while (posCB[i] != -1 && i < 10) {
      s << posCB[i] << "\n";
      ++i;
   }
   s << "-1\n";
   s << kenntext(35) << ":\n" << CBreceptor_use << "\n";
   s << kenntext(36) << ":\n" << CBreceptor_activation << "\n";
   s << kenntext(37) << ":\n" << CBreceptor_dissociation << "\n";
   s << kenntext(38) << ":\n" << CBreceptor_binding << "\n";
   s << kenntext(39) << ":\n" << CBreceptor_total << "\n";
   s << kenntext(168) << ":\n" << CXCR4down << "\n";
   s << kenntext(107) << ":\n" << CB_max_adhesion << "\n";
   s << kenntext(50) << ":\n" << D_CB << "\n";
   s << kenntext(104) << ":\n" << v_CB << "\n";
   s << kenntext(356) << ":\n" << v_CB_width << "\n";
   s << kenntext(105) << ":\n" << CB_v_modi << "\n";
   s << kenntext(131) << ":\n" << CB_n_v_states << "\n";
   s << kenntext(109) << ":\n" << v_CB_switch_deltat << "\n";
   s << kenntext(106) << ":\n" << v_CB_factor << "\n";
   s << kenntext(53) << ":\n" << distance_tolerance << "\n";
   s << kenntext(57) << ":\n" << half_tolerance_deformation << "\n";
   s << kenntext(70) << u << ":\n" << proliferate << "\n";
   s << kenntext(227) << u << ":\n" << CB_fixed_times_of_divisions << "\n";
   s << kenntext(275) << u << ":\n" << fixed_time_of_divisions_mode << "\n";
   s << kenntext(293) << u << ":\n" << CB_fixed_times_of_divisions_in_expansion << "\n";
   s << kenntext(256) << u << ":\n" << CB_dt_G1 << "\n";
   s << kenntext(257) << u << ":\n" << CB_dt_S << "\n";
   s << kenntext(258) << u << ":\n" << CB_dt_G2 << "\n";
   s << kenntext(259) << u << ":\n" << CB_dt_M << "\n";
   s << kenntext(260) << u << ":\n" << CB_dt_G0 << "\n";
   s << kenntext(261) << ":\n" << CB_dtphase_width << "\n";
   s << kenntext(263) << ":\n" << transmit_CC_delay_to_CB_cycle << "\n";
   s << kenntext(267) << ":\n" << retain_ag << "\n";
   s << kenntext(268) << ":\n" << divide_ag_asymmetric << "\n";
   s << kenntext(273) << ":\n" << asymmetric_polarity_index << "\n";
   s << kenntext(276) << ":\n" << smooth_PI << "\n";
   s << kenntext(269) << ":\n" << ag_loaded_CB_diff2output << "\n";
   s << kenntext(274) << ":\n" << ag_deleted_in_fresh_CC << "\n";
   s << kenntext(270) << ":\n" << ag_loaded_CC_directly2TFH << "\n";
   s << kenntext(271) << ":\n" << ag_loaded_CB_stop_mutation << "\n";
   s << kenntext(277) << ":\n" << BC_ag_preloaded << "\n";
   s << kenntext(19) << ":\n" << CB_maxvolume4differ2CC << "\n";
   s << kenntext(71) << ":\n" << mutation << "\n";
   s << kenntext(235) << ":\n" << mutation_after_tc << "\n";
   s << kenntext(245) << ":\n" << mutation_after_dec_tc << "\n";
   s << kenntext(236) << ":\n" << mutation_affinity_exponent << "\n";
   s << kenntext(72) << u << ":\n" << tolight << "\n";
   s << kenntext(203) << ":\n" << smooth_differentiation << "\n";
   s << kenntext(250) << ":\n" << smooth_differentiation_time << "\n";
   s << kenntext(249) << ":\n" << CB2OUT_prob << "\n";
   s << kenntext(264) << ":\n" << exit2tz << "\n";
   s << kenntext(80) << u << ":\n" << grow << "\n";
   s << kenntext(100) << ":\n" << CB_D_cytosol << "\n";
   s << kenntext(108) << ":\n" << v_CB_cytosol << "\n";
   s << kenntext(101) << ":\n" << CB_elongation << "\n";
   s << kenntext(130) << ":\n" << CB_K_elongation << "\n";
   s << kenntext(102) << ":\n" << CB_smoothmove << "\n";
   s << kenntext(103) << u << ":\n" << CB_persistence << "\n";
   s << "\n";

   s << "=========================================\n";
   s << "Blasts 2:\n";
   s << "=========================================\n";
   s << kenntext(90) << ":\n" << total_blast2 << "\n";
   s << kenntext(91) << ":\n" << blast2_radius << "\n";
   s << kenntext(92) << ":\n" << dx_blast2 << "\n";
   s << kenntext(93) << ":\n" << D_blast2 << "\n";
   s << kenntext(94) << ":\n";
   i = 0;
   while (pos_blast2[i] != -1 && i < 10) {
      s << pos_blast2[i] << "\n";
      ++i;
   }
   s << "-1\n";
   s << kenntext(95) << u << ":\n" << blast2_proliferate << "\n";
   s << kenntext(96) << u << ":\n" << blast2_grow << "\n";
   s << kenntext(97) << ":\n" << blast2_distance_tolerance << "\n";
   s << kenntext(98) << ":\n" << blast2_half_tolerance_deformation << "\n";
   s << "\n";

   s << "=========================================\n";
   s << "Centrocytes:\n";
   s << "=========================================\n";
   s << kenntext(51) << ":\n" << D_CC << "\n";
   s << kenntext(77) << u << ":\n" << apoptosis << "\n";
   s << kenntext(207) << u << ":\n" << apoptosis4FDCselected << "\n";
   s << kenntext(252) << ":\n" << ignore_apoptotic_CC << "\n";
   s << kenntext(81) << u << ":\n" << shrink << "\n";
   s << kenntext(82) << " (-1=move, in h):\n" << CC_test_delay << "\n";
   s << kenntext(205) << ":\n" << CC_ICAM_delay << "\n";
   s << kenntext(73) << u << ":\n" << selection << "\n";
   s << kenntext(206) << ":\n" << FDCsignalling << "\n";
   s << kenntext(74) << ":\n" << TCell << "\n";
   s << kenntext(75) << u << ":\n" << ccdiff << "\n";
   s << kenntext(262) << u << ":\n" << ccdiff_delay << "\n";
   s << kenntext(246) << u << ":\n" << ccdiff_delay_DEC << "\n";
   s << kenntext(76) << ":\n" << output << "\n";
   s << kenntext(244) << ":\n" << output_DEC << "\n";
   s << kenntext(204) << ":\n" << smooth_dif2out << "\n";
   s << kenntext(251) << ":\n" << smooth_dif2out_time << "\n";
   s << kenntext(255) << ":\n" << final_differentiation_rate << "\n";
   s << kenntext(78) << u << ":\n" << macrophage << "\n";
   s << kenntext(132) << ":\n" << CC_persistence << "\n";
   s << kenntext(133) << ":\n" << v_CC << "\n";
   s << kenntext(357) << ":\n" << v_CC_width << "\n";
   s << kenntext(134) << ":\n" << CC_v_modi << "\n";
   s << kenntext(135) << ":\n" << CC_n_v_states << "\n";
   s << kenntext(136) << ":\n" << v_CC_factor << "\n";
   s << kenntext(137) << ":\n" << v_CC_switch_deltat << "\n";
   s << kenntext(129) << ":\n" << CXCR5down << "\n";
   s << kenntext(283) << ":\n" << CC_apoptotic_motility_mode << "\n";
   s << kenntext(284) << ":\n" << p_apo_randomwalk << "\n";
   s << kenntext(85) << ":\n" << CC_FDC_selection << "\n";
   s << kenntext(253) << ":\n" << collectFDCsignals << "\n";
   s << kenntext(254) << ":\n" << collectFDCperiod << "\n";
   s << kenntext(351) << ":\n" << prob2kill_noFDCcontactBCs << "\n";
   s << kenntext(17) << ":\n" << present_specific_ag2TC << "\n";
   s << kenntext(294) << ":\n" << reset_antigen_after_collection << "\n";
   s << kenntext(295) << ":\n" << ignore_affinity << "\n";
   s << kenntext(126) << ":\n" << use_ab_dynamics << "\n";
   s << kenntext(208) << ":\n" << initial_ab_affinity << "\n";

   s << "----------------------------------------------------------------------\n";
   s << "CC selection at TFH\n";
   s << "----------------------------------------------------------------------\n";
   s << kenntext(49) << ":\n" << TC_CC_selection << "\n";
   s << "Possible modes of BC search for Tfh help and apoptosis:\n"
     << " 0: apoptosis rate active during search for Tfh\n"
     << " 1: fixed duration of search for Tfh, apoptosis rate applies when it ends\n"
        // Philippe : added 2017-01-05 new version from git
     << " 2: pMHC-derived duration of search for Tfh, apoptosis rate applies when it ends\n"
     << " 3: FoxO-mTORC1-regulated Tfh search period, apoptosis rate applies when FoxO>1\n";
   s << kenntext(363) << ":\n" << tc_search_duration_mode << "\n";
   s << kenntext(364) << ":\n" << tc_search_duration_fixed << "\n";
   s << kenntext(365) << ":\n" << tc_search_duration_per_FDCcontact << "\n";

   // Philippe 2018-01-05 new version from git
   s << kenntext(398) << ":\n" << dT_FoxO << "\n";
   s << kenntext(399) << ":\n" << KFoxO << "\n";
   s << kenntext(400) << ":\n" << nFoxO << "\n";
   s << kenntext(401) << ":\n" << dT_mTORC1 << "\n";

   s << kenntext(360) << ":\n" << multipleTFHcontacts << "\n";
   s << kenntext(233) << ":\n" << negativeTCselection << "\n";
   s << "Modes of how to set the interaction time between CC and Tfh:\n"
     << " 0: fixed interaction time TC_time\n"
     << " 1: Gaussian variation around TC_time with width given below\n"
     << " 2: based on the number of FDC contacts N, TC_time is the duration attributed at the\n"
     << "    number of contacts equal to the K-value of pMHC-dependent division given below\n";
   s << kenntext(362) << ":\n" << mode_of_setting_TC_time << "\n";
   s << kenntext(138) << ":\n" << TC_time << "\n";
   s << kenntext(361) << ":\n" << TC_time_width << "\n";
   s << kenntext(139) << ":\n" << TC_rescue_time << "\n";
   s << kenntext(350) << ":\n" << BCstaysonTCbyTCtime << "\n";

   s << "----------------------------------------------------------------------\n";
   s << "CD40 blocking \n";
   s << "----------------------------------------------------------------------\n";

   s << kenntext(382) << ":\n" << time_tc_selection_block << "\n";
   s << kenntext(383) << ":\n" << factor_tc_selection_block << "\n";
   ///§§§
   s << kenntext(389) << ":\n" << mode_tc_selection_block << "\n";
   s << kenntext(384) << ":\n" << time_DND_block << "\n";
   s << kenntext(385) << ":\n" << factor_DND_block << "\n";
   s << kenntext(388) << ":\n" << factor_founder_div_block << "\n";

   s << "----------------------------------------------------------------------\n";
   s << "CC: pMHC dependent division of BC induced by TFH\n";
   s << "----------------------------------------------------------------------\n";
   s << kenntext(286) << ":\n" << pMHC_dependent_division << "\n";
   s << kenntext(366) << ":\n" << signal_dependent_number_of_divisions << "\n";
   s << "Either pMHC- or Tfh-signal-dependent can be used. If Tfh-signal is used,\n"
     << "parameters below are interpreted as signal in hours instead of pMHC.\n";
   s << kenntext(292) << ":\n" << pMHC_dependent_P_standard << "\n";
   s << kenntext(291) << ":\n" << pMHC_dependent_P_min << "\n";
   s << kenntext(287) << ":\n" << pMHC_dependent_P_max << "\n";
   s << kenntext(288) << ":\n" << pMHC_dependent_nHill << "\n";
   s << kenntext(289) << ":\n" << pMHC_dependent_K << "\n";
   s << kenntext(290) << ":\n" << pMHC_dependent_pMHC_of_2divisions << "\n";
   s << kenntext(367) << ":\n" << TFHsignal_dependent_K << "\n";
   s << kenntext(368) << ":\n" << TFHsignal_of_P0divisions << "\n";

   s << "----------------------------------------------------------------------\n";
   s << "CC: ICOSL upregulation and impact on Tfh help\n";
   s << "----------------------------------------------------------------------\n";
   s << kenntext(369) << "\n" << ICOSL_dependent_Tfh_signals << "\n";
   s << kenntext(370) << "\n" << ICOSL_upregulation_mode << "\n";
   s << kenntext(371) << "\n" << ICOSL_upregulation_time << "\n";
   s << kenntext(372) << "\n" << ICOSL_memory << "\n";

   s << "\n";


   ////§§§ Philippe 2017-10-12
   s << "----------------------------------------------------------------------\n";
   s << "Predefined dynamics of T cells \n";
   s << "----------------------------------------------------------------------\n";
   s << kenntext(391) << ":\n" << ((use_predefined_tc_dynamics) ? 1 : 0) << "\n";
   s << kenntext(392) << ":\n";
   int LS = (int) TCdynT.size();
   if(TCdynNb.size() != TCdynT.size()) cerr << "ERR: setparam(391), tables TCdynt and TCdynNb don't have same size" << endl;
   s << LS << endl;
   for(unsigned int i = 0; i < LS; ++i){
       s << TCdynT[i] << "\t" << TCdynNb[i] << endl;
   }
   s << "\n";



   s << "=========================================\n";
   s << "BC class switching:\n";
   s << "=========================================\n";
   s << kenntext(278) << ":\n" << do_switch_classes << "\n";
   s << "The switch matrix contains the probabilities to switch\n"
     << "from class <line-number> to class <column-number> with\n"
     << "Ig_classes={IgM,IgG,IgE,IgA}.\n";
   s << kenntext(279) << ":\n";
   for (int j = 0; j < nIg_classes; j++) {
      for (int i = 0; i < nIg_classes; i++) {
         s << switch_matrix[i + j * nIg_classes] << " ";
      }
      s << "\n";
   }
   s << kenntext(280) << ":\n" << IgE_BCRlevel << "\n";
   s << kenntext(281) << ":\n" << IgE_factor_cellcycle << "\n";
   s << kenntext(282) << ":\n" << IgE_factor_divisions << "\n";
   s << kenntext(285) << ":\n" << CC_IgE_prob_CXCR5down << "\n";
   s << "\n";

   ///§§§ Philippe 25/03/2017
   s << kenntext(373) << ":\n" << IgG_BCRlevel<< "\n";             // 373
   s << kenntext(374) << ":\n" << IgG_factor_cellcycle<< "\n";     // 374
   s << kenntext(375) << ":\n" << IgG_factor_divisions<< "\n";     // 375
   s << kenntext(376) << ":\n";                                    // 376
   for(int i = 0; i < nIg_classes; ++i){
       s << Founder_IgX[i] << "\t";
   }
   s << "\n";
   s << kenntext(377) << ":\n" << decay_proba_switch<< "\n";        // 377 for all probabilities of switching
   s << kenntext(378) << ":\n" << IgG_factor_leaving<< "\n";        // 378
   s << kenntext(379) << ":\n" << Affinity_threshold_IgG<< "\n";    // 379
   s << kenntext(380) << ":\n" << Int_Antigen_threshold_IgG<< "\n"; // 380
   s << kenntext(381) << ":\n" << tc_help_IgG << "\n";              // 381
   s << kenntext(386) << ":\n" << stddev_initial_divisions << "\n"; // 386  if != 0, then apply it
   s << kenntext(387) << ":\n" << stddev_DND << "\n";               // 387  if != 0, then apply it

   s << "=========================================\n";
   s << "Output cells:\n";
   s << "=========================================\n";
   s << kenntext(265) << ":\n" << OUT_persistence << "\n";
   s << kenntext(266) << ":\n" << v_OUT << "\n";
   s << kenntext(358) << ":\n" << v_OUT_width << "\n";
   s << "\n";

   s << "=========================================\n";
   s << "T cells:\n";
   s << "=========================================\n";
   s << kenntext(24) << ":\n" << totalTC << "\n";
   s << kenntext(25) << ":\n" << TC_radius << "\n";
   s << kenntext(26) << ":\n" << v_TC << "\n";
   s << kenntext(359) << ":\n" << v_TC_width << "\n";
   s << kenntext(27) << ":\n" << v_TC_CC << "\n";
   s << kenntext(28) << ":\n" << TC_persistence << "\n";
   s << kenntext(242) << ":\n" << north_weight << "\n";
   s << kenntext(299) << ":\n" << do_TC_division << "\n";
   s << kenntext(300) << ":\n" << TC_doubling << "\n";
   s << kenntext(301) << ":\n" << TC_meancycle << "\n";
   s << kenntext(302) << ":\n" << TC_cyclewidth << "\n";
   s << kenntext(303) << ":\n" << TC_Ndivisions << "\n";
   s << kenntext(304) << ":\n" << dx_TC << "\n";
   s << "\n";

   s << "=========================================\n";
   s << "Follicular dendritic cells:\n";
   s << "=========================================\n";
   s << kenntext(40) << ":\n" << FDCnumber << "\n";
   s << kenntext(42) << ":\n" << FDClength << "\n";
   s << kenntext(41) << ":\n";
   i = 0;
   while (posFDC[i] != -1 && i < 100) {
      s << posFDC[i] << "\n";
      ++i;
   }
   s << "-1\n";
   s << kenntext(43) << ":\n" << FDCnetwork << "\n";
   s << kenntext(44) << ":\n" << FDCtransparent << "\n";
   s << kenntext(45) << ":\n" << FDCvesicle << "\n";
   s << kenntext(29) << ":\n" << ag_per_FDC << "\n";
   s << kenntext(46) << ":\n" << ag_saturation_FDC << "\n";
   s << kenntext(324) << ":\n" << ag_distribution_mode << "\n";
   s << kenntext(325) << ":\n" << ag_detection_mode << "\n";
   s << kenntext(68) << ":\n" << ag_threshold << "\n";

   s << kenntext(4) << ":\n";
   bool stopnow = false;
   for (i = 0; i < MAXDIM; i++) {
      if (stopnow == false) {
         s << file_output[i] << "\n";
      }
      if (file_output[i] == -1) {
         stopnow = true;
      }
   }
   s << "\n";

   s << "=========================================\n";
   s << "antibodies:\n";
   s << "=========================================\n";
   s << kenntext(83) << ":\n" << D_antibody << "\n";
   s << kenntext(47) << ":\n" << mk_ab << "\n";
   s << kenntext(48) << " (in Mol):\n" << bound_ab << "\n";
   s << kenntext(66) << ":\n" << ic_k_on << "\n";
   s << kenntext(67) << ":\n" << ic_k_off << "\n";
   s << kenntext(209) << ":\n" << antibodies_resolution << "\n";
   s << kenntext(210) << ":\n" << antibodies_production << "\n";
   s << kenntext(214) << ":\n" << antibodies_degradation << "\n";
   s << kenntext(211) << ":\n" << k_ic_exp_min << "\n";
   s << kenntext(212) << ":\n" << k_ic_exp_max << "\n";
   s << kenntext(213) << ":\n" << pm_differentiation_time << "\n";
   s << kenntext(215) << ":\n" << N_GC << "\n";
   s << kenntext(216) << ":\n" << V_blood << "\n";
   s << kenntext(217) << ":\n" << inject_antibody << "\n";
   s << kenntext(234) << ":\n" << inject_antibody_time << "\n";
   s << kenntext(218) << ":\n" << injected_antibody_affinity << "\n";
   s << kenntext(322) << ":\n" << inject_antibody_ASindex << "\n";
   s << "\n";

   s << "=========================================\n";
   s << "Arup space for sequences:\n";
   s << "=========================================\n";
   s << kenntext(328) << ":\n" << use_arup_space << "\n";
   s << kenntext(329) << ":\n" << arup_length_sequences << "\n";
   s << kenntext(330) << ":\n" << arup_N_conserved << "\n";
   s << kenntext(331) << ":\n" << arup_N_mutates << "\n";
   s << kenntext(332) << ":\n" << arup_N_shielded << "\n";
   s << kenntext(333) << ":\n" << arup_nb_ini_antigens << "\n";
   s << kenntext(334) << ":\n";
   for (int i = 0; i < (int) arup_ini_antigens.size(); ++i) {
      s << arup_ini_antigens[i] << "\n";
   }
   s << "-1";
   s << "\n";
   s << kenntext(335) << ":\n";
   for (int i = 0; i < (int) arup_ag_fraction.size(); ++i) {
      s << arup_ag_fraction[i] << "\n";
   }
   s << "-1";
   s << "\n";
   s << kenntext(336) << ":\n" << arup_nb_mutations_gen_strains << "\n";
   s << kenntext(337) << ":\n" << arup_threshold_activation << "\n";
   s << kenntext(338) << ":\n" << arup_h_min << "\n" << arup_h_max << "\n";
   s << kenntext(339) << ":\n";
   for (int i = 0; i < (int) arup_ini_bcrs.size(); ++i) {
      s << arup_ini_bcrs[i] << "\n";
   }
   s << "-1";
   s << "\n";
   s << kenntext(340) << ":\n" << arup_mutation << "\n";
   s << kenntext(341) << ":\n" << arup_proba_lethal_mut << "\n";
   s << kenntext(342) << ":\n" << arup_proba_affecting_mut << "\n";
   s << kenntext(343) << ":\n" << arup_proba_silent_mut << "\n";
   s << kenntext(344) << ":\n";
   s << arup_law_mut_Xs.size() << endl;
   if (arup_law_mut_Densities.size()
       != arup_law_mut_Xs.size()) {
      cerr
         << "ERR : from Werte, arup_law_mut_Densities and arup_law_mut_Xs do not have the same size"
         << endl;
   }
   for (int i = 0; i < (int) arup_law_mut_Xs.size(); ++i) {
      s << arup_law_mut_Xs[i] << "\t" << arup_law_mut_Densities[i] << "\n";
   }
   s << kenntext(345) << ":\n" << arup_alpha << "\n";
   s << kenntext(346) << ":\n" << arup_hprime_min << "\n" << arup_hprime_max << "\n";
   s << kenntext(347) << ":\n" << arup_hmut_min << "\n" << arup_hmut_max << "\n";

   s << "=========================================\n";
   s << "=========================================\n";
   s << "=========================================\n";
   s << "beta-cells:\n";
   s << "=========================================\n";
   s << kenntext(179) << ":\n" << BETA_Nini << "\n";
   s << kenntext(180) << ":\n" << BETA_radius << "\n";
   s << kenntext(181) << ":\n" << BETA_proliferate << "\n";
   s << kenntext(182) << ":\n" << BETA_max_pro << "\n";
   s << kenntext(183) << ":\n" << BETA_grow << "\n";
   s << kenntext(184) << ":\n" << BETA_shrink << "\n";
   s << kenntext(185) << ":\n" << BETA_max_adhesion << "\n";
   s << kenntext(186) << ":\n" << BETA_persistence << "\n";
   s << kenntext(187) << ":\n" << BETA_v << "\n";
   s << kenntext(188) << ":\n" << BETA_v_modi << "\n";
   s << kenntext(189) << ":\n" << BETA_n_v_states << "\n";
   s << kenntext(190) << ":\n" << BETA_v_factor << "\n";
   s << kenntext(191) << ":\n" << BETA_v_switch_deltat << "\n";
   s << kenntext(192) << ":\n" << BETA_v_cytosol << "\n";
   s << kenntext(193) << ":\n" << BETA_elongation << "\n";
   s << kenntext(194) << ":\n" << BETA_K_elongation << "\n";
   s << kenntext(195) << ":\n" << BETA_distance_tolerance << "\n";
   s << kenntext(196) << ":\n" << BETA_half_tolerance_deformation << "\n";
   s << kenntext(197) << ":\n" << BETA_smoothmove << "\n";
   s << kenntext(198) << ":\n";
   stopnow = false;
   for (i = 0; i < MAXDIM; i++) {
      if (stopnow == false) {
         s << BETA_pos[i] << "\n";
      }
      if (BETA_pos[i] == -1) {
         stopnow = true;
      }
   }
   s << "=========================================\n";

   return s;
}
short Werte::fFind(char * parname, ifstream &s, int n) {
   s.close();
   s.open(parname);
   s.clear();
   s.setf(ios::scientific, ios::floatfield);
   short found = 0;
   const char * find = kenntext(n);
   int lang = (int)strlen(find);
   int tmplang;
   int i;
   char d;
   while (found == 0 && s.eof() == 0) {
      i = 0;
      d = 'a';
      char tmp[MAXKENNTEXTSIZE];    // will put the line in the buffer tmp
      while (int (d) != 10 && i < MAXKENNTEXTSIZE) {
         /// Philippe (MAXKENNTEXTSIZE)
         s.get(d);
         if (int (d) != 10) {
            tmp[i] = d;
            ++i;
         }
      }
      tmplang = i;
      /* Ich verstehe hier zwar nicht warum s.eof()==0 bleibt wenn die Datei zuende
       * ist. Aber mit dem Trick Laenge==200 geht der Abbruch auch! */
      if (tmplang >= lang) {
         found = 1;
         for (i = 0; i < lang; i++) {
            if (tmp[i] != find[i]) {
               found = 0;
            }
         }
      }
   }
   if ((found == 0) && show_missing_pars) {
      cout << "      !!! --> Parameter " << n << "  >>" << kenntext(n)
           << "<<  not found! Took value in setparam.C.\n";
   } else {
      //cout << ",";
   }
   return found;
}
void Werte::fGet(char * parname, bool transform2rate) {
   int i;
   ifstream s(parname);
   if (!s) {
      ///§§§ Philippe 04-2017
      s.open((string(parname) + string(".par")).c_str());
      if(!s){
		 cerr << "ERR: " << parname << "file not found" << endl;
         return;
	  }
   }

   /// philippe changed here
   if (fFind(parname, s, 390) == 1) {
      s >>  prefixSigFiles;
      if((!prefixSigFiles.compare("none")) || (!prefixSigFiles.compare("None"))) prefixSigFiles = string("");
   }



   if (fFind(parname, s, 1) == 1) {
      s >> system;
   }
   if (fFind(parname, s, 9) == 1) {
      s >> safety_checks;
   }
   if (fFind(parname, s, 2) == 1) {
      s >> ini_random;
   }
   if (fFind(parname, s, 5) == 1) {
      s >> late_ini_random;
   }
 
   if (fFind(parname, s, 6) == 1) {
      s >> timevalues;
   }
   if (fFind(parname, s, 7) == 1) {
      s >> outputfiles;
   }
   if (fFind(parname, s, 8) == 1) {
      s >> show_Ki67;
   }

   if (fFind(parname, s, 124) == 1) {
      short show_tmp;
      s >> show_tmp;
      show_mode = representation(show_tmp);
   }

   if (fFind(parname, s, 172) == 1) {
      s >> CB_Narray;
   }
   if (fFind(parname, s, 173) == 1) {
      s >> CC_Narray;
   }
   if (fFind(parname, s, 174) == 1) {
      s >> TC_Narray;
   }
   if (fFind(parname, s, 175) == 1) {
      s >> OUT_Narray;
   }
   if (fFind(parname, s, 176) == 1) {
      s >> FDC_Narray;
   }
   if (fFind(parname, s, 177) == 1) {
      s >> STROMA_Narray;
   }
   if (fFind(parname, s, 178) == 1) {
      s >> BETA_Narray;
   }
   if (fFind(parname, s, 10) == 1) {
      s >> metrik;
   }
   if (fFind(parname, s, 11) == 1) {
      s >> DimShapeSpace;
   }
   if (fFind(parname, s, 12) == 1) {
      s >> SSStates;
   }
   if (fFind(parname, s, 20) == 1) {
      s >> DimSpace;
   }
   if (fFind(parname, s, 21) == 1) {
      s >> dx;
   }
   if (fFind(parname, s, 123) == 1) {
      s >> dx_signal;
   }
   if (dx_signal < 0.) {
      dx_signal = dx;
   }
   if (fFind(parname, s, 22) == 1) {
      s >> GC_radius;
   }
   if (fFind(parname, s, 23) == 1) {
      s >> vol_shape;
   }
   if (fFind(parname, s, 200) == 1) {
      s >> gridsize[0];
   }
   if (fFind(parname, s, 201) == 1) {
      s >> gridsize[1];
   }
   if (fFind(parname, s, 202) == 1) {
      s >> gridsize[2];
   }
   if (fFind(parname, s, 110) == 1) {
      s >> obstacles;
   }
   if (fFind(parname, s, 111) == 1) {
      s >> wall_level;
   }
   if (fFind(parname, s, 112) == 1) {
      s >> wall_width;
   }
   if (fFind(parname, s, 113) == 1) {
      s >> slit_number;
   }
   if (fFind(parname, s, 114) == 1) {
      s >> slit_width;
   }
   if (fFind(parname, s, 115) == 1) {
      s >> collagen_density;
   }
   if (fFind(parname, s, 116) == 1) {
      s >> collagen_cluster;
   }

   char c;
   long int tmp = SSRangePerDim;
   if (fFind(parname, s, 13) == 1) {
      s >> tmp;
   }
   if (tmp != int (tmp)) {
      cout << "RangePerDim-Variable is not of int-type, "
           << "adjust Dimension and SSStates!";
      cin >> c;
   }
   double dtmp = pow(SSStates, (1. / DimShapeSpace));
   SSRangePerDim = long (pow(SSStates, (1. / DimShapeSpace)));
   // cout << dtmp << "\n";
   // cout << SSRangePerDim << "\n";
   if ((SSRangePerDim - dtmp > 1E-08) || (dtmp - SSRangePerDim > 1E-08)) {
      // Adjust Values:
      ++SSRangePerDim;
      long tmpb = long (pow(double (SSRangePerDim), DimShapeSpace));
      if (tmpb != SSStates) {
         SSStates = tmpb;
         cout << "Inconsistent ShapeSpace-Parameter! "
              << "Using RangePerDim=" << SSRangePerDim << " SSStates=" << SSStates << "\n";
         cin >> c;
      }
   }
   if (tmp != SSRangePerDim) {
      cout << "Inconsistent SSRangePerDim-Value, using " << SSRangePerDim << "!";
      cin >> c;
   }

   if (fFind(parname, s, 14) == 1) {
      s >> totalA;
   }
   if (fFind(parname, s, 15) == 1) {
      s >> APeakNumber;
   }
   i = 0;
   if (fFind(parname, s, 16) == 1) {
      s >> takeA[0];
      while (takeA[i] != -1 && i < 10) {
         ++i;
         s >> takeA[i];
      }
   }
   i = 0;
   if (fFind(parname, s, 323) == 1) {
      s >> ag_fraction[0];
      while (ag_fraction[i] != -1 && i < 100) {
         ++i;
         s >> ag_fraction[i];
      }
   }

   // "Sequence Space:\n";
   if (fFind(parname, s, 321) == 1) {
      s >> use_sequence_space;
   }
   if (fFind(parname, s, 327) == 1) {
      s >> type_affinity_function;
   }
   if (fFind(parname, s, 315) == 1) {
      s >> size_sequences;
   }
   if (fFind(parname, s, 349) == 1) {
      s >> use_logarithmic_seq_affinity;
   }
   if (fFind(parname, s, 348) == 1) {
      s >> sequence_mut_per_base;
   }
   if (fFind(parname, s, 307) == 1) {
      s >> init_antigen_sequences;
   }
   i = 0;
   if (fFind(parname, s, 310) == 1) {
      s >> initAntigenSeqs[0];
      while (initAntigenSeqs[i].compare(std::string("-1")) && i < MAXDIM) {
         ++i;
         s >> initAntigenSeqs[i];
      }
   }
   if (fFind(parname, s, 308) == 1) {
      s >> max_hamming_antigens;
   }
   if (fFind(parname, s, 309) == 1) {
      s >> min_hamming_antigens;
   }
   i = 0;
   if (fFind(parname, s, 311) == 1) {
      s >> initBCRSeqs[0];
      while (initBCRSeqs[i].compare(std::string("-1")) && i < MAXDIM) {
         ++i;
         s >> initBCRSeqs[i];
      }
   }
   if (fFind(parname, s, 312) == 1) {
      s >> max_hamming_BCRs;
   }
   if (fFind(parname, s, 313) == 1) {
      s >> min_initial_affinity_BCRs;
   }
   if (fFind(parname, s, 314) == 1) {
      s >> max_initial_affinity_BCRs;
   }
   i = 0;
   if (fFind(parname, s, 317) == 1) {
      s >> initTCRSeqs[0];
      while (initTCRSeqs[i].compare(std::string("-1")) && i < MAXDIM) {
         ++i;
         s >> initTCRSeqs[i];
      }
   }
   if (fFind(parname, s, 318) == 1) {
      s >> max_hamming_TCRs;
   }
   if (fFind(parname, s, 319) == 1) {
      s >> min_initial_affinity_TCRs;
   }
   if (fFind(parname, s, 320) == 1) {
      s >> max_initial_affinity_TCRs;
   }
   if (fFind(parname, s, 316) == 1) {
      s >> R_affinity;
   }
   if (fFind(parname, s, 326) == 1) {
      s >> max_affinity_cluster;
   }
   if (fFind(parname, s, 117) == 1) {
      s >> chemo_max;
   }
   if (fFind(parname, s, 118) == 1) {
      s >> chemo_steep;
   }
   if (fFind(parname, s, 119) == 1) {
      s >> chemo_half;
   }


   /// Philippe 2017-10-22
   if (fFind(parname, s, 393) == 1) {
      s >> chemo_max_tc;
   }
   if (fFind(parname, s, 394) == 1) {
      s >> chemo_steep_tc;
   }
   if (fFind(parname, s, 395) == 1) {
      s >> chemo_half_tc;
   }
   if (fFind(parname, s, 396) == 1) {
      s >> use_specific_tc_chemotaxis;
   }

   /// Philippe 26-10-2017
   if (fFind(parname, s, 397) == 1) {
      s >> proba_TC_CC_interaction;
   }




   if (fFind(parname, s, 120) == 1) {
      s >> adhesion_time;
   }
   if (fFind(parname, s, 122) == 1) {
      s >> p_macrophage;
   }
   if (fFind(parname, s, 125) == 1) {
      s >> allow_exchange;
   }
   if (fFind(parname, s, 128) == 1) {
      s >> use_specific_turning_angles;
   }
   if (fFind(parname, s, 30) == 1) {
      s >> totalB;
   }
   if (fFind(parname, s, 99) == 1) {
      s >> totalBss;
   }
   if (fFind(parname, s, 297) == 1) {
      s >> newBCinflux_rate;
   }
   if (fFind(parname, s, 298) == 1) {
      s >> smooth_stopBCinflux;
   }
   if (fFind(parname, s, 305) == 1) {
      s >> min_seeder_dist;
   }
   if (fFind(parname, s, 306) == 1) {
      s >> max_seeder_dist;
   }
   if (fFind(parname, s, 31) == 1) {
      s >> CB_radius;
   }
   if (fFind(parname, s, 32) == 1) {
      s >> dx_CB;
   }
   i = 0;
   if (fFind(parname, s, 33) == 1) {
      s >> takeB[0];
      while (takeB[i] != -1 && i < 10) {
         ++i;
         s >> takeB[i];
      }
   }
   i = 0;
   if (fFind(parname, s, 34) == 1) {
      s >> posCB[0];
      while (posCB[i] != -1 && i < 10) {
         ++i;
         s >> posCB[i];
      }
   }
   if (fFind(parname, s, 35) == 1) {
      s >> CBreceptor_use;
   }
   if (CBreceptor_use == 1) {
      cout << "\nCBreceptor_use=1 is currently not working! --> use 2\n";
   }
   if (fFind(parname, s, 36) == 1) {
      s >> CBreceptor_activation;
   }
   if (fFind(parname, s, 37) == 1) {
      s >> CBreceptor_dissociation;
   }
   if (fFind(parname, s, 38) == 1) {
      s >> CBreceptor_binding;
   }
   if (fFind(parname, s, 39) == 1) {
      s >> CBreceptor_total;
   }
   if (fFind(parname, s, 168) == 1) {
      s >> CXCR4down;
   }
   if (fFind(parname, s, 101) == 1) {
      s >> CB_elongation;
   }
   if (fFind(parname, s, 130) == 1) {
      s >> CB_K_elongation;
   }
   if (fFind(parname, s, 102) == 1) {
      s >> CB_smoothmove;
   }
   if (fFind(parname, s, 103) == 1) {
      s >> CB_persistence;
   }
   if (fFind(parname, s, 19) == 1) {
      s >> CB_maxvolume4differ2CC;
   }
   double DCB = -1;
   if (D_CB > 0) { DCB = D_CB; }
   if (fFind(parname, s, 50) == 1) {
      s >> DCB;
   }
   if (fFind(parname, s, 104) == 1) {
      s >> v_CB;
   }
   // This may be simplified, as usage of D_CB is depricated (mmh 6.9.2016)
   if (v_CB < 0.) {
      if (DCB < 0.) {
         cout << "D_CB not found! Took standard value.\n";
      } else {
         D_CB = DCB;
      }
   } else {
      if (DCB >= 0.) {
         cout << "Error: Diffusion constant and cell velocity provided. Took diffusion!\n";
         v_CB = -1.;
         D_CB = DCB;
      } else {
         D_CB = DCB;
      }
   }
   if (fFind(parname, s, 356) == 1) {
      s >> v_CB_width;
   }
   if (fFind(parname, s, 105) == 1) {
      s >> CB_v_modi;
   }
   if (fFind(parname, s, 131) == 1) {
      s >> CB_n_v_states;
   }
   if (fFind(parname, s, 109) == 1) {
      s >> v_CB_switch_deltat;
   }
   if (fFind(parname, s, 106) == 1) {
      s >> v_CB_factor;
   }
   if (fFind(parname, s, 107) == 1) {
      s >> CB_max_adhesion;
   }
   double vCBc = v_CB_cytosol, CBDc = CB_D_cytosol;
   if (fFind(parname, s, 108) == 1) {
      s >> vCBc;
   }
   if (fFind(parname, s, 100) == 1) {
      s >> CBDc;
   }
   /* MMH removed this 6.9.2016.
    * This part is only useful for very old versions.
    * Without this, still it is taken care that either v_CB_cytosol or CB_D_cytosol is set.
    * Here, the case of no D_CB and no CB_D_cytosol and now v_CB_cytosol was rescued.
    *  else {
    *  CBDc = D_CB;    // this may also be -1
    *  if ((CBDc == -1.) && (vCBc == -1.)) {
    *     vCBc = v_CB;
    *  }
    *  }
    */
   if (vCBc < 0.) {
      if (CBDc < 0.) {
         cout << "\nCB_D_cytosol not found! Took standard value.\n";
      } else {
         CB_D_cytosol = CBDc;
         v_CB_cytosol = vCBc;
      }
   } else {
      if (CBDc >= 0.) {
         cout << "Error: Diffusion constant and cell velocity provided for fragments.\n"
              << "Took velocity!\n";
         CB_D_cytosol = -1.;
         v_CB_cytosol = vCBc;
      } else {
         CB_D_cytosol = CBDc;
         v_CB_cytosol = vCBc;
      }
   }
   double DCC = -1;
   if (D_CC > 0) { DCC = D_CC; }
   // This may be simplified, as usage of D_CB is depricated (mmh 6.9.2016)
   if (fFind(parname, s, 51) == 1) {
      s >> DCC;
   }
   if (fFind(parname, s, 133) == 1) {
      s >> v_CC;
   }
   if (v_CC < 0.) {
      if (DCC < 0.) {
         cout << "D_CC not found! Took standard value.\n";
      } else {
         D_CC = DCC;
      }
   } else {
      if (DCC >= 0.) {
         cout << "Error: CC Diffusion constant and CC velocity provided. Took diffusion!\n";
         v_CC = -1.;
         D_CC = DCC;
      } else {
         D_CC = DCC;
      }
   }
   if (fFind(parname, s, 357) == 1) {
      s >> v_CC_width;
   }
   if (fFind(parname, s, 134) == 1) {
      s >> CC_v_modi;
   }
   if (fFind(parname, s, 135) == 1) {
      s >> CC_n_v_states;
   }
   if (fFind(parname, s, 136) == 1) {
      s >> v_CC_factor;
   }
   if (fFind(parname, s, 137) == 1) {
      s >> v_CC_switch_deltat;
   }
   if (fFind(parname, s, 129) == 1) {
      s >> CXCR5down;
   }
   if (fFind(parname, s, 283) == 1) {
      s >> CC_apoptotic_motility_mode;
   }
   if (fFind(parname, s, 284) == 1) {
      s >> p_apo_randomwalk;
   }
   if (fFind(parname, s, 132) == 1) {
      s >> CC_persistence;
   }
   if (fFind(parname, s, 85) == 1) {
      s >> CC_FDC_selection;
   }
   if (fFind(parname, s, 253) == 1) {
      s >> collectFDCsignals;
   }
   if ((CC_FDC_selection == 0) && (collectFDCsignals != 0)) {
      cout << "WARNING: CC_FDC_selection=" << CC_FDC_selection
           << " is in conflict to collectFDCsignals=" << collectFDCsignals << "\n"
           << "Reset collectFDCsignals to 0!\n";
      collectFDCsignals = 0;
   }
   if (fFind(parname, s, 254) == 1) {
      s >> collectFDCperiod;
   }
   if (fFind(parname, s, 351) == 1) {
      s >> prob2kill_noFDCcontactBCs;
   }
   if (fFind(parname, s, 17) == 1) {
      s >> present_specific_ag2TC;
   }
   if (fFind(parname, s, 126) == 1) {
      s >> use_ab_dynamics;
      if (((use_ab_dynamics == 1) || (use_ab_dynamics == 2)) && (APeakNumber > 1)) {
         cout << "ERROR! Inconsistent setting. Number of antigens is " << APeakNumber << "\n"
              << "       while adaptation of Ag-binding is based on average antibody affinity\n"
              << "       to the best antigen instead of to the antigen presented at this site.\n"
              << "       Multi-antigen is compatible with use_ab_dynamics in 0,3.\n"
              << "       --> Abort.\n";
         exit(1);
      }
   }
   if (fFind(parname, s, 208) == 1) {
      s >> initial_ab_affinity;
   }
   if (fFind(parname, s, 49) == 1) {
      s >> TC_CC_selection;
   }
   if (fFind(parname, s, 363) == 1) {
      s >> tc_search_duration_mode;
   }

   // Philippe 2018-01-05 New version from git
   if (tc_search_duration_mode == 3 && not(collectFDCsignals)) {
     cerr << "ERROR in setparam:\n"
      << "tc_search_duration==3 (par 363) has to be combined \n"
      << "with collectFDCsignals==TRUE (par 253).\n";
     exit(1);
   }

   if (fFind(parname, s, 364) == 1) {
      s >> tc_search_duration_fixed;
   }
   if (fFind(parname, s, 365) == 1) {
      s >> tc_search_duration_per_FDCcontact;
   }

   // Philippe 2018-01-05 New version from git
   if (fFind(parname, s, 398) == 1) {
     s >> dT_FoxO;
   }
   if (fFind(parname, s, 399) == 1) {
     s >> KFoxO;
   }
   if (fFind(parname, s, 400) == 1) {
     s >> nFoxO;
   }
   if (fFind(parname, s, 401) == 1) {
     s >> dT_mTORC1;
   }



   if (fFind(parname, s, 360) == 1) {
      s >> multipleTFHcontacts;
   }
   if (fFind(parname, s, 233) == 1) {
      s >> negativeTCselection;
   }
   if (negativeTCselection && multipleTFHcontacts) {
      cerr << "ERROR: It is not possible to run the code with negative selection by Tfh\n"
           << "       and signal integration from multiple contacts with Tfh.\n"
           << "       Set either negativeTCselection or multipleTFHcontacts false.\n";
      exit(1);
   }
   if (negativeTCselection && (tc_search_duration_mode > 0)) {
      cerr << "ERROR: It is not possible to run the code with negative selection by Tfh\n"
           << "       and the duration of search for Tfh determined otherwise.\n"
           << "       Set either negativeTCselection false or tc_search_duration_mode=0.\n";
      exit(1);
   }
   if (fFind(parname, s, 362) == 1) {
      s >> mode_of_setting_TC_time;
   }
   if ((mode_of_setting_TC_time == 2) && not (collectFDCsignals)) {
      cerr << "ERROR: It is not possible to use pMHC-dependent T-B-interaction times\n"
           << "       when collectFDCsignals is set false.\n";
      exit(1);
   }
   if (fFind(parname, s, 138) == 1) {
      s >> TC_time;
   }
   if (fFind(parname, s, 361) == 1) {
      s >> TC_time_width;
   }
   if (fFind(parname, s, 139) == 1) {
      s >> TC_rescue_time;
   }
   if (fFind(parname, s, 350) == 1) {
      s >> BCstaysonTCbyTCtime;
   }
   if (fFind(parname, s, 265) == 1) {
      s >> OUT_persistence;
   }
   // else { OUT_persistence = CC_persistence; }
   if (fFind(parname, s, 266) == 1) {
      s >> v_OUT;
   }
   // else { v_OUT = v_CC; }
   if (fFind(parname, s, 358) == 1) {
      s >> v_OUT_width;
   }
   if (fFind(parname, s, 52) == 1) {
      s >> D_differ2CC;
   }
   if (fFind(parname, s, 165) == 1) {
      s >> D_CXCL12;
   }
   if (fFind(parname, s, 69) == 1) {
      s >> D_CXCL13;
   }
   if (fFind(parname, s, 83) == 1) {
      s >> D_antibody;
   }
   if (fFind(parname, s, 86) == 1) {
      s >> D_antigen;
   }
   if (fFind(parname, s, 89) == 1) {
      s >> D_SEMA4D;
   }
   if (fFind(parname, s, 169) == 1) {
      for (i = 0; i < MAXDIMSMALL; i++) {
         s >> fix_signals[i];
      }
   }
   if (fFind(parname, s, 54) == 1) {
      s >> signal_mode;
   }
   if ((signal_mode == 0) && (CBreceptor_use > 0)) {
      cout << "\nError: CBreceptor_use=" << CBreceptor_use
           << " does not work with Quanta-diffusion mode.\n"
           << "Use CBreceptor_use=0 (no use) instead!\n";
      CBreceptor_use = 0;
   }
   /*
    * if (DimSpace==2 && signal_mode==2) {
    * cout<<"ADI works in 3D only -> use Euler method for diffusion!\n";
    * signal_mode=1;
    * }
    */

   if (fFind(parname, s, 53) == 1) {
      s >> distance_tolerance;
   }
   if (fFind(parname, s, 57) == 1) {
      s >> half_tolerance_deformation;
   }
   if (fFind(parname, s, 79) == 1) {
      s >> mksignal;
   }
   if (fFind(parname, s, 166) == 1) {
      s >> mkCXCL12;
   }
   if (fFind(parname, s, 59) == 1) {
      s >> mkCXCL13;
   }
   if (fFind(parname, s, 47) == 1) {
      s >> mk_ab;
   }
   if (fFind(parname, s, 87) == 1) {
      s >> mk_SEMA4D;
   }
   if (fFind(parname, s, 55) == 1) {
      s >> bound_differ2CC;
   }
   if (fFind(parname, s, 167) == 1) {
      s >> bound_CXCL12;
   }
   if (fFind(parname, s, 58) == 1) {
      s >> bound_CXCL13;
   }
   if (fFind(parname, s, 48) == 1) {
      s >> bound_ab;
   }
   if (fFind(parname, s, 84) == 1) {
      s >> bound_ag;
   }
   if (fFind(parname, s, 88) == 1) {
      s >> bound_SEMA4D;
   }
   if (fFind(parname, s, 170) == 1) {
      s >> CXCL12crit;
   }
   if (fFind(parname, s, 171) == 1) {
      s >> CXCL13crit;
   }
   if (fFind(parname, s, 237) == 1) {
      s >> CXCL12recrit;
   }
   if (fFind(parname, s, 238) == 1) {
      s >> CXCL13recrit;
   }
   if (fFind(parname, s, 56) == 1) {
      s >> objects_transparent;
   }

   if (fFind(parname, s, 24) == 1) {
      s >> totalTC;
   }
   if (fFind(parname, s, 25) == 1) {
      s >> TC_radius;
   }
   if (fFind(parname, s, 26) == 1) {
      s >> v_TC;
   }
   if (fFind(parname, s, 359) == 1) {
      s >> v_TC_width;
   }
   if (fFind(parname, s, 27) == 1) {
      s >> v_TC_CC;
   }
   if (fFind(parname, s, 28) == 1) {
      s >> TC_persistence;
   }
   if (fFind(parname, s, 242) == 1) {
      s >> north_weight;
   }
   if (fFind(parname, s, 299) == 1) {
      s >> do_TC_division;
   }
   if (fFind(parname, s, 300) == 1) {
      s >> TC_doubling;
   }
   if (fFind(parname, s, 301) == 1) {
      s >> TC_meancycle;
   }
   if (fFind(parname, s, 302) == 1) {
      s >> TC_cyclewidth;
   }
   if (fFind(parname, s, 303) == 1) {
      s >> TC_Ndivisions;
   }
   if (fFind(parname, s, 304) == 1) {
      s >> dx_TC;
   }

   if (fFind(parname, s, 40) == 1) {
      s >> FDCnumber;
   }
   i = 0;
   if (fFind(parname, s, 41) == 1) {
      s >> posFDC[0];
      while (posFDC[i] != -1 && i < 100) {
         ++i;
         s >> posFDC[i];
      }
   }
   if (fFind(parname, s, 42) == 1) {
      s >> FDClength;
   }
   if (fFind(parname, s, 43) == 1) {
      s >> FDCnetwork;
   }
   if (fFind(parname, s, 44) == 1) {
      s >> FDCtransparent;
   }
   if ((FDCtransparent == 0) && (objects_transparent == 0)) {
      cout << "Objects non-transparent for signal diffusion are inconsistent with\n"
           << "non-transparent FDCs -> treat FDC-dendrites as transparent for objects!\n";
      FDCtransparent = 1;
      // ### Problem zu loesen: Produktion in Zellen und Austritt der Molekuele aus diesen!!!
   }
   if (fFind(parname, s, 45) == 1) {
      s >> FDCvesicle;
   }
   if (fFind(parname, s, 29) == 1) {
      s >> ag_per_FDC;
   }
   if (fFind(parname, s, 46) == 1) {
      s >> ag_saturation_FDC;
   }
   if (fFind(parname, s, 324) == 1) {
      s >> ag_distribution_mode;
   }
   if (fFind(parname, s, 325) == 1) {
      s >> ag_detection_mode;
   }
   if (fFind(parname, s, 68) == 1) {
      s >> ag_threshold;
   }
   if (fFind(parname, s, 66) == 1) {
      s >> ic_k_on;
   }
   if (fFind(parname, s, 67) == 1) {
      s >> ic_k_off;
   }
   if (fFind(parname, s, 209) == 1) {
      s >> antibodies_resolution;
   }
   if (fFind(parname, s, 210) == 1) {
      s >> antibodies_production;
   }
   if ((antibodies_production > 0) && (mk_ab > 0)) {
      cerr << "Two modes of antibody production by plasma cells activated.\n"
           << "antibodies_production>0 induces Ab-production on the affinity space.\n"
           << "mk_ab>0 induces Ab-production of one affinity outside affinity space.\n"
           << "Make either mk_ab or antibodies_production negative!\n";
      exit(1);
   }
   if ((mk_ab > 0) && (antibodies_resolution > 0)) {
      cerr << "Antibody production and diffusion on the lattice is in conflict\n"
           << "to using antibody classes in bins.\n"
           << "Either set mk_ab<0 or set antibodies_resolution=0.\n";
      exit(1);
   }
   if ((antibodies_resolution > 0) && (ag_threshold <= 0)) {
      cerr << "Antibody feedback is set on (antibodies_resolution>0)\n"
           << "but no Ag-portion is set (ag_threshold<=0).\n"
           << "Ab feedback requires Ag dynamics. Set ag_threshold>0.\n";
      exit(1);
   }
   if (fFind(parname, s, 214) == 1) {
      s >> antibodies_degradation;
   }
   if (fFind(parname, s, 211) == 1) {
      s >> k_ic_exp_min;
   }
   if (fFind(parname, s, 212) == 1) {
      s >> k_ic_exp_max;
   }
   if (fFind(parname, s, 213) == 1) {
      s >> pm_differentiation_time;
   }
   if (fFind(parname, s, 215) == 1) {
      s >> N_GC;
   }
   if (fFind(parname, s, 216) == 1) {
      s >> V_blood;
   }
   if (fFind(parname, s, 217) == 1) {
      s >> inject_antibody;
   }
   if (fFind(parname, s, 234) == 1) {
      s >> inject_antibody_time;
   }
   if (fFind(parname, s, 218) == 1) {
      s >> injected_antibody_affinity;
   }
   if (fFind(parname, s, 322) == 1) {
      s >> inject_antibody_ASindex;
   }
   if (fFind(parname, s, 60) == 1) {
      s >> deltat;
   }
   if (fFind(parname, s, 61) == 1) {
      s >> tmin;
   }
   if (fFind(parname, s, 62) == 1) {
      s >> tmax;
   }
   if (fFind(parname, s, 63) == 1) {
      s >> ToFileStep;
   }
   if (fFind(parname, s, 296) == 1) {
      s >> newBCinflux_stop;
   }
   if (fFind(parname, s, 64) == 1) {
      s >> Start_Differentiation;
   }
   if (fFind(parname, s, 127) == 1) {
      s >> Start_Mutation;
   }
   // else { Start_Mutation = Start_Differentiation; }
   // removed this 6.9.2016 (MMH)
   if (fFind(parname, s, 65) == 1) {
      s >> StartOutput;
   }
   if (fFind(parname, s, 70) == 1) {
      s >> proliferate;
   }
   if (fFind(parname, s, 227) == 1) {
      s >> CB_fixed_times_of_divisions;
   }
   if (fFind(parname, s, 275) == 1) {
      s >> fixed_time_of_divisions_mode;
   }
   if (fFind(parname, s, 293) == 1) {
      s >> CB_fixed_times_of_divisions_in_expansion;
   }
   if (fFind(parname, s, 256) == 1) {
      s >> CB_dt_G1;
   }
   if (fFind(parname, s, 257) == 1) {
      s >> CB_dt_S;
   }
   if (fFind(parname, s, 258) == 1) {
      s >> CB_dt_G2;
   }
   if (fFind(parname, s, 259) == 1) {
      s >> CB_dt_M;
   }
   if (fFind(parname, s, 260) == 1) {
      s >> CB_dt_G0;
   }
   if ((CB_fixed_times_of_divisions > 0) && (CB_dt_G1 + CB_dt_G2 + CB_dt_S + CB_dt_M <= 0.)) {
      cout << "ERROR: From hyphasma11.05.4 and higher\n"
           << "       fixed times of CB divisions must be combined\n"
           << "       with explicit durations of the cell cycle phases.\n\n"
           << "       The combination of fixed division numbers with\n"
           << "       probabilistic rates of division events was inconsistent\n"
           << "       and, thus, deleted from hyphasma.\n\n"
           << "ABORT.\n\n";
      exit(1);
   }
   if (fFind(parname, s, 261) == 1) {
      s >> CB_dtphase_width;
   }
   if (fFind(parname, s, 263) == 1) {
      s >> transmit_CC_delay_to_CB_cycle;
   }
   if (fFind(parname, s, 352) == 1) {
      s >> t_inject_BrdU;
   }
   if ((t_inject_BrdU >= 0) && (CB_fixed_times_of_divisions <= 0)) {
      cout << "WARNING: BrdU injections require explicit distinction of cell cycle phases.\n"
           << "         BrdU injections are suppressed in this simulation. To activate:\n"
           << "         set CB_fixed_times_of_divisions>0 and provide phase durations.\n\n";
      t_inject_BrdU = -1;
   }
   if (fFind(parname, s, 353) == 1) {
      s >> deltat_inject_BrdU;
   }
   if (t_inject_BrdU < 0) { deltat_inject_BrdU = -1; }
   if (fFind(parname, s, 354) == 1) {
      s >> n_inject_BrdU;
   }
   if (fFind(parname, s, 355) == 1) {
      s >> BrdU_detection_threshold;
   }
   if (fFind(parname, s, 267) == 1) {
      s >> retain_ag;
   }
   if (fFind(parname, s, 268) == 1) {
      s >> divide_ag_asymmetric;
   }
   if (fFind(parname, s, 273) == 1) {
      s >> asymmetric_polarity_index;
   }
   if (fFind(parname, s, 276) == 1) {
      s >> smooth_PI;
   }
   if (fFind(parname, s, 269) == 1) {
      s >> ag_loaded_CB_diff2output;
   }
   if (fFind(parname, s, 274) == 1) {
      s >> ag_deleted_in_fresh_CC;
   }
   if (fFind(parname, s, 270) == 1) {
      s >> ag_loaded_CC_directly2TFH;
   }
   if (fFind(parname, s, 271) == 1) {
      s >> ag_loaded_CB_stop_mutation;
   }
   if (fFind(parname, s, 277) == 1) {
      s >> BC_ag_preloaded;
   }
   if (ag_loaded_CB_diff2output && not (retain_ag)) {
      cout << "ERROR! ag_loaded_CB_diff2outpt=true must be combined with\n"
           << "       retain_ag=true. It makes no sense not to retain antigen\n"
           << "       but to make differentiation dependent on retained antigen. EXIT.\n";
      exit(1);
   }
   if (ag_loaded_CB_diff2output && (divide_ag_asymmetric == 0.)) {
      cout << "ERROR! Symmetric division (divide_ag_asymmetric=0) is not comptatible\n"
           << "       with ag_loaded_CB_diff2output=true. When the differentiation signal\n"
           << "       is coupled to antigen retention, antigen needs to be distributed\n"
           << "       asymmetrically. EXIT.\n";
      exit(1);
   }
   if (retain_ag && not (collectFDCsignals)) {
      cout << "Ag-retention, asymmetric distribution of Ag on daughters and\n"
           << "differentiation of Ag-loaded CB to output can only be combined\n"
           << "with the mode of serial collection of Ag from FDC being on.\n"
           << "EXIT program. Please correct this in the parameter file!\n\n";
      exit(1);
   }
   if ((asymmetric_polarity_index < 1.0) && ag_loaded_CB_stop_mutation) {
      cout << "WARNING!!!\n"
           << "With incomplete asymmetric division (polarity_index=" << asymmetric_polarity_index
           << "%), ag_loaded_CB_stop_mutation==true leads to an artificial setting:\n"
           << "Every cell with a little bit of antigen retained will stop mutating.\n"
           << "Recommendation: Either use complete asymmetry of keep mutations running.\n\n";
   }
   if (fFind(parname, s, 71) == 1) {
      s >> mutation;
   }
   if (fFind(parname, s, 235) == 1) {
      s >> mutation_after_tc;
   }
   if (fFind(parname, s, 245) == 1) {
      s >> mutation_after_dec_tc;
   }
   if (fFind(parname, s, 236) == 1) {
      s >> mutation_affinity_exponent;
   }
   if (fFind(parname, s, 72) == 1) {
      s >> tolight;
   }
   if (fFind(parname, s, 203) == 1) {
      s >> smooth_differentiation;
   }
   if (fFind(parname, s, 250) == 1) {
      s >> smooth_differentiation_time;
   }
   if (fFind(parname, s, 249) == 1) {
      s >> CB2OUT_prob;
   }
   if (fFind(parname, s, 264) == 1) {
      s >> exit2tz;
   }
   if (fFind(parname, s, 73) == 1) {
      s >> selection;
   }
   if (fFind(parname, s, 206) == 1) {
      s >> FDCsignalling;
   }
   if (fFind(parname, s, 74) == 1) {
      s >> TCell;
   }
   if (fFind(parname, s, 75) == 1) {
      s >> ccdiff;
   }
   if (fFind(parname, s, 262) == 1) {
      s >> ccdiff_delay;
   }
   if (fFind(parname, s, 246) == 1) {
      s >> ccdiff_delay_DEC;
   }
   if (fFind(parname, s, 76) == 1) {
      s >> output;
   }
   if (fFind(parname, s, 244) == 1) {
      s >> output_DEC;
   }
   if (fFind(parname, s, 204) == 1) {
      s >> smooth_dif2out;
   }
   if (fFind(parname, s, 251) == 1) {
      s >> smooth_dif2out_time;
   }
   if (fFind(parname, s, 255) == 1) {
      s >> final_differentiation_rate;
   }
   if (fFind(parname, s, 77) == 1) {
      s >> apoptosis;
   }
   if (fFind(parname, s, 207) == 1) {
      s >> apoptosis4FDCselected;
   }
   if (fFind(parname, s, 252) == 1) {
      s >> ignore_apoptotic_CC;
   }
   if (fFind(parname, s, 78) == 1) {
      s >> macrophage;
   }
   if (fFind(parname, s, 80) == 1) {
      s >> grow;
   } else {
      if (show_missing_pars) {
         cout << "Korrigiere CB_radius von " << CB_radius;
         if (2 * CB_radius >= dx) {
            CB_radius = dx / 2.0001;
         }
         cout << " nach " << CB_radius << "!\n";
         /* Falls grow nicht angegeben ist im file, ist dafuer zu sorgen,
          * dass die CB alle Volumen 1 haben, denn die Daten stammen dann
          * von einer Parameterdatei, die Zellen auf mehreren Gitterpunkten
          * noch nicht kennen. */
      }
   }
   if (fFind(parname, s, 81) == 1) {
      s >> shrink;
   }
   if (fFind(parname, s, 82) == 1) {
      s >> CC_test_delay;
   }
   // else { CC_test_delay = -1.0; } // removed 6.9.2016 (MMH)
   if (fFind(parname, s, 205) == 1) {
      s >> CC_ICAM_delay;
   }
   if (fFind(parname, s, 18) == 1) {
      s >> GammaGauss;
   }
   if (fFind(parname, s, 121) == 1) {
      s >> amplitudeGauss;
   }
   // blast2 cells:
   if (fFind(parname, s, 90) == 1) {
      s >> total_blast2;
   }
   if (fFind(parname, s, 91) == 1) {
      s >> blast2_radius;
   }
   if (fFind(parname, s, 92) == 1) {
      s >> dx_blast2;
   }
   if (fFind(parname, s, 93) == 1) {
      s >> D_blast2;
   }
   i = 0;
   if (fFind(parname, s, 94) == 1) {
      s >> pos_blast2[0];
      while (pos_blast2[i] != -1 && i < 10) {
         ++i;
         s >> pos_blast2[i];
      }
   }
   if (fFind(parname, s, 95) == 1) {
      s >> blast2_proliferate;
   }
   if (fFind(parname, s, 96) == 1) {
      s >> blast2_grow;
   }
   if (fFind(parname, s, 97) == 1) {
      s >> blast2_distance_tolerance;
   }
   if (fFind(parname, s, 98) == 1) {
      s >> blast2_half_tolerance_deformation;
   }
   if (fFind(parname, s, 140) == 1) {
      s >> use_glucose;
   }
   if (fFind(parname, s, 141) == 1) {
      s >> use_oxygen;
   }
   if (fFind(parname, s, 142) == 1) {
      s >> use_glucose_pro;
   }
   if (fFind(parname, s, 143) == 1) {
      s >> use_oxygen_pro;
   }
   if (fFind(parname, s, 144) == 1) {
      s >> bound_glucose;
   }
   if (fFind(parname, s, 145) == 1) {
      s >> bound_oxygen;
   }
   if (fFind(parname, s, 146) == 1) {
      s >> D_glucose_H2O;
   }
   if (fFind(parname, s, 147) == 1) {
      s >> D_oxygen_H2O;
   }
   if (fFind(parname, s, 148) == 1) {
      s >> D_glucose;
   }
   if (fFind(parname, s, 149) == 1) {
      s >> D_oxygen;
   }
   if (fFind(parname, s, 150) == 1) {
      s >> critical_nutrient;
   }
   if (fFind(parname, s, 223) == 1) {
      s >> fix_glucose_gradient;
   }
   if (fFind(parname, s, 225) == 1) {
      s >> fix_glucose_gradient_min;
   }
   if (fFind(parname, s, 226) == 1) {
      s >> fix_glucose_gradient_max;
   }
   if (fFind(parname, s, 224) == 1) {
      s >> const_dynamic_glucose_field;
   }
   if (fix_glucose_gradient && (const_dynamic_glucose_field > 0)) {
      cout << "ERROR in parameter setting:\n"
           << "It is not possible to define a constant glucose field and a gradient!\n"
           << "Correct this setting in parameters 223 or 224 in the parameter file.\n";
      exit(1);
   }
   if (fFind(parname, s, 151) == 1) {
      s >> tALL;
   }
   if (fFind(parname, s, 152) == 1) {
      s >> tCB;
   }
   if (fFind(parname, s, 153) == 1) {
      s >> tCC;
   }
   if (fFind(parname, s, 154) == 1) {
      s >> tOUT;
   }
   if (fFind(parname, s, 155) == 1) {
      s >> tTC;
   }
   if (fFind(parname, s, 156) == 1) {
      s >> trackfrom;
   }
   if (fFind(parname, s, 157) == 1) {
      s >> trackuntil;
   }
   if (fFind(parname, s, 158) == 1) {
      s >> track_delta_t;
   }
   if (fFind(parname, s, 159) == 1) {
      s >> v_resolution;
   }
   if (fFind(parname, s, 160) == 1) {
      s >> delta_v;
   }
   if (fFind(parname, s, 161) == 1) {
      s >> s_resolution;
   }
   if (fFind(parname, s, 162) == 1) {
      s >> delta_s;
   }
   if (fFind(parname, s, 163) == 1) {
      s >> alpha_resolution;
   }
   if (fFind(parname, s, 164) == 1) {
      s >> delta_alpha;
   }

   if (fFind(parname, s, 219) == 1) {
      s >> photoactivation;
   }
   if (fFind(parname, s, 220) == 1) {
      s >> photoactivation_t0;
   }
   if (fFind(parname, s, 221) == 1) {
      s >> photoactivation_x0;
      s >> photoactivation_y0;
      s >> photoactivation_z0;
   }
   if (fFind(parname, s, 222) == 1) {
      s >> photoactivation_delta_x;
      s >> photoactivation_delta_y;
      s >> photoactivation_delta_z;
   }
   if (fFind(parname, s, 228) == 1) {
      s >> def_DEC205;
   }
   if (fFind(parname, s, 229) == 1) {
      s >> def_DEC205_t0;
   }
   if (fFind(parname, s, 230) == 1) {
      s >> p_DEC205;
   }
   if (fFind(parname, s, 231) == 1) {
      s >> inject_antiDEC205OVA;
   }
   if (fFind(parname, s, 232) == 1) {
      s >> inject_antiDEC205OVA_t0;
   }
   if (fFind(parname, s, 239) == 1) {
      s >> antiDEC205OVA_tend;
   }
   if (fFind(parname, s, 240) == 1) {
      s >> TC_dec205ova_time;
   }
   if (fFind(parname, s, 241) == 1) {
      s >> TC_factor_dec205ova;
   }
   if (fFind(parname, s, 243) == 1) {
      s >> DEC205_p_factor;
   }
   if (fFind(parname, s, 247) == 1) {
      s >> DEC205_induce_CBdifferentiation;
   }
   if (fFind(parname, s, 248) == 1) {
      s >> DEC205_forces_output;
   }
   if (fFind(parname, s, 272) == 1) {
      s >> retain_DEC205_ag;
   }
   if (not (retain_ag)) {
      if (retain_DEC205_ag) {
         cout << "\nINCONSISTENT PARAMETER!\n"
              << "retain-DEC205-acquired antigen is set off.\n"
              << "set retain_ag=true to activate.\n\n";
         retain_DEC205_ag = false;
      }
   }
   if (fFind(parname, s, 278) == 1) {
      s >> do_switch_classes;
   }
   if (do_switch_classes > 0) {
      if (fFind(parname, s, 279) == 1) {
         bool inrange = true;
         for (int i = 0; i < switch_dimension; i++) {
            s >> switch_matrix[i];
            if ((switch_matrix[i] > 1.0) || (switch_matrix[i] < 0.0)) {
               inrange = false;
            }
            // ### one might even check that the sum of each line is 1.
         }
         if (not (inrange)) {
            cout << "\nERROR!\n"
                 << "Values of switch_matrix are not all probabilities.\n"
                 << "As class switch is set ON, the programme is aborted.\n\n";
         }
      }
   }   // else the values of switch_matrix are just the standard values and not used anyway.
   if (fFind(parname, s, 280) == 1) {
      s >> IgE_BCRlevel;
   }
   if (fFind(parname, s, 281) == 1) {
      s >> IgE_factor_cellcycle;
   }
   if (fFind(parname, s, 282) == 1) {
      s >> IgE_factor_divisions;
   }
   if (fFind(parname, s, 285) == 1) {
      s >> CC_IgE_prob_CXCR5down;
   }

   ///§§§ Philippe 26/03/2017
   if (fFind(parname, s, 376) == 1) {
      for(int i = 0; i < nIg_classes; ++i){
         s >> Founder_IgX[i];
      }
   }
   if (fFind(parname, s, 373) == 1) {
      s >> IgG_BCRlevel;
   }
   if (fFind(parname, s, 374) == 1) {
      s >> IgG_factor_cellcycle;
   }
   if (fFind(parname, s, 375) == 1) {
      s >> IgG_factor_divisions;
   }
   if (fFind(parname, s, 376) == 1) {
      for(int i = 0; i < nIg_classes; ++i){
         s >> Founder_IgX[i];
      }
   }
   if (fFind(parname, s, 377) == 1) {
      s >> decay_proba_switch;
   }
   if (fFind(parname, s, 378) == 1) {
      s >> IgG_factor_leaving;
   }
   if (fFind(parname, s, 379) == 1) {
      s >> Affinity_threshold_IgG;
   }
   if (fFind(parname, s, 380) == 1) {
      s >> Int_Antigen_threshold_IgG;
   }
   if (fFind(parname, s, 381) == 1) {
      s >> tc_help_IgG;
   }

   ///§§§ Philippe 21-04-2017 Very important
   if (fFind(parname, s, 382) == 1) {
      s >> time_tc_selection_block;
   }
   if (fFind(parname, s, 383) == 1) {
      s >> factor_tc_selection_block;
   }
   if (fFind(parname, s, 389) == 1) {
      s >> mode_tc_selection_block;
   }

   if (fFind(parname, s, 384) == 1) {
      s >> time_DND_block;
   }
   if (fFind(parname, s, 385) == 1) {
      s >> factor_DND_block;
   }
   if (fFind(parname, s, 388) == 1) {
      s >> factor_founder_div_block;
   }
   if (fFind(parname, s, 386) == 1) {
      s >> stddev_initial_divisions;
   }
   if (fFind(parname, s, 387) == 1) {
      s >> stddev_DND;
   }





   if (fFind(parname, s, 179) == 1) {
      s >> BETA_Nini;
   }
   if (fFind(parname, s, 180) == 1) {
      s >> BETA_radius;
   }
   if (fFind(parname, s, 181) == 1) {
      s >> BETA_proliferate;
   }
   if (fFind(parname, s, 182) == 1) {
      s >> BETA_max_pro;
   }
   if (fFind(parname, s, 183) == 1) {
      s >> BETA_grow;
   }
   if (fFind(parname, s, 184) == 1) {
      s >> BETA_shrink;
   }
   if (fFind(parname, s, 185) == 1) {
      s >> BETA_max_adhesion;
   }
   if (fFind(parname, s, 186) == 1) {
      s >> BETA_persistence;
   }
   if (fFind(parname, s, 187) == 1) {
      s >> BETA_v;
   }
   if (fFind(parname, s, 188) == 1) {
      s >> BETA_v_modi;
   }
   if (fFind(parname, s, 189) == 1) {
      s >> BETA_n_v_states;
   }
   if (fFind(parname, s, 190) == 1) {
      s >> BETA_v_factor;
   }
   if (fFind(parname, s, 191) == 1) {
      s >> BETA_v_switch_deltat;
   }
   if (fFind(parname, s, 192) == 1) {
      s >> BETA_v_cytosol;
   }
   if (fFind(parname, s, 193) == 1) {
      s >> BETA_elongation;
   }
   if (fFind(parname, s, 194) == 1) {
      s >> BETA_K_elongation;
   }
   if (fFind(parname, s, 195) == 1) {
      s >> BETA_distance_tolerance;
   }
   if (fFind(parname, s, 196) == 1) {
      s >> BETA_half_tolerance_deformation;
   }
   if (fFind(parname, s, 197) == 1) {
      s >> BETA_smoothmove;
   }
   if (fFind(parname, s, 199) == 1) {
      s >> tBETA;
   }
   i = 0;
   if (fFind(parname, s, 198) == 1) {
      s >> BETA_pos[0];
      while (BETA_pos[i] != -1) {
         ++i;
         s >> BETA_pos[i];
      }
   }
   if (fFind(parname, s, 286) == 1) {
      s >> pMHC_dependent_division;
   }
   if (fFind(parname, s, 366) == 1) {
      s >> signal_dependent_number_of_divisions;
   }
   if ((pMHC_dependent_division || signal_dependent_number_of_divisions)
       && not (collectFDCsignals)) {
      cout << "WARNING: pMHC_ or signal_dependent_division is on.\n"
           << "This is in conflict to collectFDCsignals=" << collectFDCsignals << "\n"
           << "Reset pMHC_ and signal_dependent_division to 0!\n";
      pMHC_dependent_division = false;
      signal_dependent_number_of_divisions = false;
   }
   if (pMHC_dependent_division && signal_dependent_number_of_divisions) {
      cerr << "ERROR: The Tfh-derived number of divisions can only be based on\n"
           << "       either amount of collected pMHC or amount of collected signals.\n"
           << "Set parameters correspondingly.\n";
      exit(1);
   }
   if (fFind(parname, s, 369) == 1) {
      s >> ICOSL_dependent_Tfh_signals;
   }
   if (fFind(parname, s, 370) == 1) {
      s >> ICOSL_upregulation_mode;
   }
   if (fFind(parname, s, 371) == 1) {
      s >> ICOSL_upregulation_time;
   }
   if (fFind(parname, s, 372) == 1) {
      s >> ICOSL_memory;
   }









   ////§§§ Philippe 2017-10-12
   if (fFind(parname, s, 391) == 1) {
      s >> use_predefined_tc_dynamics;
   }
   if (fFind(parname, s, 392) == 1) {
       int LS = 0;
       s >> LS;
       if(LS > 10000) cerr << "ERR, reading " << kenntext(374) << ", nb of lines not specified properly" << endl;
       for(int i = 0; i < LS; ++i){
           double tc_t = 0;
           double tc_nb = 0;
           s >> tc_t >> tc_nb;
           TCdynT.push_back(tc_t);
           TCdynNb.push_back(tc_nb);
       }
   }








   if (fFind(parname, s, 287) == 1) {
      s >> pMHC_dependent_P_max;
   }
   if (fFind(parname, s, 288) == 1) {
      s >> pMHC_dependent_nHill;
   }
   if (fFind(parname, s, 289) == 1) {
      s >> pMHC_dependent_K;
   }
   if (fFind(parname, s, 290) == 1) {
      s >> pMHC_dependent_pMHC_of_2divisions;
   }
   if (fFind(parname, s, 367) == 1) {
      s >> TFHsignal_dependent_K;
   }
   if (fFind(parname, s, 368) == 1) {
      s >> TFHsignal_of_P0divisions;
   }
   if (fFind(parname, s, 291) == 1) {
      s >> pMHC_dependent_P_min;
   }
   if (fFind(parname, s, 292) == 1) {
      s >> pMHC_dependent_P_standard;
   }
   if (fFind(parname, s, 294) == 1) {
      s >> reset_antigen_after_collection;
   }
   if (fFind(parname, s, 295) == 1) {
      s >> ignore_affinity;
   }
   if (fFind(parname, s, 328) == 1) {
      s >> use_arup_space;
   }
   if (fFind(parname, s, 329) == 1) {
      s >> arup_length_sequences;
   }
   if (fFind(parname, s, 330) == 1) {
      s >> arup_N_conserved;
   }
   if (fFind(parname, s, 331) == 1) {
      s >> arup_N_mutates;
   }
   if (fFind(parname, s, 332) == 1) {
      s >> arup_N_shielded;
   }
   if (fFind(parname, s, 333) == 1) {
      s >> arup_nb_ini_antigens;
   }
   if (fFind(parname, s, 334) == 1) {
      string buffer;
      i = 0;
      s >> buffer;
      while ((buffer.compare(string("-1"))) && (i < 1000)) {
         arup_ini_antigens.push_back(buffer);
         ++i;
         s >> buffer;
      }
   }
   if (fFind(parname, s, 335) == 1) {
      double bufferd;
      i = 0;
      s >> bufferd;
      while ((bufferd != -1) && (i < 1000)) {
         arup_ag_fraction.push_back(bufferd);
         ++i;
         s >> bufferd;
      }
   }
   if (fFind(parname, s, 336) == 1) {
      s >> arup_nb_mutations_gen_strains;
   }
   if (fFind(parname, s, 337) == 1) {
      s >> arup_threshold_activation;
   }
   if (fFind(parname, s, 338) == 1) {
      s >> arup_h_min;
      s >> arup_h_max;
   }
   if (fFind(parname, s, 339) == 1) {
      string buffer;
      i = 0;
      s >> buffer;
      while ((buffer.compare(string("-1"))) && (i < 1000)) {
         arup_ini_bcrs.push_back(buffer);
         ++i;
         s >> buffer;
      }
   }
   if (fFind(parname, s, 340) == 1) {
      s >> arup_mutation;
   }
   if (fFind(parname, s, 341) == 1) {
      s >> arup_proba_lethal_mut;
   }
   if (fFind(parname, s, 342) == 1) {
      s >> arup_proba_affecting_mut;
   }
   if (fFind(parname, s, 343) == 1) {
      s >> arup_proba_silent_mut;
   }
   if (fFind(parname, s, 344) == 1) {
      int nbLinesToRead;
      s >> nbLinesToRead;
      arup_law_mut_Xs.resize(nbLinesToRead);
      arup_law_mut_Densities.resize(nbLinesToRead);
      for (i = 0; i < nbLinesToRead; ++i) {
         double b1, b2;
         s >> b1 >> b2;
         arup_law_mut_Xs[i] = b1;
         arup_law_mut_Densities[i] = b2;
      }
   }
   if (fFind(parname, s, 345) == 1) {
      s >> arup_alpha;
   }
   if (fFind(parname, s, 346) == 1) {
      s >> arup_hprime_min;
      s >> arup_hprime_max;
   }
   if (fFind(parname, s, 347) == 1) {
      s >> arup_hmut_min;
      s >> arup_hmut_max;
   }

   double ln2 = log(2.);
   // correction of rates by ln(2) to get the right half value times
   // note: this correction was missing up to hyphasma3.12.1 !
   /* Now all values in kinetic equations (k_off and k_on), and durations
    * (like persistence times) are given as inverse rates (without ln2).
    * Only those times that are related to growth or decay remain to be
    * weighted with ln2. However, some exceptions are mentioned below.
    * 18.6.2007
    */
   if ((timevalues == 1) && transform2rate) {
      // cell cycle time is doubling time
      proliferate = ln2 / proliferate;
      // no use of half times for the following items:
//      tolight = 1. / tolight;
       
      selection = 1. / selection;
      ccdiff = 1. / ccdiff;
      if (final_differentiation_rate > 0.) {
         final_differentiation_rate = 1. / final_differentiation_rate;
      }
      // if (ccdiff_delay_DEC>0.) ccdiff_DEC=1./ccdiff_DEC; // now it is a delay time (no rate
      // anymore)!
      CXCR4down = 1. / CXCR4down;
      CXCR5down = 1. / CXCR5down;
      // half value decay times are useful for decay processes:
      apoptosis = ln2 / apoptosis;
      if (apoptosis4FDCselected > 0) {
         apoptosis4FDCselected = ln2 / apoptosis4FDCselected;
      }
      if (p_apo_randomwalk > 0.) {
         p_apo_randomwalk = ln2 / p_apo_randomwalk;
      } else {
         p_apo_randomwalk = 0.;
      }
       //      macrophage = ln2 / macrophage; //danial: commented to put in conversion section of parameters
      p_macrophage = ln2 / p_macrophage;
      // and also for growth and shrink processes:
      grow = ln2 / grow;
      shrink = ln2 / shrink;
      // same for blast2:
      blast2_proliferate = ln2 / blast2_proliferate;
      blast2_grow = ln2 / blast2_grow;

      // TC doubling time
      TC_doubling = ln2 / TC_doubling;

      if (BETA_proliferate < 0) {
         BETA_proliferate = 0;
      } else {
         BETA_proliferate = ln2 / BETA_proliferate;
      }
      BETA_grow = ln2 / BETA_grow;
      BETA_shrink = ln2 / BETA_shrink;

      /*
       * // Unclear to me why this was corrected with ln2 in former versions!?
       * // if receptors are used:
       * if (CBreceptor_use==2) {
       * CBreceptor_binding=ln2/CBreceptor_binding;
       * CBreceptor_dissociation=ln2/CBreceptor_dissociation;
       * }
       */
      // no correction of signal-production necessary!
      // mksignal=ln2*mksignal;
   }

   int tmpi;
   if (fFind(parname, s, 4) == 1) {
      for (i = 0; i < MAXDIM; i++) {
         s >> tmpi;
         file_output[i] = tmpi;
      }
   }
   s.close();
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++ betaWerte +++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void betaWerte::ini() {
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // The present standard values are derived from bc0072.par
   // Deviation from other values that are considered to be realistic are marked with ###
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // Some flags (do not change these values to make old parameter-files work properly):
   use_Nernst = 1;
   set_leakage_zero = 0;
   use_inactivation = 1;

   use_dynamic_tau_K_V = 1;       // note that if =1 p.val.tau_K_V is used as half max tau
   use_dynamic_tau_Na_V = 1;      // note that if =1 p.val.tau_Na_V is not used
   use_dynamic_tau_fNa_V = 1;     // note that if =1 p.val.tau_Na_V is not used
   use_dynamic_H_K_Ca = 1;        // note that if =1 p.val.H_K_Ca is not used
   use_voltage_gating_K_Ca = 1;   // switches on/off voltage-gating of K,Ca channel
   Vbar_Ca_delta = 78.0;          // subtracts this value [mV] from Nernst-Vbar_Ca
   use_dynamic_IP3 = 0;           // if =0 IP3=IP3_0 all the time

   // time is calculated in seconds !!!
   dt = 0.001;         // step size dt = 1 msec (internally this is half of the value given here)
   dy = 1.e-5;         // ### maximum tolerance of double t-step deviation (in cell fractions)
   t_0 = 0.;           // inital t
   t_max = 60;         // final t  (1 hour)
   dt_output = 0.02;   // step size of writing in output file (every 0.1 seconds)

   T = 310.;         // K (body temperature)
   V_0 = -70.;       // mV (generally agreed on)
   C_m = 0.0009;     // pF/micron^2 (Michele fact sheet, Gentet et al. Biophys. J. 2000)
                     //             (0.01 in Chay 1997) (Bertram04: 0.001)
   R_bc = 6.1;       // micron (Michele fact sheet, Straub et al. Diabetes 2004)
   Sur_ER = 63.62;   // micron^2 value consisten with spherical volume
   Vol_ER = 47.71;   // micron^3 ER-volume is 5% of cell volume with R_bc=6.1 micron
   // Vbar_K=-75;           // mV (Chay 1997)
   // Vbar_Na=80;           // mV (Chay 1997)
   // Vbar_Ca=50;           // mV (Erler 2004)
   cal_0 = 0.0594;      // mM (adapted to 1% free calcium with 0.1microM resting calcium,
                        //     0.3 for Sherman88, Erler 2004 0.1)
   K_cal = 0.0005;      // mM (Erler 2004)
   buf_0 = 0.;          // switched off
   K_buf = 0.0006;      // fluorescent marker (Erler 2004)
   buf_ER_0 = 0.0594;   // mM same as calmodulin
   K_buf_ER = 0.0005;   // mM same as calmodulin
   glu_0 = 1.;          // mM (Rorsman 2003) (is that extracellular?)
                        // below 5mM no Ca osciallations in Beauvois 2006
   IP3_0 = 0.00033;     // mM (Fridlyand 2003, theo-paper)

   K_0 = 95.0;       // mM (Atwater 1978)
   Na_0 = 20.0;      // mM (Atwater 1978)
   Ca_0 = 0.0001;    // mM (Erler 2004)
   K_ext = 5.7;      // mM (Vbar=-75mV (Nernst) with K=95mM in resting state)
   Na_ext = 400.0;   // mM (Vbar=80mV (Nernst) with Na=20mM in resting state)
   Ca_ext = 1.5;     // mM (Erler 2004, Vbar=128mV with Ca=0.0001mM in resting state)
   Ca_ER_0 = 0.02;   // mM (Tengholm 2001 at 1mM glucose)

   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // Na,K
   rho[NaK] = 30.;      // #/micron^2 (unknown)
   Ihat_NaK = 3.e-05;   // pA: value of 200 ATP-molecules per sec (Maixent93), i.e. 0.00003 pA
   H_NaK = 33.333;      // [K] mM (Figure 7D, Chapman83)
   n_NaK = 2.0;         // (estimate, based on reasonable recover after perturbations)
   H2_NaK = 0.1;        // [Na] mM (Figure 7A, Chapman83: 20mM;
                        //          standard value 0.1 used for compatibility)
   n2_NaK = 2.0;        // (Chapman83)
   alpha_NaK = 1.5;     // 3Na+(out):2K+(in) (fixed)
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // K,ATP

   rho[K_ATP] = 0.00144475;   // ###
                              // /micron^2 (Ashcroft et al. 1984: 0.3)
   gbar_K_ATP = 54.;          // pS (Cook 1984, see also Dunne 2001 and Ashcroft 1984)
                              // product rho*gbar=2pS/micron^2 is different from range in
                              // betacell.tex
   tau_K_ATP = 1.0;           // ###
                              // second (estimate for time scale of metabolism is 10 seconds)
   s_h_K_ATP = 1.2;           // mM (estimate sligthly larger than glu_0)
   kappa_K_ATP = 6.0;         // mM (unknown)
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // K,V
   rho[K_V] = 0.15;     // ###
                        // /micron^2 (Chay97 and Bertram04 (product) 0.6; Kelly 1991 0.5/mum^2)
   gbar_K_V = 10.;      // pS (Dunne 2001)
   tau_K_V = 0.030;     // seconds (Bertram04: 16ms; Kelly 1991 8-37ms; Chay97 says 0.001sec)
                        // (shall be slower than other channels (Dunne)!?)
   V_h_K_V = 1.0;       // mV (Kelly 1991; but Chay 1997 App.D -18mV, Bertram04 says -16mV)
   kappa_K_V = 8.5;     // mV (Kelly 1991; but Chay 1997 App.D 14mV, Bertram04 says 5mV)
   theta_K_V = 0.400;   // seconds (Kelly 1991 lower limit)
   W_h_K_V = -25;       // mV (Kelly 1991)
   lambda_K_V = 7.3;    // mV (Kelly 1991)

   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // K,Ca
   rho[K_Ca] = 0.976251;   // /micron^2 (product rho*gbar~Bertram04: 0.002, 10fold of Chay97 App.E)
   gbar_K_Ca = 220;        // pS (Dunne)
   H_K_Ca = 0.0025;        // mM (Barrett82 Fig.8 0.0025, Chay 1997 0.001, Bertram04 0.0003)
   n_K_Ca = 2.;            // (fit of Barrett82 Fig.8 =2,
                           //  Chay/Keizer model and Chay1997 =3, Bertram04 uses 5)
   V_h_K_Ca = -40.;        // mV (fit of Barrett82 Fig.8 -40mV)
   kappa_K_Ca = 25.;       // mV (fit of Barrett82 Fig.8 +25mV)
   tau_K_Ca = 0.1;         // seconds (no known value, functional estimate)
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // sK,Ca
   rho[sK_Ca] = 3.5;        // /micron^2 (whole cell: 0.8nS Goepel 1999 J.Gen.Physiol. with
                            // R=6.1micron)
   gbar_sK_Ca = 0.5;        // pS (Goepel 1999 J.Gen.Physiol.)
   C_sK_Ca = 0.00064;       // mM (estimated from in Goepel 1999 J.Gen.Physiol. see betacell.tex)
   kappa_sK_Ca = 0.00039;   // mM (estimated from Goepel 1999 J.Gen.Physiol. see betacell.tex)
   tau_sK_Ca = 0.075;       // seconds (Hirschberg 1998)
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // Na,V
   rho[Na_V] = 0.01;      // /micron^2 (estimated to get reasonable response after K-injection)
   gbar_Na_V = 14.0;      // pS (Bezanilla87)
                          // (product rho*gbar Chay97 App.G g_Na,L (not used here))
   tau_Na_V = 0.003;      // seconds (average value of the dynamic value in Hille 1992)
   V_h_Na_V = -35.;       // mV (Plant 1988 Pflug. Arch.)
   kappa_Na_V = 8.;       // mV (Plant 1988 Pflug. Arch.)
   theta_Na_V = 0.0046;   // seconds (Hille 1992)
   W_h_Na_V = -100.0;     // mV (Plant 1988 Pflug. Arch. -120mV, Dunne 2001: -70 to -40mV)
   lambda_Na_V = 20.0;    // mV (Plant 1988 Pflug. Arch. and Hiriart 1988:
                          //     estimated from the inactivation-range -150 to -40)
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // fNa,V (non-inactivating part of voltage dependent sodium channel)
   rho[fNa_V] = 0.0;    // /micron^2 (estimated to get reasonable response after K-injection)
   gbar_fNa_V = 14.0;   // pS (Bezanilla87)
                        // (product rho*gbar Chay97 App.G g_Na,L (not used here))
   tau_fNa_V = 0.003;   // seconds (average value of the dynamic value in Hille 1992)
   V_h_fNa_V = -35.;    // mV (Plant 1988 Pflug. Arch.)
   kappa_fNa_V = 8.;    // mV (Plant 1988 Pflug. Arch.)
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // NCX
   rho[NCX] = 0.700002;   // /micron^2 (unknown)
   Ihat_NCX = -0.0005;    // pA (Juhaszova00)
   H_NCX = 0.0018;        // mM (Blaustein99)
   n_NCX = 1.0;           // (Erler 2004, no justification)
   alpha_NCX = 3.0;       // 3Na+(in):1Ca2+(out)
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // PMCA
   rho[PMCA] = 35.3025;   // /micron^2
   Ihat_PMCA = 0.00001;   // pA Juhaszova 2000
   H_PMCA = 0.0001;       // mM Elwess et al 1997
   n_PMCA = 2.0;          // Caride et al 2001 JBiolChem
   alpha_PMCA = 1.0;      // 1 ATP(turned to ADP):1Ca2+(out) (Carafoli 2001)
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // Ca,L
   rho[Ca_L] = 0.015;   // /micron^2
                        // 0.2 would be consistent with product in Bertram04 (2pS/micron^2)
   gbar_Ca_L = 27.;     // pS (Magee 1995; N-type channels Erler 2004 Tab.1 are 14pS)
   tau_Ca_L = 0.001;    // sec (short in Bertram 2004, N-type channels Erler 2004 Tab.1 1ms)
   V_h_Ca_L = 0.;       // mV (Magee95 9mV, Vinet99 -18mV, Dunne01 estimate -20mV, and Bertram04
                        // -20mV)
   kappa_Ca_L = 12.;    // mV (Magee95 6mV, Sherman88 14mV, Dunne01 and Bertram04 12mV)
   theta_Ca_L = 10.0;   // seconds (inactivation not observed)
   W_h_Ca_L = 100.;     // mV (inactivation not observed)
   lambda_Ca_L = 10.;   // mV (inactivation not observed)
   C_Ca_L = 0.004;      // mMol (Hoefer 1997)
   n_Ca_L = 1.0;        // (Hoefer 1997)
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // Ca,T
   rho[Ca_T] = 0.128408;   // /micron^2 (unknown)
   gbar_Ca_T = 10.;        // pS (Magee 1995)
   tau_Ca_T = 0.010;       // sec (Vinet 1999)
   V_h_Ca_T = -25.;        // ###
                           // mV (Magee 1995 -32mV, Dunne 2001 -55mV)
   kappa_Ca_T = 7.;        // mV (Magee 1995)
   theta_Ca_T = 0.018;     // seconds (Vinet 1999)
   W_h_Ca_T = -67.;        // mV (Magee 1995)
   lambda_Ca_T = 6.5;      // mV (Magee 1995)
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
   rho[SERCA] = 0.;         // switched off
   Ihat_SERCA = 0.000003;   // pA (Lytton 1992)
   H_SERCA = 0.0004;        // mM calcium (Lytton 1992, 0.001mM for SERCA3)
   n_SERCA = 2.;            // (Wolosker 1998)
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
   rho[IP3] = 0.;      // switched off for the moment
   gbar_IP3 = 60.0;    // pS (Bezprozvanny 1991)
   g_IP3_max = 0.81;   // # (Mak 1998)
   use_dynamic_tau_IP3 = 0;
   tau_IP3 = 0.1;          // sec (contradicting statements Mak 1997, Marchant 1998, Meyer 1990)
   C_IP3_act = 0.0002;     // mM (Mak 1998)
   n_IP3_act = 1.9;        // (Mak 1998, Chay/Keizer: 3)
   theta_IP3 = 0.3;        // sec (Marchant 1998)
   Cbar_IP3_inh = 0.055;   // mM (Mak 1998)
   n_IP3_inh = 3.9;        // (Mak 1998, Chay/Keizer: 3)
   P_IP3 = 0.000046;       // mM (data from Mak 1998, curve corrected by MMH)
   kappa_IP3 = 0.000006;   // mM (data from Mak 1998, curve corrected by MMH)
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
   k_IP3_plus = 0.0003;   // mM/sec (Fridlyand 2003, 10fold value in Shen 1995)
   k_IP3_minus = 0.04;    // /sec (Fridlyand 2003)
   C_P = 0.001;           // mM (Shen 1995, Fridlyand 2003: 0.0004mM but different model)
   n_P = 1;               // (Shen 1995, Fridlyand 2003: 2 but different model)
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // +++++++++++++++ gap-junctions +++++++++++++++++++++++++
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
   rho[gap] = 30.;   // number of gap-junction per cell-cell contact
   gbar_gap = 5.1;   // pS
   tau_gap = 0.002;   // seconds
   gap_dynamic = false;
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

   randomise_beta_proteins = false;
   randomisation_type = equal;
   randomisation_range = 0.5;
}
betaWerte::betaWerte() {
   rho = new double[N_beta_proteins];
   ini();
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Parameter File-Formate:
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const char* betaWerte::kenntext(int nummer) {
   switch (nummer) {
      // Flags
      // use_Nernst
      case 120: {
         return "Use the Nernst-equation for dynamic reversal potential [0,1]";
      }
      break;

      // set_leakage_zero
      case 121: {
         return "Set leakage currents to zero [0,1]";
      }
      break;

      // use_inactivation
      case 122: {
         return "Activate membrane protein inactivation [0,1]";
      }
      break;

      // Vbar_Ca_delta
      case 123: {
         return "Correct Ca-reversal potential by subtracting this value [mV]";
      }
      break;

      // buf_ER_0
      case 124: {
         return "ER-buffer equilibrium concentration [mMol]";
      }
      break;

      // K_buf_ER
      case 125: {
         return "ER-buffer--Ca2+ dissociation contant [mMol]";
      }
      break;

      // Ca_ER_0
      case 126: {
         return "ER-calcium (Ca2+) equilibrium concentration [mMol]";
      }
      break;

      // use_dynamic_IP3
      case 127: {
         return "IP3 follows dynamic eqs (=1) or is constant (=0)";
      }
      break;

      // dt
      case 1: {
         return "time step size [seconds]";
      }
      break;

      // dy
      case 2: {
         return "maximum tolerance of double t-step deviation [value fraction]";
      }
      break;

      // t_0
      case 3: {
         return "initial time [seconds]";
      }
      break;

      // t_max
      case 4: {
         return "final time [seconds]";
      }
      break;

      // dt_output
      case 5: {
         return "step size of writing to output file [seconds]";
      }
      break;

      //
      // V_0 equilibrium potential
      case 6: {
         return "equilibrium potential [mV]";
      }
      break;

      // R_bc cell radius
      case 7: {
         return "cell radius [micron]";
      }
      break;

      // Sur_ER total surface of ER
      case 8: {
         return "total ER surface [micron^2]";
      }
      break;

      // Vol_ER total volume of ER
      case 9: {
         return "total ER volume [micron^3]";
      }
      break;

      // C_m capacitance per surface
      case 19: {
         return "membrane capacitance per membrane surface [pF/micron^2]";
      }
      break;

      //
      // external ion concentrations:
      // K_ext
      case 10: {
         return "K+ external concentration [mMol]";
      }
      break;

      // Na_ext
      case 11: {
         return "Na+ external concentration [mMol]";
      }
      break;

      // Ca_ext
      case 12: {
         return "Ca2+ external concentration [mMol]";
      }
      break;

      // Calcium buffers
      // cal_0 calmodulin equilibrium concentration
      case 13: {
         return "calmodulin equilibrium concentration [mMol]";
      }
      break;

      // K_cal calmodulin-Ca2+ dissociation constant
      case 14: {
         return "calmodulin-Ca2+ dissociation contant [mMol]";
      }
      break;

      // buf_0 buffer equilibrium concentration
      case 15: {
         return "buffer equilibrium concentration [mMol]";
      }
      break;

      // K_buf buffer-Ca2+ dissociation constant
      case 16: {
         return "buffer-Ca2+ dissociation contant [mMol]";
      }
      break;

      //
      // glu_0 glucose equilibrium concentration
      case 17: {
         return "glucose equilibrium concentration [mMol]";
      }
      break;

      // IP3_0 glucose equilibrium concentration
      case 18: {
         return "IP3 equilibrium concentration [mMol]";
      }
      break;

      // K_0 glucose equilibrium concentration
      case 26: {
         return "potassium (K+) equilibrium concentration [mMol]";
      }
      break;

      // Na_0 glucose equilibrium concentration
      case 27: {
         return "sodium (Na+) equilibrium concentration [mMol]";
      }
      break;

      // Ca_0 glucose equilibrium concentration
      case 28: {
         return "calcium (Ca2+) equilibrium concentration [mMol]";
      }
      break;

      //
      // Temperature
      case 29: {
         return "Temperature [K]";
      }
      break;

      //
      // sodium/potassium (Na+/K+)-exchanger:
      // rho_NaK
      case 20: {
         return "(Na+/K+)-exchanger: surface density [#/micron^2]";
      }
      break;

      // Ihat_NaK (<0)
      // current corresponds to electrical charge transported by K+ ions into of the cell
      case 21: {
         return "(Na+/K+)-exchanger: single channel maximum K+ current (>0) [pA]";
      }
      break;

      // H_NaK
      case 22: {
         return "(Na+/K+)-exchanger: half open probability [K+] [mMol]";
      }
      break;

      // n_NaK
      case 23: {
         return "(Na+/K+)-exchanger: open probability Hill coefficient";
      }
      break;

      // alpha_NaK (3:2)
      case 24: {
         return "(Na+/K+)-exchanger: stoichiometry [Na:K]";
      }
      break;

      // H2_NaK
      case 35: {
         return "(Na+/K+)-exchanger: half open probability [Na+] [mMol]";
      }
      break;

      // n2_NaK
      case 36: {
         return "(Na+/K+)-exchanger: open probability Na-Hill coefficient";
      }
      break;

      //
      // ATP-sensitive K+ channel
      // rho_K_ATP
      case 30: {
         return "ATP-sensitive (K+)-channel: surface density [#/micron^2]";
      }
      break;

      // gbar_K_ATP
      case 31: {
         return "ATP-sensitive (K+)-channel: max conductance [pS]";
      }
      break;

      // tau_K_ATP
      case 32: {
         return "ATP-sensitive (K+)-channel: glucose metabolism time scale [seconds]";
      }
      break;

      // s_h_K_ATP
      case 33: {
         return "ATP-sensitive (K+)-channel: half open prob glucose concentration [mMol]";
      }
      break;

      // kappa_K_ATP
      case 34: {
         return "ATP-sensitive (K+)-channel: open prob steepness [mMol]";
      }
      break;
      // eventual add m^3 h like in hodgkin-huxley as alternative description

      //
      // Delayed rectifier (K+)-channel, V-gated
      // rho_K_V
      case 40: {
         return "V-gated (K+)-channel: surface density [#/micron^2]";
      }
      break;

      // gbar_K_V
      case 41: {
         return "V-gated (K+)-channel: max conductance [pS]";
      }
      break;

      // tau_K_V
      case 42: {
         return "V-gated (K+)-channel: delay time scale [seconds]";
      }
      break;

      // V_h_K_V
      case 43: {
         return "V-gated (K+)-channel: half open prob potential [mV]";
      }
      break;

      // kappa_K_V
      case 44: {
         return "V-gated (K+)-channel: open prob steepness [mV]";
      }
      break;

      // theta_K_V
      case 45: {
         return "V-gated (K+)-channel: delay inactivation time scale [seconds]";
      }
      break;

      // W_h_K_V
      case 46: {
         return "V-gated (K+)-channel: half inactivation prob potential [mV]";
      }
      break;

      // lambda_K_V
      case 47: {
         return "V-gated (K+)-channel: inactivation prob steepness [mV]";
      }
      break;

      // use_dynamic_tau_K_V
      case 48: {
         return "Use fit for V-dependence of K,V-activation time (c=2*tau) [0,1]";
      }
      break;

      //
      // (Ca2+)- and V-gated (K+)-channel
      // rho_K_Ca
      case 50: {
         return "(Ca2+)- and V-gated (K+)-channel: surface density [#/micron^2]";
      }
      break;

      // gbar_K_Ca
      case 51: {
         return "(Ca2+)- and V-gated (K+)-channel: max conductance [pS]";
      }
      break;

      // H_K_Ca
      case 52: {
         return "(Ca2+)- and V-gated (K+)-channel: half open probability [Ca2+] [mMol]";
      }
      break;

      // n_K_Ca (=3 Chay/Keizer-model) (=5 Bertram/Sherman-model)
      case 53: {
         return "(Ca2+)- and V-gated (K+)-channel: open probability Hill coefficient";
      }
      break;

      // V_h_K_Ca
      case 54: {
         return "(Ca2+)- and V-gated (K+)-channel: half open potential [mV]";
      }
      break;

      // kappa_K_Ca
      case 55: {
         return "(Ca2+)- and V-gated (K+)-channel: open prob steepness [mV]";
      }
      break;

      // tau_K_Ca
      case 56: {
         return "(Ca2+)- and V-gated (K+)-channel: V-open prob time scale [seconds]";
      }
      break;

      // use_dynamic_H_K_Ca
      case 57: {
         return "Use dynamics half activation concentration (ignores C_K,Ca) [0,1]";
      }
      break;

      // use_voltage_gating_K_Ca
      case 58: {
         return "Use voltage-gating of of K,Ca-channel [0,1]";
      }
      break;

      //
      // (Ca2+)-gated small conductance (K+)-channel
      // rho_sK_Ca
      case 150: {
         return "(Ca2+)-gated S(K+)-channel: surface density [#/micron^2]";
      }
      break;

      // gbar_sK_Ca
      case 151: {
         return "(Ca2+)-gated S(K+)-channel: max conductance [pS]";
      }
      break;

      // C_sK_Ca
      case 152: {
         return "(Ca2+)-gated S(K+)-channel: half open probability [Ca2+] [mMol]";
      }
      break;

      // kappa_sK_Ca
      case 153: {
         return "(Ca2+)-gated S(K+)-channel: open prob steepness [mMol]";
      }
      break;

      // tau_sK_Ca
      case 154: {
         return "(Ca2+)-gated S(K+)-channel: Ca-open prob time scale [seconds]";
      }
      break;

      //
      // sodium (Na+)-channels, V-gated
      // rho_Na_V
      case 60: {
         return "V-gated (Na+)-channels: surface density [#/micron^2]";
      }
      break;

      // gbar_Na_V
      case 61: {
         return "V-gated (Na+)-channel: max conductance [pS]";
      }
      break;

      // tau_Na_V
      case 62: {
         return "V-gated (Na+)-channel: delay time scale [seconds]";
      }
      break;

      // V_h_Na_V
      case 63: {
         return "V-gated (Na+)-channel: half open prob potential [mV]";
      }
      break;

      // kappa_Na_V
      case 64: {
         return "V-gated (Na+)-channel: open prob steepness [mV]";
      }
      break;

      // theta_Na_V
      case 65: {
         return "V-gated (Na+)-channel: delay inactivation time scale [seconds]";
      }
      break;

      // W_h_Na_V
      case 66: {
         return "V-gated (Na+)-channel: half inactivation prob potential [mV]";
      }
      break;

      // lambda_Na_V
      case 67: {
         return "V-gated (Na+)-channel: inactivation prob steepness [mV]";
      }
      break;

      // use_dynamic_tau_Na_V
      case 68: {
         return "Use fit for V-dependence of Na,V-activation time (ignores tau) [0,1]";
      }
      break;

      //
      // sodium (Na+)-channels, V-gated, non-inactivating part
      // rho_fNa_V
      case 160: {
         return "non-inactivating V-gated (Na+)-channels: surface density [#/micron^2]";
      }
      break;

      // gbar_fNa_V
      case 161: {
         return "non-inactivating V-gated (Na+)-channel: max conductance [pS]";
      }
      break;

      // tau_fNa_V
      case 162: {
         return "non-inactivating V-gated (Na+)-channel: delay time scale [seconds]";
      }
      break;

      // V_h_fNa_V
      case 163: {
         return "non-inactivating V-gated (Na+)-channel: half open prob potential [mV]";
      }
      break;

      // kappa_fNa_V
      case 164: {
         return "non-inactivating V-gated (Na+)-channel: open prob steepness [mV]";
      }
      break;

      // use_dynamic_tau_fNa_V
      case 168: {
         return
            "Use fit for V-dependence of non-inactivating Na,V-activation time (ignores tau) [0,1]";
      }
      break;

      //
      // sodium-calcium (Na+/Ca2+)-echanger
      // rho_NCX
      case 70: {
         return "(Na+/Ca2+)-exchanger: surface density [#/micron^2]";
      }
      break;

      // Ihat_NCX (<0)
      // current corresponds to electrical charge transported by Ca2+ ions out of the cell
      case 71: {
         return "(Na+/Ca2+)-exchanger: single channel maximum Ca2+ current (<0) [pA]";
      }
      break;

      // H_NCX
      case 72: {
         return "(Na+/Ca2+)-exchanger: half open probability [Ca2+] [mMol]";
      }
      break;

      // n_NCX (=1)
      case 73: {
         return "(Na+/Ca2+)-exchanger: open probability Hill coefficient";
      }
      break;

      // alpha_NCX (=3:1)
      case 74: {
         return "(Na+/Ca2+)-exchanger: stoichiometry [Na:Ca]";
      }
      break;

      //
      // PMCA (Ca2+)-ATPase
      // rho_PMCA
      case 130: {
         return "PMCA: surface density [#/micron^2]";
      }
      break;

      // Ihat_PMCA (>0)
      // current corresponds to electrical charge transported by Ca2+ ions out of the cell
      case 131: {
         return "PMCA: single channel maximum Ca2+ current (>0) [pA]";
      }
      break;

      // H_PMCA
      case 132: {
         return "PMCA: half activity [Ca2+] [mMol]";
      }
      break;

      // n_PMCA
      case 133: {
         return "PMCA: activity Hill coefficient";
      }
      break;

      // alpha_PMCA
      case 134: {
         return "PMCA: stoichiometry [Ca:ATP]";
      }
      break;

      //
      // L-type V-gated (Ca2+)-channel
      // rho_Ca_L
      case 80: {
         return "V-gated (Ca2+)-channel L-type: surface density [#/micron^2]";
      }
      break;

      // gbar_Ca_L
      case 81: {
         return "V-gated (Ca2+)-channel L-type: max conductance [pS]";
      }
      break;

      // tau_Ca_L
      case 82: {
         return "V-gated (Ca2+)-channel L-type: delay time scale [seconds]";
      }
      break;

      // V_h_Ca_L
      case 83: {
         return "V-gated (Ca2+)-channel L-type: half open prob potential [mV]";
      }
      break;

      // kappa_Ca_L
      case 84: {
         return "V-gated (Ca2+)-channel L-type: open prob steepness [mV]";
      }
      break;

      // theta_Ca_L
      case 85: {
         return "V-gated (Ca2+)-channel L-type: delay inactivation time scale [seconds]";
      }
      break;

      // W_h_Ca_L
      case 86: {
         return "V-gated (Ca2+)-channel L-type: half inactivation prob potential [mV]";
      }
      break;

      // lambda_Ca_L
      case 87: {
         return "V-gated (Ca2+)-channel L-type: inactivation prob steepness [mV]";
      }
      break;

      // C_Ca_L
      case 88: {
         return "V-gated (Ca2+)-channel L-type: Ca half-inactivation [mM]";
      }
      break;

      // n_Ca_L
      case 89: {
         return "V-gated (Ca2+)-channel L-type: Ca-inactivation Hill-coefficient";
      }
      break;

      //
      // T-type V-gated (Ca2+)-channel
      // rho_Ca_T
      case 90: {
         return "V-gated (Ca2+)-channel T-type: surface density [#/micron^2]";
      }
      break;

      // gbar_Ca_T
      case 91: {
         return "V-gated (Ca2+)-channel T-type: max conductance [pS]";
      }
      break;

      // tau_Ca_T
      case 92: {
         return "V-gated (Ca2+)-channel T-type: delay time scale [seconds]";
      }
      break;

      // V_h_Ca_T
      case 93: {
         return "V-gated (Ca2+)-channel T-type: half open prob potential [mV]";
      }
      break;

      // kappa_Ca_T
      case 94: {
         return "V-gated (Ca2+)-channel T-type: open prob steepness [mV]";
      }
      break;

      // theta_Ca_T
      case 95: {
         return "V-gated (Ca2+)-channel T-type: delay inactivation time scale [seconds]";
      }
      break;

      // W_h_Ca_T
      case 96: {
         return "V-gated (Ca2+)-channel T-type: half inactivation prob potential [mV]";
      }
      break;

      // lambda_Ca_T
      case 97: {
         return "V-gated (Ca2+)-channel T-type: inactivation prob steepness [mV]";
      }
      break;

      //
      // SERCA (Ca2+)-pump into ER
      // rho_SERCA
      case 100: {
         return "SERCA (Ca2+)-pump into ER: surface density [#/micron^2]";
      }
      break;

      // Ihat_SERCA (>0)
      // current corresponds to electrical charge transported by Ca2+ ions into of the ER
      case 101: {
         return "SERCA (Ca2+)-pump into ER: single pump max Ca2+ current (>0) [pA]";
      }
      break;

      // H_SERCA
      case 102: {
         return "SERCA (Ca2+)-pump into ER: half open probability [Ca2+] [mMol]";
      }
      break;

      // n_SERCA (=2)
      case 103: {
         return "SERCA (Ca2+)-pump into ER: open probability Hill coefficient";
      }
      break;

      //
      // IP3-gated (Ca2+)-channels in the ER
      // rho_IP3
      case 110: {
         return "IP3-gated (Ca2+)-channels: surface density [#/micron^2]";
      }
      break;

      // gbar_IP3
      case 111: {
         return "IP3-gated (Ca2+)-channels: conductivity [pS]";
      }
      break;

      // g_IP3_max
      case 112: {
         return "IP3-gated (Ca2+)-channels: Max-activation [#]";
      }
      break;

      // use_dynamic_tau_IP3
      case 113: {
         return "IP3-gated (Ca2+)-channels: activation time =0: const; =1: 1/IP3-linear";
      }
      break;

      // tau_IP3
      case 114: {
         return "IP3-gated (Ca2+)-channels: activation time scale [sec]";
      }
      break;

      // C_IP3_act
      case 115: {
         return "IP3-gated (Ca2+)-channels: half activation [Ca2+] [mMol]";
      }
      break;

      // n_IP3_act
      case 116: {
         return "IP3-gated (Ca2+)-channels: activation Hill coefficient";
      }
      break;

      // theta_IP3
      case 117: {
         return "IP3-gated (Ca2+)-channels: inactivation time scale [sec]";
      }
      break;

      // Cbar_IP3_inh
      case 118: {
         return "IP3-gated (Ca2+)-channels: steady state half inactivation [Ca2+] [mMol]";
      }
      break;

      // n_IP3_inh
      case 119: {
         return "IP3-gated (Ca2+)-channels: inactivation Hill coefficition";
      }
      break;

      // P_IP3
      case 140: {
         return "IP3-gated (Ca2+)-channels: half inactivation [IP3] [mMol]";
      }
      break;

      // kappa_IP3
      case 141: {
         return "IP3-gated (Ca2+)-channels: inactivation steepness [IP3] [mMol]";
      }
      break;

      // k_IP3_plus
      case 142: {
         return "IP3-dynamics: IP3 production rate [mMol/sec]";
      }
      break;

      // k_IP3_minus
      case 143: {
         return "IP3-dynamics: IP3 degradation rate [/sec]";
      }
      break;

      // C_P
      case 144: {
         return "IP3-dynamics: half IP3-production [Ca2+] [mMol]";
      }
      break;

      // n_P
      case 145: {
         return "IP3-dynamics: IP3-production Hill coefficient";
      }
      break;

      // gap-junction density
      case 170: {
         return "gap-junction: number per cell-cell connection [#]";
      }
      break;

      // gap-junction conductance
      case 171: {
         return "gap-junction: conductance [pS]";
      }
      break;

      // treat gap-junction conductance dynamic
      case 172: {
         return "gap-junction: dynamic gap-junction (0=no; 1=yes)";
      }
      break;

      // time constant of conductance adaption
      case 173: {
         return "gap-junction: time constant [s]";
      }
      break;

      // randomise expression
      case 180: {
         return "Randomise expression of proteins (0=no; 1=yes)";
      }
      break;

      case 181: {
         return "Type of randomisation (0=equal; 1=poisson; 2=gauss)";
      }
      break;

      case 182: {
         return "Range of randomisation (%, width, width)";
      }
      break;
   }
   return "error";
}
ofstream&betaWerte::fPut(ofstream &s) {
   /*
    * char* u;
    * if (timevalues==1) u=" (h)";
    * else u=" (/h)";
    */

   // int i;

   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
   s << kenntext(1) << ":\n";
   s << dt << "\n";
   s << kenntext(2) << ":\n";
   s << dy << "\n";
   s << kenntext(3) << ":\n";
   s << t_0 << "\n";
   s << kenntext(4) << ":\n";
   s << t_max << "\n";
   s << kenntext(5) << ":\n";
   s << dt_output << "\n";
   s << kenntext(7) << ":\n";
   s << R_bc << "\n";
   s << kenntext(8) << ":\n";
   s << Sur_ER << "\n";
   s << kenntext(9) << ":\n";
   s << Vol_ER << "\n";
   s << kenntext(19) << ":\n";
   s << C_m << "\n";
   s << kenntext(29) << ":\n";
   s << T << "\n";
   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
   s << kenntext(6) << ":\n";
   s << V_0 << "\n";
   s << kenntext(26) << ":\n";
   s << K_0 << "\n";
   s << kenntext(27) << ":\n";
   s << Na_0 << "\n";
   s << kenntext(28) << ":\n";
   s << Ca_0 << "\n";
   s << kenntext(126) << ":\n";
   s << Ca_ER_0 << "\n";
   s << kenntext(17) << ":\n";
   s << glu_0 << "\n";
   s << kenntext(120) << ":\n";
   s << use_Nernst << "\n";
   s << kenntext(10) << ":\n";
   s << K_ext << "\n";
   s << kenntext(11) << ":\n";
   s << Na_ext << "\n";
   s << kenntext(12) << ":\n";
   s << Ca_ext << "\n";
   s << kenntext(123) << ":\n";
   s << Vbar_Ca_delta << "\n";
   s << kenntext(121) << ":\n";
   s << set_leakage_zero << "\n";
   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
   s << kenntext(13) << ":\n";
   s << cal_0 << "\n";
   s << kenntext(14) << ":\n";
   s << K_cal << "\n";
   s << kenntext(15) << ":\n";
   s << buf_0 << "\n";
   s << kenntext(16) << ":\n";
   s << K_buf << "\n";
   s << kenntext(124) << ":\n";
   s << buf_ER_0 << "\n";
   s << kenntext(125) << ":\n";
   s << K_buf_ER << "\n";
   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
   s << kenntext(18) << ":\n";
   s << IP3_0 << "\n";
   s << kenntext(127) << ":\n";
   s << use_dynamic_IP3 << "\n";
   s << kenntext(142) << ":\n";
   s << k_IP3_plus << "\n";
   s << kenntext(143) << ":\n";
   s << k_IP3_minus << "\n";
   s << kenntext(144) << ":\n";
   s << C_P << "\n";
   s << kenntext(145) << ":\n";
   s << n_P << "\n";
   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
   s << kenntext(122) << ":\n";
   s << use_inactivation << "\n";
   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
   s << kenntext(20) << ":\n";
   s << rho[NaK] << "\n";
   s << kenntext(21) << ":\n";
   s << Ihat_NaK << "\n";
   s << kenntext(22) << ":\n";
   s << H_NaK << "\n";
   s << kenntext(23) << ":\n";
   s << n_NaK << "\n";
   s << kenntext(35) << ":\n";
   s << H2_NaK << "\n";
   s << kenntext(36) << ":\n";
   s << n2_NaK << "\n";
   s << kenntext(24) << ":\n";
   s << alpha_NaK << "\n";
   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
   s << kenntext(30) << ":\n";
   s << rho[K_ATP] << "\n";
   s << kenntext(31) << ":\n";
   s << gbar_K_ATP << "\n";
   s << kenntext(32) << ":\n";
   s << tau_K_ATP << "\n";
   s << kenntext(33) << ":\n";
   s << s_h_K_ATP << "\n";
   s << kenntext(34) << ":\n";
   s << kappa_K_ATP << "\n";
   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
   s << kenntext(40) << ":\n";
   s << rho[K_V] << "\n";
   s << kenntext(41) << ":\n";
   s << gbar_K_V << "\n";
   s << kenntext(42) << ":\n";
   s << tau_K_V << "\n";
   s << kenntext(48) << ":\n";
   s << use_dynamic_tau_K_V << "\n";
   s << kenntext(43) << ":\n";
   s << V_h_K_V << "\n";
   s << kenntext(44) << ":\n";
   s << kappa_K_V << "\n";
   s << kenntext(45) << ":\n";
   s << theta_K_V << "\n";
   s << kenntext(46) << ":\n";
   s << W_h_K_V << "\n";
   s << kenntext(47) << ":\n";
   s << lambda_K_V << "\n";
   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
   s << kenntext(50) << ":\n";
   s << rho[K_Ca] << "\n";
   s << kenntext(51) << ":\n";
   s << gbar_K_Ca << "\n";
   s << kenntext(52) << ":\n";
   s << H_K_Ca << "\n";
   s << kenntext(57) << ":\n";
   s << use_dynamic_H_K_Ca << "\n";
   s << kenntext(53) << ":\n";
   s << n_K_Ca << "\n";
   s << kenntext(58) << ":\n";
   s << use_voltage_gating_K_Ca << "\n";
   s << kenntext(54) << ":\n";
   s << V_h_K_Ca << "\n";
   s << kenntext(55) << ":\n";
   s << kappa_K_Ca << "\n";
   s << kenntext(56) << ":\n";
   s << tau_K_Ca << "\n";
   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
   s << kenntext(150) << ":\n";
   s << rho[sK_Ca] << "\n";
   s << kenntext(151) << ":\n";
   s << gbar_sK_Ca << "\n";
   s << kenntext(152) << ":\n";
   s << C_sK_Ca << "\n";
   s << kenntext(153) << ":\n";
   s << kappa_sK_Ca << "\n";
   s << kenntext(154) << ":\n";
   s << tau_sK_Ca << "\n";
   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
   s << kenntext(60) << ":\n";
   s << rho[Na_V] << "\n";
   s << kenntext(61) << ":\n";
   s << gbar_Na_V << "\n";
   s << kenntext(62) << ":\n";
   s << tau_Na_V << "\n";
   s << kenntext(68) << ":\n";
   s << use_dynamic_tau_Na_V << "\n";
   s << kenntext(63) << ":\n";
   s << V_h_Na_V << "\n";
   s << kenntext(64) << ":\n";
   s << kappa_Na_V << "\n";
   s << kenntext(65) << ":\n";
   s << theta_Na_V << "\n";
   s << kenntext(66) << ":\n";
   s << W_h_Na_V << "\n";
   s << kenntext(67) << ":\n";
   s << lambda_Na_V << "\n";
   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
   s << kenntext(160) << ":\n";
   s << rho[fNa_V] << "\n";
   s << kenntext(161) << ":\n";
   s << gbar_fNa_V << "\n";
   s << kenntext(162) << ":\n";
   s << tau_fNa_V << "\n";
   s << kenntext(168) << ":\n";
   s << use_dynamic_tau_fNa_V << "\n";
   s << kenntext(163) << ":\n";
   s << V_h_fNa_V << "\n";
   s << kenntext(164) << ":\n";
   s << kappa_fNa_V << "\n";
   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
   s << kenntext(70) << ":\n";
   s << rho[NCX] << "\n";
   s << kenntext(71) << ":\n";
   s << Ihat_NCX << "\n";
   s << kenntext(72) << ":\n";
   s << H_NCX << "\n";
   s << kenntext(73) << ":\n";
   s << n_NCX << "\n";
   s << kenntext(74) << ":\n";
   s << alpha_NCX << "\n";
   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
   s << kenntext(130) << ":\n";
   s << rho[PMCA] << "\n";
   s << kenntext(131) << ":\n";
   s << Ihat_PMCA << "\n";
   s << kenntext(132) << ":\n";
   s << H_PMCA << "\n";
   s << kenntext(133) << ":\n";
   s << n_PMCA << "\n";
   s << kenntext(134) << ":\n";
   s << alpha_PMCA << "\n";
   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
   s << kenntext(80) << ":\n";
   s << rho[Ca_L] << "\n";
   s << kenntext(81) << ":\n";
   s << gbar_Ca_L << "\n";
   s << kenntext(82) << ":\n";
   s << tau_Ca_L << "\n";
   s << kenntext(83) << ":\n";
   s << V_h_Ca_L << "\n";
   s << kenntext(84) << ":\n";
   s << kappa_Ca_L << "\n";
   s << kenntext(85) << ":\n";
   s << theta_Ca_L << "\n";
   s << kenntext(86) << ":\n";
   s << W_h_Ca_L << "\n";
   s << kenntext(87) << ":\n";
   s << lambda_Ca_L << "\n";
   s << kenntext(88) << ":\n";
   s << C_Ca_L << "\n";
   s << kenntext(89) << ":\n";
   s << n_Ca_L << "\n";
   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
   s << kenntext(90) << ":\n";
   s << rho[Ca_T] << "\n";
   s << kenntext(91) << ":\n";
   s << gbar_Ca_T << "\n";
   s << kenntext(92) << ":\n";
   s << tau_Ca_T << "\n";
   s << kenntext(93) << ":\n";
   s << V_h_Ca_T << "\n";
   s << kenntext(94) << ":\n";
   s << kappa_Ca_T << "\n";
   s << kenntext(95) << ":\n";
   s << theta_Ca_T << "\n";
   s << kenntext(96) << ":\n";
   s << W_h_Ca_T << "\n";
   s << kenntext(97) << ":\n";
   s << lambda_Ca_T << "\n";
   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
   s << kenntext(100) << ":\n";
   s << rho[SERCA] << "\n";
   s << kenntext(101) << ":\n";
   s << Ihat_SERCA << "\n";
   s << kenntext(102) << ":\n";
   s << H_SERCA << "\n";
   s << kenntext(103) << ":\n";
   s << n_SERCA << "\n";
   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
   s << kenntext(110) << ":\n";
   s << rho[IP3] << "\n";
   s << kenntext(111) << ":\n";
   s << gbar_IP3 << "\n";
   s << kenntext(112) << ":\n";
   s << g_IP3_max << "\n";
   s << kenntext(113) << ":\n";
   s << use_dynamic_tau_IP3 << "\n";
   s << kenntext(114) << ":\n";
   s << tau_IP3 << "\n";
   s << kenntext(115) << ":\n";
   s << C_IP3_act << "\n";
   s << kenntext(116) << ":\n";
   s << n_IP3_act << "\n";
   s << kenntext(117) << ":\n";
   s << theta_IP3 << "\n";
   s << kenntext(118) << ":\n";
   s << Cbar_IP3_inh << "\n";
   s << kenntext(119) << ":\n";
   s << n_IP3_inh << "\n";
   s << kenntext(140) << ":\n";
   s << P_IP3 << "\n";
   s << kenntext(141) << ":\n";
   s << kappa_IP3 << "\n";
   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
   s << kenntext(170) << ":\n";
   s << rho[gap] << "\n";
   s << kenntext(171) << ":\n";
   s << gbar_gap << "\n";
   s << kenntext(172) << ":\n";
   s << gap_dynamic << "\n";
   s << kenntext(173) << ":\n";
   s << tau_gap << "\n";
   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
   s << kenntext(180) << ":\n";
   s << randomise_beta_proteins << "\n";
   s << kenntext(181) << ":\n";
   s << randomisation_type << "\n";
   s << kenntext(182) << ":\n";
   s << randomisation_range << "\n";
   s << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

   return s;
}
short betaWerte::fFind(char * parname, ifstream &s, int n) {
   s.close();
   s.open(parname);
   s.clear();
   s.setf(ios::scientific, ios::floatfield);
   short found = 0;
   const char * find = kenntext(n);
   int lang = (int)strlen(find);
   int tmplang;
   int i;
   char d;
   // short stop=0;
   // while (found==0 && s.eof()==0 && stop==0) {
   while (found == 0 && s.eof() == 0) {
      // while (found==0 && stop==0) {
      i = 0;
      d = 'a';
      char tmp[200];
      while (int (d) != 10 && i < 200) {
         s.get(d);
         if (int (d) != 10) {
            tmp[i] = d;
            i++;
         }
      }
      tmplang = i;
      /* Ich verstehe hier zwar nicht warum s.eof()==0 bleibt wenn die Datei zuende
       * ist. Aber mit dem Trick Laenge==200 geht der Abbruch auch! */
      //	if (tmplang==200) stop=1; else {
      if (tmplang >= lang) {
         found = 1;
         for (i = 0; i < lang; i++) {
            if (tmp[i] != find[i]) {
               found = 0;
            }
         }
      }
      // }
   }
   if (found == 0) {
      cout << "X" << n << ",";
   }
   // cout << "\n -> val for "<<kenntext(n)<<" not found! Took standard.\n";
   else {
      cout << ",";
   }
   return found;
}
void betaWerte::fGet(char * parname) {
   // int i;
   ifstream s(parname);

   if (fFind(parname, s, 1) == 1) {
      s >> dt;
   }
   if (fFind(parname, s, 2) == 1) {
      s >> dy;
   }
   if (fFind(parname, s, 3) == 1) {
      s >> t_0;
   }
   if (fFind(parname, s, 4) == 1) {
      s >> t_max;
   }
   if (fFind(parname, s, 5) == 1) {
      s >> dt_output;
   }
   if (fFind(parname, s, 6) == 1) {
      s >> V_0;
   }
   if (fFind(parname, s, 7) == 1) {
      s >> R_bc;
   }
   if (fFind(parname, s, 8) == 1) {
      s >> Sur_ER;
   }
   if (fFind(parname, s, 9) == 1) {
      s >> Vol_ER;
   }
   if (fFind(parname, s, 10) == 1) {
      s >> K_ext;
   }
   if (fFind(parname, s, 11) == 1) {
      s >> Na_ext;
   }
   if (fFind(parname, s, 12) == 1) {
      s >> Ca_ext;
   }
   if (fFind(parname, s, 13) == 1) {
      s >> cal_0;
   }
   if (fFind(parname, s, 14) == 1) {
      s >> K_cal;
   }
   if (fFind(parname, s, 15) == 1) {
      s >> buf_0;
   }
   if (fFind(parname, s, 16) == 1) {
      s >> K_buf;
   }
   if (fFind(parname, s, 17) == 1) {
      s >> glu_0;
   }
   if (fFind(parname, s, 18) == 1) {
      s >> IP3_0;
   }
   if (fFind(parname, s, 19) == 1) {
      s >> C_m;
   }
   if (fFind(parname, s, 26) == 1) {
      s >> K_0;
   }
   if (fFind(parname, s, 27) == 1) {
      s >> Na_0;
   }
   if (fFind(parname, s, 28) == 1) {
      s >> Ca_0;
   }
   if (fFind(parname, s, 29) == 1) {
      s >> T;
   }
   if (fFind(parname, s, 20) == 1) {
      s >> rho[NaK];
   }
   if (fFind(parname, s, 21) == 1) {
      s >> Ihat_NaK;
   }
   if (fFind(parname, s, 22) == 1) {
      s >> H_NaK;
   }
   if (fFind(parname, s, 23) == 1) {
      s >> n_NaK;
   }
   if (fFind(parname, s, 35) == 1) {
      s >> H2_NaK;
   }
   if (fFind(parname, s, 36) == 1) {
      s >> n2_NaK;
   }
   if (fFind(parname, s, 24) == 1) {
      s >> alpha_NaK;
   }
   if (fFind(parname, s, 30) == 1) {
      s >> rho[K_ATP];
   }
   if (fFind(parname, s, 31) == 1) {
      s >> gbar_K_ATP;
   }
   if (fFind(parname, s, 32) == 1) {
      s >> tau_K_ATP;
   }
   if (fFind(parname, s, 33) == 1) {
      s >> s_h_K_ATP;
   }
   if (fFind(parname, s, 34) == 1) {
      s >> kappa_K_ATP;
   }
   if (fFind(parname, s, 40) == 1) {
      s >> rho[K_V];
   }
   if (fFind(parname, s, 41) == 1) {
      s >> gbar_K_V;
   }
   if (fFind(parname, s, 42) == 1) {
      s >> tau_K_V;
   }
   if (fFind(parname, s, 48) == 1) {
      s >> use_dynamic_tau_K_V;
   }
   if (fFind(parname, s, 43) == 1) {
      s >> V_h_K_V;
   }
   if (fFind(parname, s, 44) == 1) {
      s >> kappa_K_V;
   }
   if (fFind(parname, s, 45) == 1) {
      s >> theta_K_V;
   }
   if (fFind(parname, s, 46) == 1) {
      s >> W_h_K_V;
   }
   if (fFind(parname, s, 47) == 1) {
      s >> lambda_K_V;
   }
   if (fFind(parname, s, 50) == 1) {
      s >> rho[K_Ca];
   }
   if (fFind(parname, s, 51) == 1) {
      s >> gbar_K_Ca;
   }
   if (fFind(parname, s, 52) == 1) {
      s >> H_K_Ca;
   }
   if (fFind(parname, s, 53) == 1) {
      s >> n_K_Ca;
   }
   if (fFind(parname, s, 54) == 1) {
      s >> V_h_K_Ca;
   }
   if (fFind(parname, s, 55) == 1) {
      s >> kappa_K_Ca;
   }
   if (fFind(parname, s, 56) == 1) {
      s >> tau_K_Ca;
   }
   if (fFind(parname, s, 57) == 1) {
      s >> use_dynamic_H_K_Ca;
   }
   if (fFind(parname, s, 58) == 1) {
      s >> use_voltage_gating_K_Ca;
   }
   if (fFind(parname, s, 60) == 1) {
      s >> rho[Na_V];
   }
   if (fFind(parname, s, 61) == 1) {
      s >> gbar_Na_V;
   }
   if (fFind(parname, s, 62) == 1) {
      s >> tau_Na_V;
   }
   if (fFind(parname, s, 63) == 1) {
      s >> V_h_Na_V;
   }
   if (fFind(parname, s, 64) == 1) {
      s >> kappa_Na_V;
   }
   if (fFind(parname, s, 65) == 1) {
      s >> theta_Na_V;
   }
   if (fFind(parname, s, 66) == 1) {
      s >> W_h_Na_V;
   }
   if (fFind(parname, s, 67) == 1) {
      s >> lambda_Na_V;
   }
   if (fFind(parname, s, 68) == 1) {
      s >> use_dynamic_tau_Na_V;
   }
   if (fFind(parname, s, 70) == 1) {
      s >> rho[NCX];
   }
   if (fFind(parname, s, 71) == 1) {
      s >> Ihat_NCX;
   }
   if (fFind(parname, s, 72) == 1) {
      s >> H_NCX;
   }
   if (fFind(parname, s, 73) == 1) {
      s >> n_NCX;
   }
   if (fFind(parname, s, 74) == 1) {
      s >> alpha_NCX;
   }
   if (fFind(parname, s, 80) == 1) {
      s >> rho[Ca_L];
   }
   if (fFind(parname, s, 81) == 1) {
      s >> gbar_Ca_L;
   }
   if (fFind(parname, s, 82) == 1) {
      s >> tau_Ca_L;
   }
   if (fFind(parname, s, 83) == 1) {
      s >> V_h_Ca_L;
   }
   if (fFind(parname, s, 84) == 1) {
      s >> kappa_Ca_L;
   }
   if (fFind(parname, s, 85) == 1) {
      s >> theta_Ca_L;
   }
   if (fFind(parname, s, 86) == 1) {
      s >> W_h_Ca_L;
   }
   if (fFind(parname, s, 87) == 1) {
      s >> lambda_Ca_L;
   }
   if (fFind(parname, s, 88) == 1) {
      s >> C_Ca_L;
   }
   if (fFind(parname, s, 89) == 1) {
      s >> n_Ca_L;
   }
   if (fFind(parname, s, 90) == 1) {
      s >> rho[Ca_T];
   }
   if (fFind(parname, s, 91) == 1) {
      s >> gbar_Ca_T;
   }
   if (fFind(parname, s, 92) == 1) {
      s >> tau_Ca_T;
   }
   if (fFind(parname, s, 93) == 1) {
      s >> V_h_Ca_T;
   }
   if (fFind(parname, s, 94) == 1) {
      s >> kappa_Ca_T;
   }
   if (fFind(parname, s, 95) == 1) {
      s >> theta_Ca_T;
   }
   if (fFind(parname, s, 96) == 1) {
      s >> W_h_Ca_T;
   }
   if (fFind(parname, s, 97) == 1) {
      s >> lambda_Ca_T;
   }
   if (fFind(parname, s, 100) == 1) {
      s >> rho[SERCA];
   }
   if (fFind(parname, s, 101) == 1) {
      s >> Ihat_SERCA;
   }
   if (fFind(parname, s, 102) == 1) {
      s >> H_SERCA;
   }
   if (fFind(parname, s, 103) == 1) {
      s >> n_SERCA;
   }
   if (fFind(parname, s, 110) == 1) {
      s >> rho[IP3];
   }
   if (fFind(parname, s, 111) == 1) {
      s >> gbar_IP3;
   }
   if (fFind(parname, s, 112) == 1) {
      s >> g_IP3_max;
   }
   if (fFind(parname, s, 113) == 1) {
      s >> use_dynamic_tau_IP3;
   }
   if (fFind(parname, s, 114) == 1) {
      s >> tau_IP3;
   }
   if (fFind(parname, s, 115) == 1) {
      s >> C_IP3_act;
   }
   if (fFind(parname, s, 116) == 1) {
      s >> n_IP3_act;
   }
   if (fFind(parname, s, 117) == 1) {
      s >> theta_IP3;
   }
   if (fFind(parname, s, 118) == 1) {
      s >> Cbar_IP3_inh;
   }
   if (fFind(parname, s, 119) == 1) {
      s >> n_IP3_inh;
   }
   if (fFind(parname, s, 120) == 1) {
      s >> use_Nernst;
   }
   if (fFind(parname, s, 121) == 1) {
      s >> set_leakage_zero;
   }
   if (fFind(parname, s, 122) == 1) {
      s >> use_inactivation;
   }
   if (fFind(parname, s, 123) == 1) {
      s >> Vbar_Ca_delta;
   }
   if (fFind(parname, s, 124) == 1) {
      s >> buf_ER_0;
   }
   if (fFind(parname, s, 125) == 1) {
      s >> K_buf_ER;
   }
   if (fFind(parname, s, 126) == 1) {
      s >> Ca_ER_0;
   }
   if (fFind(parname, s, 127) == 1) {
      s >> use_dynamic_IP3;
   }
   if (fFind(parname, s, 130) == 1) {
      s >> rho[PMCA];
   }
   if (fFind(parname, s, 131) == 1) {
      s >> Ihat_PMCA;
   }
   if (fFind(parname, s, 132) == 1) {
      s >> H_PMCA;
   }
   if (fFind(parname, s, 133) == 1) {
      s >> n_PMCA;
   }
   if (fFind(parname, s, 134) == 1) {
      s >> alpha_PMCA;
   }
   if (fFind(parname, s, 140) == 1) {
      s >> P_IP3;
   }
   if (fFind(parname, s, 141) == 1) {
      s >> kappa_IP3;
   }
   if (fFind(parname, s, 142) == 1) {
      s >> k_IP3_plus;
   }
   if (fFind(parname, s, 143) == 1) {
      s >> k_IP3_minus;
   }
   if (fFind(parname, s, 144) == 1) {
      s >> C_P;
   }
   if (fFind(parname, s, 145) == 1) {
      s >> n_P;
   }
   if (fFind(parname, s, 150) == 1) {
      s >> rho[sK_Ca];
   }
   if (fFind(parname, s, 151) == 1) {
      s >> gbar_sK_Ca;
   }
   if (fFind(parname, s, 152) == 1) {
      s >> C_sK_Ca;
   }
   if (fFind(parname, s, 153) == 1) {
      s >> kappa_sK_Ca;
   }
   if (fFind(parname, s, 154) == 1) {
      s >> tau_sK_Ca;
   }
   if (fFind(parname, s, 160) == 1) {
      s >> rho[fNa_V];
   }
   if (fFind(parname, s, 161) == 1) {
      s >> gbar_fNa_V;
   }
   if (fFind(parname, s, 162) == 1) {
      s >> tau_fNa_V;
   }
   if (fFind(parname, s, 163) == 1) {
      s >> V_h_fNa_V;
   }
   if (fFind(parname, s, 164) == 1) {
      s >> kappa_fNa_V;
   }
   if (fFind(parname, s, 168) == 1) {
      s >> use_dynamic_tau_fNa_V;
   }
   if (fFind(parname, s, 170) == 1) {
      s >> rho[gap];
   }
   if (fFind(parname, s, 171) == 1) {
      s >> gbar_gap;
   }
   if (fFind(parname, s, 172) == 1) {
      s >> gap_dynamic;
   }
   if (fFind(parname, s, 173) == 1) {
      s >> tau_gap;
   }
   if (fFind(parname, s, 180) == 1) {
      s >> randomise_beta_proteins;
   }
   short aa = 0;
   if (fFind(parname, s, 181) == 1) {
      s >> aa;
   }
   randomisation_type = random_laws(aa);
   if (fFind(parname, s, 182) == 1) {
      s >> randomisation_range;
   }

   s.close();
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Parameter-Schreib-und Lese Prozeduren:
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void hyphasmaParameter::goon() {
   // Das Auskommentierte wird unter Unix nicht compiliert!!!
   //   int hit=0;
   cout << "Weiter mit der Eingabetaste ...";
   //   while (hit==0) { hit=kbhit(); }
   //   hit=getch();
   char bla = 0;
   cin >> bla;
   cout << "\n";
}
void hyphasmaParameter::read(bool transform2rate) {

   //cout << "Data File " << namesuffix << " ... ";
   Value.fGet(namesuffix, transform2rate);

   cout << "   ... Reading finished!\n\n";
   if (Value.show_mode == islet) {
      char name2suffix[namelength] = "betacell.par";
      cout << "Data " << name2suffix << " ... ";
      betaValue.fGet(name2suffix);
      cout << "   ... Reading finished!\n";
   }
   // ##### Note that no write procedure is provided here.
   // ##### The handling is not adapted to betacells.
   // ##### using transform2rate is not thought through for betacells
}

///§§§ Philippe
void hyphasmaParameter::write(string filename) {
   if(filename.size() ==0) filename = string(namesuffix);
   ofstream oParam(filename.c_str());
   oParam.setf(ios::scientific, ios::floatfield);
   Value.fPut(oParam);
   oParam.close();
   cout << "Parametersatz in der Datei " << filename << " gespeichert.\n";
}

void hyphasmaParameter::save() {
   strcpy(ntmp, name);
   vtmp = Value;
}
void hyphasmaParameter::recover() {
   strcpy(name, ntmp);
   Value = vtmp;
}
void hyphasmaParameter::in_datei(suffix suff, suffix log) {
   cout << "Dateiname der Parameterdatei (ohne Suffix): ";
   cin >> name;
   strcpy(namesuffix, name);
   strcat(namesuffix, suff);
   strcpy(logfile, name);
   strcat(logfile, log);
   cout << "Datei: " << namesuffix << "\n";
   // goon();
}
//void hyphasmaParameter::reset_random_ini(const char * namecommand, int value) {
//   // Erweiterungen:
//   suffix suff = ".par";
//   // name="";
//   // strcpy(name,namecommand);
//   strcpy(namesuffix, namecommand);
//   strcat(namesuffix, suff);
//   read(false);
//   Value.ini_random = value;
//   cout << "Reset the random number generator seed value to " << value << ".\n";
//   write();
//}
int hyphasmaParameter::wahl(const char * namecommand, bool transform2rate, bool show_missing) {
   // Lokale Variablen:
   int menu = 1;
   char runfree = 0;
   // Default Dateiname:
   char namedefault[namelength] = "standard";
   // Erweiterungen:
   suffix suff = ".par";
   suffix log = ".log";
   Value.show_missing_pars = show_missing;

   while (menu != 0 && menu != 5) {
      if (namecommand[0] == '#') {
         cout << "\n";
         cout << "=========================================\n";
         cout << "-------- Parameterdefinition: -----------\n";
         cout << "=========================================\n\n";
         cout << "\n";
         cout << "0: Ende\n";
         cout << "1: 2D-Parameterdatei erzeugen\n";
         cout << "2: 3D-Parameterdatei erzeugen\n";
         cout << "4: Parameterdatei laden\n";
         cout << "5: Rechnung starten\n";
         //   cout << "6: Rechnung starten (mit log-file)\n";
         cin >> menu;
      } else {
         menu = 5;
      }
       
      if (menu != 0) {
         switch (menu) {
            case 1: {
               // write 2D parameter file with default values
               in_datei(suff, log);
               Value.ini2d();
               write();
               runfree = 1;
                break;
            }

            case 2: {
               // write 3D parameter file with default values
               in_datei(suff, log);
               Value.ini3d();
               write();
               runfree = 1;
                break;
            }

            case 3: {
               // what is this for???
               save(); // memorize current parameters
               in_datei(suff, log); // get a filename from console
               read(true); // read this file
               recover(); // restore previous parameters
                break;
            }

            case 4: {
               // read a parameter file determined at the console
               in_datei(suff, log);
               read(true);
               runfree = 1;
                break;
            }

            case 5: {
               // read a parameter file (either namedefault or namecommand)
               if (runfree == 0) {
                  if (namecommand[0] == '#') {
                     strcpy(name, namedefault);
                     strcpy(namesuffix, namedefault);
                     strcat(namesuffix, suff);
                     strcpy(logfile, namedefault);
                     strcat(logfile, log);
                  } else {
                     if (string(namecommand).size() >= namelength - 1) {
                        cerr << "ERR: file path/name exceeds namelength : " << endl
                             << namecommand << endl;
                        exit(-1);
                     } /// Philippe, got a weird seg fault when this happens
                     cout << "   -> Reading Parameter file :" << name << " " << namecommand << "(.par)\n";
                     strcpy(name, namecommand);
                     strcpy(namesuffix, namecommand);
                     strcat(namesuffix, suff);
                     strcpy(logfile, namecommand);
                     strcat(logfile, log);
                  }
                  read(transform2rate);
               }
               cout << "   -> Now starting the simulation\n";
               //cout << "Programmlauf mit dem Parametersatz  ";
               //cout << namesuffix << "  !\n";
               // cout << "Keine Dokumentation in einer log-Datei.\n";
               // Value.doku=0;
                break;

            }
               /*     case 6:
                *     if (runfree==0) {
                *       strcpy(name,namedefault);
                *       strcpy(namesuffix,namedefault);
                *       strcat(namesuffix,suff);
                *       strcpy(logfile,namedefault);
                *       strcat(logfile,log);
                *       read(transform2rate);
                *     }
                *     cout << "Programmlauf mit dem Parametersatz  ";
                *     cout << namesuffix << "  !\n";
                *     cout << "Mit Dokumentation in der Datei:  ";
                *     cout << logfile << "  .\n";
                *     // Standardeinstellung fuer doku verwenden
                *     menu=5;
                *     break; */
         }     // switch
      }        // menu!=0
   }           // menu!=0 && menu!=5
   return menu;
}









///§§§ Philippe !!!!!!!!!!!!!!!!!!!!



void compare(short s1, short s2, const char* message){
    if(s1 != s2) cout << message << ":" << s1 << " != " << s2 << endl;
}
void compare(int s1, int s2, const char* message){
    if(s1 != s2) cout << message<< ":" << s1 << " != " << s2<< endl;
}
void compare(double s1, double s2, const char* message){
    if(s1 != s2) cout << message<< ":" << s1 << " != " << s2<< endl;
}
void compare(bool s1, bool s2, const char* message){
    if(s1 != s2) cout << message<< ":" << s1 << " != " << s2<< endl;
}
void compare(long s1, long s2, const char* message){
    if(s1 != s2) cout << message<< ":" << s1 << " != " << s2<< endl;
}
void compare(unsigned short s1, unsigned short s2, const char* message){
    if(s1 != s2) cout << message<< ":" << s1 << " != " << s2<< endl;
}

void compare(vector<int> s1, vector<int> s2, const char* message){
    int cpt = 0;
    if(s1.size() != s2.size()) cout << message<< "unmatched size" << endl;
    for(int i = 0; i < (int) s1.size(); ++i){
        if(s1[i] != s2[i]){
            cout << message<< "[" << i << "]:" << s1[i] << " != " << s2[i]<< endl;
            cpt++;
        }
        if(cpt >= 10) return;
    }
}

void compare(vector<double> s1, vector<double> s2, const char* message){
    int cpt = 0;
    if(s1.size() != s2.size()) cout << message<< "unmatched size" << endl;
    for(int i = 0; i < (int) s1.size(); ++i){
        if(s1[i] != s2[i]){
            cout << message<< "[" << i << "]:" << s1[i] << " != " << s2[i]<< endl;
            cpt++;
        }
        if(cpt >= 10) return;
    }
}

void compare(vector<string> &s1, vector<string> &s2, const char* message){
    int cpt = 0;
    if(s1.size() != s2.size()) cout << message<< "unmatched size" << endl;
    for(int i = 0; i < (int) s1.size(); ++i){
        if(s1[i] != s2[i]){
            cout << message<< "[" << i << "]:" << s1[i] << " != " << s2[i]<< endl;
            cpt++;
        }
        if(cpt >= 10) return;
    }
}

void compare(dynarray<long int> &s1, dynarray<long int> &s2, const char* message){
    if(s1.Vergleich(s2)) return;
    int cpt = 0;
    if(s1.benutzt() != s2.benutzt()) cout << message<< "unmatched size" << endl;
    for(int i = 0; i < s1.benutzt(); ++i){
        if(s1[i] != s2[i]){
            cout << message<< "[" << i << "]:" << s1[i] << " != " << s2[i]<< endl;
            cpt++;
        }
        if(cpt >= 10) return;
    }
}

void compare(dynarray<bool> &s1, dynarray<bool> &s2, const char* message){
    if(s1.Vergleich(s2)) return;
    int cpt = 0;
    if(s1.benutzt() != s2.benutzt()) cout << message<< "unmatched size" << endl;
    for(int i = 0; i < s1.benutzt(); ++i){
        if(s1[i] != s2[i]){
            cout << message<< "[" << i << "]:" << s1[i] << " != " << s2[i]<< endl;
            cpt++;
        }
         if(cpt >= 10) return;
    }
}

void Werte::compareParameterSets( Werte &w1, Werte &w2){
   compare(w1.show_missing_pars, w2.show_missing_pars, "show_missing_pars");
   compare(w1.system, w2.system, "system");
   compare(w1.outputfiles, w2.outputfiles, "outputfiles");
   compare(w1.timevalues, w2.timevalues, "timevalues");
   compare(w1.show_Ki67, w2.show_Ki67, "show_Ki67");
   compare(w1.safety_checks, w2.safety_checks, "safety_checks");
   compare(w1.show_mode, w2.show_mode, "show_mode");
   compare(w1.ini_random, w2.ini_random, "ini_random");
   compare(w1.late_ini_random, w2.late_ini_random, "late_ini_random");
   compare(w1.file_output, w2.file_output, "file_output");
   compare(w1.CB_Narray, w2.CB_Narray, "CB_Narray");
   compare(w1.CC_Narray, w2.CC_Narray, "CC_Narray");
   compare(w1.TC_Narray, w2.TC_Narray, "TC_Narray");
   compare(w1.FDC_Narray, w2.FDC_Narray, "FDC_Narray");
   compare(w1.OUT_Narray, w2.OUT_Narray, "OUT_Narray");
   compare(w1.STROMA_Narray, w2.STROMA_Narray, "STROMA_Narray");
   compare(w1.BETA_Narray, w2.BETA_Narray, "BETA_Narray");

   // Shape space:
   // ============
   compare(w1.DimShapeSpace, w2.DimShapeSpace, "DimShapeSpace");
   compare(w1.metrik, w2.metrik, "metrik");
   compare(w1.SSStates, w2.SSStates, "SSStates");
   compare(w1.SSRangePerDim, w2.SSRangePerDim, "SSRangePerDim");
   compare(w1.totalA, w2.totalA, "totalA");
   compare(w1.APeakNumber, w2.APeakNumber, "APeakNumber");
   compare(w1.takeA, w2.takeA, "takeA");
   compare(w1.ag_fraction, w2.ag_fraction, "ag_fraction");
   compare(w1.GammaGauss,w2.GammaGauss, "GammaGauss");
   compare(w1.amplitudeGauss, w2.amplitudeGauss, "amplitudeGauss");
   compare(w1.type_affinity_function,w2.type_affinity_function,"type_affinity_function");
   compare(w1.use_logarithmic_seq_affinity,w2.use_logarithmic_seq_affinity,"use_logarithmic_seq_affinity");
   compare(w1.use_sequence_space,w2.use_sequence_space,"use_sequence_space");
   compare(w1.size_sequences,w2.size_sequences,"size_sequences");
   compare(w1.sequence_mut_per_base, w2.sequence_mut_per_base,"sequence_mut_per_base");
   compare(w1.init_antigen_sequences,w2.init_antigen_sequences,"init_antigen_sequences");
   compare(w1.initAntigenSeqs, w2.initAntigenSeqs,"initAntigenSeqs");
   compare(w1.max_hamming_antigens,w2.max_hamming_antigens,"max_hamming_antigens");
   compare(w1.min_hamming_antigens,w2.min_hamming_antigens,"min_hamming_antigens");
   compare(w1.initBCRSeqs, w2.initBCRSeqs,"initBCRSeqs");
   compare(w1.max_hamming_BCRs,w2.max_hamming_BCRs,"max_hamming_BCRs");
   compare(w1.min_initial_affinity_BCRs,w2.min_initial_affinity_BCRs,"min_initial_affinity_BCRs");
   compare(w1.max_initial_affinity_BCRs,w2.max_initial_affinity_BCRs,"max_initial_affinity_BCRs");
   compare(w1.initTCRSeqs,w2.initTCRSeqs,"initTCRSeqs");
   compare(w1.max_hamming_TCRs,w2.max_hamming_TCRs,"max_hamming_TCRs");
   compare(w1.min_initial_affinity_TCRs,w2.min_initial_affinity_TCRs,"min_initial_affinity_TCRs");
   compare(w1.max_initial_affinity_TCRs,w2.max_initial_affinity_TCRs,"max_initial_affinity_TCRs");
   compare(w1.R_affinity,w2.R_affinity,"R_affinity");
   compare(w1.max_affinity_cluster,w2.max_affinity_cluster,"max_affinity_cluster");

   // Arup space:        // Philippe 20-03-2016
   // ============
   compare(w1.use_arup_space,w2.use_arup_space,"use_arup_space");
   compare(w1.arup_length_sequences,w2.arup_length_sequences,"arup_length_sequences");
   compare(w1.arup_N_conserved,w2.arup_N_conserved,"arup_N_conserved");
   compare(w1.arup_N_mutates,w2.arup_N_mutates,"arup_N_mutates");
   compare(w1.arup_N_shielded,w2.arup_N_shielded,"arup_N_shielded");
   compare(w1.arup_nb_ini_antigens,w2.arup_nb_ini_antigens,"arup_nb_ini_antigens");
   compare(w1.arup_ini_antigens,w2.arup_ini_antigens,"arup_ini_antigens");
   compare(w1.arup_ag_fraction,w2.arup_ag_fraction,"arup_ag_fraction");
   compare(w1.arup_nb_mutations_gen_strains,w2.arup_nb_mutations_gen_strains,"arup_nb_mutations_gen_strains");
   compare(w1.arup_threshold_activation,w2.arup_threshold_activation,"arup_threshold_activation");
   compare(w1.arup_h_min,w2.arup_h_min,"arup_h_min");
   compare(w1.arup_h_max,w2.arup_h_max,"arup_h_max");
   compare(w1.arup_ini_bcrs,w2.arup_ini_bcrs,"arup_ini_bcrs");
   compare(w1.arup_mutation,w2.arup_mutation,"arup_mutation");
   compare(w1.arup_proba_lethal_mut,w2.arup_proba_lethal_mut,"arup_proba_lethal_mut");
   compare(w1.arup_proba_affecting_mut,w2.arup_proba_affecting_mut,"arup_proba_affecting_mut");
   compare(w1.arup_proba_silent_mut,w2.arup_proba_silent_mut,"arup_proba_silent_mut");
   compare(w1.arup_law_mut_Xs,w2.arup_law_mut_Xs,"arup_law_mut_Xs");
   compare(w1.arup_law_mut_Densities,w2.arup_law_mut_Densities,"arup_law_mut_Densities");
   compare(w1.arup_alpha,w2.arup_alpha,"arup_alpha");
   compare(w1.arup_hprime_min,w2.arup_hprime_min,"arup_hprime_min");
   compare(w1.arup_hprime_max,w2.arup_hprime_max,"arup_hprime_max");
   compare(w1.arup_hmut_min,w2.arup_hmut_min,"arup_hmut_min");
   compare(w1.arup_hmut_max,w2.arup_hmut_max,"arup_hmut_max");

   // Signals:
   // ========
   compare(w1.signal_mode,w2.signal_mode,"signal_mode");
   compare(w1.objects_transparent,w2.objects_transparent,"objects_transparent");
   compare(w1.bound_differ2CC,w2.bound_differ2CC,"bound_differ2CC");
   compare(w1.bound_CXCL12,w2.bound_CXCL12,"bound_CXCL12");
   compare(w1.bound_CXCL13,w2.bound_CXCL13,"bound_CXCL13");
   compare(w1.bound_ab,w2.bound_ab,"bound_ab");
   compare(w1.bound_ag,w2.bound_ag,"bound_ag");
   compare(w1.bound_SEMA4D,w2.bound_SEMA4D,"bound_SEMA4D");
   compare(w1.CXCL12crit,w2.CXCL12crit,"CXCL12crit");
   compare(w1.CXCL13crit,w2.CXCL13crit,"CXCL13crit");
   compare(w1.CXCL12recrit,w2.CXCL12recrit,"CXCL12recrit");
   compare(w1.CXCL13recrit,w2.CXCL13recrit,"CXCL13recrit");
   compare(w1.D_differ2CC,w2.D_differ2CC,"D_differ2CC");
   compare(w1.D_CXCL12,w2.D_CXCL12,"D_CXCL12");
   compare(w1.D_CXCL13,w2.D_CXCL13,"D_CXCL13");
   compare(w1.D_antibody,w2.D_antibody,"D_antibody");
   compare(w1.D_antigen,w2.D_antigen,"D_antigen");
   compare(w1.D_SEMA4D,w2.D_SEMA4D,"D_SEMA4D");
   compare(w1.fix_signals,w2.fix_signals,"fix_signals");
   compare(w1.fix_glucose_gradient,w2.fix_glucose_gradient,"fix_glucose_gradient");
   compare(w1.const_dynamic_glucose_field,w2.const_dynamic_glucose_field,"const_dynamic_glucose_field");
   compare(w1.fix_glucose_gradient_min,w2.fix_glucose_gradient_min,"fix_glucose_gradient_min");
   compare(w1.fix_glucose_gradient_max,w2.fix_glucose_gradient_max,"fix_glucose_gradient_max");

   // Space properties:
   // =================
   compare(w1.DimSpace,w2.DimSpace,"DimSpace");
   compare(w1.GC_radius,w2.GC_radius,"GC_radius");
   compare(w1.gridsize[0],w2.gridsize[0],"gridsize[0]");
   compare(w1.gridsize[1],w2.gridsize[1],"gridsize[1]");
   compare(w1.gridsize[2],w2.gridsize[2],"gridsize[2]");
   compare(w1.dx,w2.dx,"dx");
   compare(w1.dx_signal,w2.dx_signal,"dx_signal");
   compare(w1.vol_shape,w2.vol_shape,"vol_shape");
   compare(w1.obstacles,w2.obstacles,"obstacles");
   compare(w1.wall_level,w2.wall_level,"wall_level");
   compare(w1.collagen_density,w2.collagen_density,"collagen_density");
   compare(w1.collagen_cluster,w2.collagen_cluster,"collagen_cluster");
   compare(w1.wall_width,w2.wall_width,"wall_width");
   compare(w1.slit_number,w2.slit_number,"slit_number");
   compare(w1.slit_width,w2.slit_width,"slit_width");

   // Time and phases:
   // ================
   compare(w1.Start_Mutation,w2.Start_Mutation,"Start_Mutation");
   compare(w1.Start_Differentiation,w2.Start_Differentiation,"Start_Differentiation");
   compare(w1.StartOutput,w2.StartOutput,"StartOutput");
   compare(w1.newBCinflux_stop,w2.newBCinflux_stop,"newBCinflux_stop");
   compare(w1.deltat,w2.deltat,"deltat");
   compare(w1.tmin,w2.tmin,"tmin");
   compare(w1.tmax,w2.tmax,"tmax");
   compare(w1.ToFileStep,w2.ToFileStep,"ToFileStep");

   // Cells in general
   // ================
   compare(w1.adhesion_time,w2.adhesion_time,"adhesion_time");
   compare(w1.chemo_max,w2.chemo_max,"chemo_max");
   compare(w1.chemo_steep,w2.chemo_steep,"chemo_steep");
   compare(w1.chemo_half,w2.chemo_half,"chemo_half");
   compare(w1.use_glucose,w2.use_glucose,"use_glucose");
   compare(w1.use_oxygen,w2.use_oxygen,"use_oxygen");
   compare(w1.use_glucose_pro,w2.use_glucose_pro,"use_glucose_pro");
   compare(w1.use_oxygen_pro,w2.use_oxygen_pro,"use_oxygen_pro");
   compare(w1.critical_nutrient,w2.critical_nutrient,"critical_nutrient");
   compare(w1.bound_glucose,w2.bound_glucose,"bound_glucose");
   compare(w1.bound_oxygen,w2.bound_oxygen,"bound_oxygen");
   compare(w1.D_glucose,w2.D_glucose,"D_glucose");
   compare(w1.D_glucose_H2O,w2.D_glucose_H2O,"D_glucose_H2O");
   compare(w1.D_oxygen,w2.D_oxygen,"D_oxygen");
   compare(w1.D_oxygen_H2O,w2.D_oxygen_H2O,"D_oxygen_H2O");
   compare(w1.p_macrophage,w2.p_macrophage,"p_macrophage");
   compare(w1.allow_exchange,w2.allow_exchange,"allow_exchange");    // Makes exchange of contact inhibited cells possible >=v7.05.1
   compare(w1.use_specific_turning_angles,w2.use_specific_turning_angles,"use_specific_turning_angles");

   // Centroblasts or blast1:
   // =======================
   compare(w1.totalB,w2.totalB,"totalB");
   compare(w1.totalBss,w2.totalBss,"totalBss");
   compare(w1.newBCinflux_rate,w2.newBCinflux_rate,"newBCinflux_rate");
   compare(w1.smooth_stopBCinflux,w2.smooth_stopBCinflux,"smooth_stopBCinflux");
   compare(w1.CB_radius,w2.CB_radius,"CB_radius");
   compare(w1.takeB,w2.takeB,"takeB");
   compare(w1.min_seeder_dist,w2.min_seeder_dist,"min_seeder_dist");
   compare(w1.max_seeder_dist,w2.max_seeder_dist,"max_seeder_dist");
   compare(w1.posCB,w2.posCB,"posCB");
   compare(w1.proliferate,w2.proliferate,"proliferate");
   compare(w1.dx_CB,w2.dx_CB,"dx_CB");       // maximal distance for CB-proliferation from dividing cell
   compare(w1.grow,w2.grow,"grow");
   compare(w1.tolight,w2.tolight,"tolight");
   compare(w1.CB_maxvolume4differ2CC,w2.CB_maxvolume4differ2CC,"CB_maxvolume4differ2CC");
   compare(w1.mutation,w2.mutation,"mutation");
   compare(w1.mutation_after_tc,w2.mutation_after_tc,"mutation_after_tc");
   compare(w1.mutation_after_dec_tc,w2.mutation_after_dec_tc,"mutation_after_dec_tc");
   compare(w1.mutation_affinity_exponent,w2.mutation_affinity_exponent,"mutation_affinity_exponent");
   compare(w1.CB_fixed_times_of_divisions,w2.CB_fixed_times_of_divisions,"CB_fixed_times_of_divisions");
   compare(w1.CB_fixed_times_of_divisions_in_expansion,w2.CB_fixed_times_of_divisions_in_expansion,"CB_fixed_times_of_divisions_in_expansion");
   compare(w1.CB2OUT_prob,w2.CB2OUT_prob,"CB2OUT_prob");
   compare(w1.reset_antigen_after_collection,w2.reset_antigen_after_collection,"reset_antigen_after_collection");
   compare(w1.present_specific_ag2TC,w2.present_specific_ag2TC,"present_specific_ag2TC");
   compare(w1.fixed_time_of_divisions_mode,w2.fixed_time_of_divisions_mode,"fixed_time_of_divisions_mode");
   compare(w1.smooth_differentiation,w2.smooth_differentiation,"smooth_differentiation");
   compare(w1.smooth_dif2out,w2.smooth_dif2out,"smooth_dif2out");
   compare(w1.exit2tz,w2.exit2tz,"exit2tz");
   compare(w1.smooth_differentiation_time,w2.smooth_differentiation_time,"smooth_differentiation_time");
   compare(w1.smooth_dif2out_time,w2.smooth_dif2out_time,"smooth_dif2out_time");
   compare(w1.CB_dt_G0,w2.CB_dt_G0,"CB_dt_G0");
   compare(w1.CB_dt_G1,w2.CB_dt_G1,"CB_dt_G1");
   compare(w1.CB_dt_G2,w2.CB_dt_G2,"CB_dt_G2");
   compare(w1.CB_dt_S,w2.CB_dt_S,"CB_dt_S");
   compare(w1.CB_dt_M,w2.CB_dt_M,"CB_dt_M");
   compare(w1.CB_dtphase_width,w2.CB_dtphase_width,"CB_dtphase_width");
   compare(w1.t_inject_BrdU, w2.t_inject_BrdU,"t_inject_BrdU");
   compare(w1.deltat_inject_BrdU,w2.deltat_inject_BrdU,"deltat_inject_BrdU");
   compare(w1.BrdU_detection_threshold,w2.BrdU_detection_threshold,"BrdU_detection_threshold");
   compare(w1.n_inject_BrdU,w2.n_inject_BrdU,"n_inject_BrdU");
   compare(w1.transmit_CC_delay_to_CB_cycle,w2.transmit_CC_delay_to_CB_cycle,"transmit_CC_delay_to_CB_cycle");
   compare(w1.retain_ag,w2.retain_ag,"retain_ag");
   compare(w1.ag_loaded_CB_diff2output,w2.ag_loaded_CB_diff2output,"ag_loaded_CB_diff2output");
   compare(w1.ag_loaded_CC_directly2TFH,w2.ag_loaded_CC_directly2TFH,"ag_loaded_CC_directly2TFH");
   compare(w1.ag_loaded_CB_stop_mutation,w2.ag_loaded_CB_stop_mutation,"ag_loaded_CB_stop_mutation");
   compare(w1.ag_deleted_in_fresh_CC,w2.ag_deleted_in_fresh_CC,"ag_deleted_in_fresh_CC");
   compare(w1.divide_ag_asymmetric,w2.divide_ag_asymmetric,"divide_ag_asymmetric");
   compare(w1.asymmetric_polarity_index,w2.asymmetric_polarity_index,"asymmetric_polarity_index");
   compare(w1.smooth_PI,w2.smooth_PI,"smooth_PI");
   compare(w1.BC_ag_preloaded,w2.BC_ag_preloaded,"BC_ag_preloaded");
   compare(w1.CBreceptor_use,w2.CBreceptor_use,"CBreceptor_use");
   compare(w1.CBreceptor_dissociation,w2.CBreceptor_dissociation,"CBreceptor_dissociation");
   compare(w1.CBreceptor_binding,w2.CBreceptor_binding,"CBreceptor_binding");
   compare(w1.CBreceptor_total,w2.CBreceptor_total,"CBreceptor_total");
   compare(w1.CBreceptor_activation,w2.CBreceptor_activation,"CBreceptor_activation");
   compare(w1.CB_max_adhesion,w2.CB_max_adhesion,"CB_max_adhesion");    // maximum adhesion force in % of full stickness
   compare(w1.D_CB,w2.D_CB,"D_CB");  // Diffusion (alternative to cell velocity -- if v is used take -1)
   compare(w1.v_CB,w2.v_CB,"v_CB");             // CB-velocity (alternative to diffusion -- if D is used take -1)
   compare(w1.v_CB_width,w2.v_CB_width,"v_CB_width");       // defines a width of Gauss distributed v_CB values (-1 for fixed)
   compare(w1.CB_smoothmove,w2.CB_smoothmove,"CB_smoothmove");    // Distributes a barycenter movement thought to overcome one
   compare(w1.CB_persistence,w2.CB_persistence,"CB_persistence");   // average time gap in minutes between changes of direction of the
   compare(w1.CB_v_modi,w2.CB_v_modi,"CB_v_modi");         // modus of velocity state treatment
   compare(w1.CB_n_v_states,w2.CB_n_v_states,"CB_n_v_states");     // # of velocity states
   compare(w1.v_CB_switch_deltat,w2.v_CB_switch_deltat,"v_CB_switch_deltat");    // Mean duration in a v-state in minutes
   compare(w1.v_CB_factor,w2.v_CB_factor,"v_CB_factor");      // for 2 velocities: the factor by which the velocity is reduced
   compare(w1.distance_tolerance,w2.distance_tolerance,"distance_tolerance");
   compare(w1.half_tolerance_deformation,w2.half_tolerance_deformation,"half_tolerance_deformation");
   compare(w1.CB_D_cytosol,w2.CB_D_cytosol,"CB_D_cytosol");     // Diffusion constant for fragments in the cytosol.
   compare(w1.v_CB_cytosol,w2.v_CB_cytosol,"v_CB_cytosol");     // Alternative to CB_D_cytosol (one of both has to be set to -1).
   compare(w1.CB_elongation,w2.CB_elongation,"CB_elongation");    // Cell elongation by active movement
   compare(w1.CB_K_elongation,w2.CB_K_elongation,"CB_K_elongation");    // Elongation in units of spherical cell radius, at which

   // blast2:
   // =======================
   compare(w1.total_blast2,w2.total_blast2,"total_blast2");
   compare(w1.blast2_radius,w2.blast2_radius,"blast2_radius");
   compare(w1.pos_blast2,w2.pos_blast2,"pos_blast2");
   compare(w1.dx_blast2,w2.dx_blast2,"dx_blast2");
   compare(w1.blast2_proliferate,w2.blast2_proliferate,"blast2_proliferate");
   compare(w1.blast2_grow,w2.blast2_grow,"blast2_grow");
   compare(w1.blast2_distance_tolerance,w2.blast2_distance_tolerance,"blast2_distance_tolerance");
   compare(w1.blast2_half_tolerance_deformation,w2.blast2_half_tolerance_deformation,"blast2_half_tolerance_deformation");
   compare(w1.D_blast2,w2.D_blast2,"D_blast2");

   // Centrocytes:
   // ============
   compare(w1.CC_test_delay,w2.CC_test_delay,"CC_test_delay");
   compare(w1.CC_ICAM_delay,w2.CC_ICAM_delay,"CC_ICAM_delay");
   compare(w1.CC_FDC_selection,w2.CC_FDC_selection,"CC_FDC_selection");
   compare(w1.collectFDCsignals,w2.collectFDCsignals,"collectFDCsignals");
   compare(w1.collectFDCperiod,w2.collectFDCperiod,"collectFDCperiod");
   compare(w1.prob2kill_noFDCcontactBCs,w2.prob2kill_noFDCcontactBCs,"prob2kill_noFDCcontactBCs");
   compare(w1.TCell,w2.TCell,"TCell");
   compare(w1.output,w2.output,"output");
   compare(w1.output_DEC,w2.output_DEC,"output_DEC");
   compare(w1.FDCsignalling,w2.FDCsignalling,"FDCsignalling");
   compare(w1.shrink,w2.shrink,"shrink");
   compare(w1.apoptosis,w2.apoptosis,"apoptosis");
   compare(w1.apoptosis4FDCselected,w2.apoptosis4FDCselected,"apoptosis4FDCselected");
   compare(w1.macrophage,w2.macrophage,"macrophage");
   compare(w1.ignore_affinity,w2.ignore_affinity,"ignore_affinity");
   compare(w1.selection,w2.selection,"selection");
   compare(w1.ccdiff,w2.ccdiff,"ccdiff");
   compare(w1.ccdiff_delay,w2.ccdiff_delay,"ccdiff_delay");
   compare(w1.ccdiff_delay_DEC,w2.ccdiff_delay_DEC,"ccdiff_delay_DEC");
   compare(w1.final_differentiation_rate,w2.final_differentiation_rate,"final_differentiation_rate");
   compare(w1.tc_search_duration_per_FDCcontact,w2.tc_search_duration_per_FDCcontact,"tc_search_duration_per_FDCcontact"); // add this to the search time for each collected
   compare(w1.tc_search_duration_fixed,w2.tc_search_duration_fixed,"tc_search_duration_fixed"); // BC search duration for Tfh for mode==1
   compare(w1.TC_time,w2.TC_time,"TC_time");          // Duration of TC-CC interaction in hours
   compare(w1.TC_time_width,w2.TC_time_width,"TC_time_width");    // width of this duration in hours
   compare(w1.TC_rescue_time,w2.TC_rescue_time,"TC_rescue_time");   // Minimum duration for selection
   compare(w1.mode_of_setting_TC_time,w2.mode_of_setting_TC_time,"mode_of_setting_TC_time");
   compare(w1.tc_search_duration_mode,w2.tc_search_duration_mode,"tc_search_duration_mode");
   compare(w1.negativeTCselection,w2.negativeTCselection,"negativeTCselection");
   compare(w1.BCstaysonTCbyTCtime,w2.BCstaysonTCbyTCtime,"BCstaysonTCbyTCtime");
   compare(w1.multipleTFHcontacts,w2.multipleTFHcontacts,"multipleTFHcontacts");
   compare(w1.D_CC,w2.D_CC,"D_CC");             // Diffusion (alternative to cell velocity -- if v is used take -1)
   compare(w1.CXCR5down,w2.CXCR5down,"CXCR5down");        // rate of CXCR5 downregulation (-1 for none)
   compare(w1.CXCR4down,w2.CXCR4down,"CXCR4down");        // rate of CXCR4 downregulation (-1 for none)
   compare(w1.v_CC,w2.v_CC,"v_CC");             // CB-velocity (alternative to diffusion -- if D is used take -1)
   compare(w1.v_CC_width,w2.v_CC_width,"v_CC_width");       // defines a width of Gauss distributed v_CC values (-1 for fixed)
   compare(w1.CC_persistence,w2.CC_persistence,"CC_persistence");   // average time gap in minutes between changes of direction of the
   compare(w1.CC_v_modi,w2.CC_v_modi,"CC_v_modi");         // modus of velocity state treatment
   compare(w1.CC_n_v_states,w2.CC_n_v_states,"CC_n_v_states");     // # of velocity states
   compare(w1.v_CC_switch_deltat,w2.v_CC_switch_deltat,"v_CC_switch_deltat");    // Mean duration in a v-state in minutes
   compare(w1.v_CC_factor,w2.v_CC_factor,"v_CC_factor");      // for 2 velocities: the factor by which the velocity is reduced
   compare(w1.use_ab_dynamics,w2.use_ab_dynamics,"use_ab_dynamics");    // =0 for old affinity model;
   compare(w1.initial_ab_affinity,w2.initial_ab_affinity,"initial_ab_affinity");    // -1 for take seeder cell average
   compare(w1.ignore_apoptotic_CC,w2.ignore_apoptotic_CC,"ignore_apoptotic_CC");    // 0: include them; 1: ignore them for cell number analysis
   compare(w1.CC_apoptotic_motility_mode,w2.CC_apoptotic_motility_mode,"CC_apoptotic_motility_mode");
   compare(w1.p_apo_randomwalk,w2.p_apo_randomwalk,"p_apo_randomwalk");

   // class switch
   // ============
   compare(w1.do_switch_classes,w2.do_switch_classes,"do_switch_classes");
   compare(w1.switch_dimension,w2.switch_dimension,"switch_dimension");
   for(int i = 0; i < w1.switch_dimension; ++i){
       compare(w1.switch_matrix[i],w2.switch_matrix[i],"w1.switch_matrix[i=?]");
   }
   compare(w1.IgE_BCRlevel,w2.IgE_BCRlevel,"IgE_BCRlevel");

   ///§§§§ Philippe
   compare(w1.IgG_BCRlevel,w2.IgG_BCRlevel,"IgG_BCRlevel");


   compare(w1.IgE_factor_cellcycle,w2.IgE_factor_cellcycle,"IgE_factor_cellcycle");
   compare(w1.IgE_factor_divisions,w2.IgE_factor_divisions,"IgE_factor_divisions");
   compare(w1.CC_IgE_prob_CXCR5down,w2.CC_IgE_prob_CXCR5down,"CC_IgE_prob_CXCR5down");

   ///§§§Philippe
   compare(w1.IgG_BCRlevel,w2.IgG_BCRlevel,"IgG_BCRlevel");             // 373
   compare(w1.IgG_factor_cellcycle,w2.IgG_factor_cellcycle,"IgG_factor_cellcycle");     // 374
   compare(w1.IgG_factor_divisions,w2.IgG_factor_divisions,"IgG_factor_divisions");     // 375
   for(int i = 0; i < nIg_classes; ++i){
       compare(w1.Founder_IgX[i],w2.Founder_IgX[i],"w1.Founder_IgX[i=?]");
   }
   compare(w1.decay_proba_switch, w2.decay_proba_switch,"decay_proba_switch");       // 377 for all probabilities of switching
   compare(w1.IgG_factor_leaving,w2.IgG_factor_leaving,"IgG_factor_leaving");       // 378
   compare(w1.Affinity_threshold_IgG,w2.Affinity_threshold_IgG,"Affinity_threshold_IgG");   // 379
   compare(w1.Int_Antigen_threshold_IgG,w2.Int_Antigen_threshold_IgG,"Int_Antigen_threshold_IgG");// 380
   compare(w1.tc_help_IgG,w2.tc_help_IgG,"tc_help_IgG");              // 381
   compare(w1.stddev_initial_divisions,w2.stddev_initial_divisions,"stddev_initial_divisions"); // 386  if != 0, then apply it
   compare(w1.stddev_DND,w2.stddev_DND,"stddev_DND");               // 387  if != 0, then apply it



   compare(w1.pMHC_dependent_division,w2.pMHC_dependent_division,"pMHC_dependent_division");
   compare(w1.signal_dependent_number_of_divisions,w2.signal_dependent_number_of_divisions,"signal_dependent_number_of_divisions");
   compare(w1.pMHC_dependent_P_max,w2.pMHC_dependent_P_max,"pMHC_dependent_P_max");
   compare(w1.pMHC_dependent_K,w2.pMHC_dependent_K,"pMHC_dependent_K");
   compare(w1.pMHC_dependent_nHill,w2.pMHC_dependent_nHill,"pMHC_dependent_nHill");
   compare(w1.pMHC_dependent_P_min,w2.pMHC_dependent_P_min,"pMHC_dependent_P_min");
   compare(w1.pMHC_dependent_P_standard,w2.pMHC_dependent_P_standard,"pMHC_dependent_P_standard");
   compare(w1.pMHC_dependent_pMHC_of_2divisions,w2.pMHC_dependent_pMHC_of_2divisions,"pMHC_dependent_pMHC_of_2divisions");
   compare(w1.TFHsignal_dependent_K,w2.TFHsignal_dependent_K,"TFHsignal_dependent_K");
   compare(w1.TFHsignal_of_P0divisions,w2.TFHsignal_of_P0divisions,"TFHsignal_of_P0divisions");


   ///§§§ Philippe 10-04-2017
   compare(w1.time_tc_selection_block,w2.time_tc_selection_block,"time_tc_selection_block");  // (hours) -> no blocking    382
   compare(w1.factor_tc_selection_block,w2.factor_tc_selection_block,"factor_tc_selection_block");    //  -> no blocking           383
   compare(w1.time_DND_block,w2.time_DND_block,"time_DND_block");               //                           384
   compare(w1.factor_DND_block,w2.factor_DND_block,"factor_DND_block");             //  -> no blocking           385
   compare(w1.factor_founder_div_block, w2.factor_founder_div_block, "factor_founder_div_block"); // 388


   compare(w1.ICOSL_dependent_Tfh_signals,w2.ICOSL_dependent_Tfh_signals,"ICOSL_dependent_Tfh_signals");
   compare(w1.ICOSL_memory,w2.ICOSL_memory,"ICOSL_memory");
   compare(w1.ICOSL_upregulation_mode,w2.ICOSL_upregulation_mode,"ICOSL_upregulation_mode");
   compare(w1.ICOSL_upregulation_time,w2.ICOSL_upregulation_time,"ICOSL_upregulation_time");

   // TC:
   // ===============
   compare(w1.totalTC,w2.totalTC,"totalTC");
   compare(w1.TC_radius,w2.TC_radius,"TC_radius");
   compare(w1.v_TC,w2.v_TC,"v_TC");
   compare(w1.v_TC_width,w2.v_TC_width,"v_TC_width");       // defines a width of Gauss distributed v_TC values (-1 for fixed)
   compare(w1.v_TC_CC,w2.v_TC_CC,"v_TC_CC");          // TC-velocity if encountering a CC
   compare(w1.TC_persistence,w2.TC_persistence,"TC_persistence");   // average time gap in minutes between changes of direction of the
   compare(w1.north_weight,w2.north_weight,"north_weight");     // tendency "a" to walk north: p=(1-a)r+an with r: random, n: north
   compare(w1.TC_CC_selection,w2.TC_CC_selection,"TC_CC_selection");    // 1 if TC rescue CC from apoptosis
   compare(w1.TFH_CC_selection_mode,w2.TFH_CC_selection_mode,"TFH_CC_selection_mode");
   compare(w1.TFH_CC_selection_gauss_width,w2.TFH_CC_selection_gauss_width,"TFH_CC_selection_gauss_width");
   compare(w1.TFH_ASpos,w2.TFH_ASpos,"TFH_ASpos");
   compare(w1.do_TC_division,w2.do_TC_division,"do_TC_division");
   compare(w1.TC_doubling,w2.TC_doubling,"TC_doubling");
   compare(w1.TC_meancycle,w2.TC_meancycle,"TC_meancycle");
   compare(w1.TC_cyclewidth,w2.TC_cyclewidth,"TC_cyclewidth");
   compare(w1.dx_TC,w2.dx_TC,"dx_TC");
   compare(w1.TC_Ndivisions,w2.TC_Ndivisions,"TC_Ndivisions");

   // OUT:
   // ===============
   compare(w1.mk_ab,w2.mk_ab,"mk_ab");
   compare(w1.pm_differentiation_time,w2.pm_differentiation_time,"pm_differentiation_time");
   compare(w1.v_OUT,w2.v_OUT,"v_OUT");             // OUT-velocity
   compare(w1.v_OUT_width,w2.v_OUT_width,"v_OUT_width");       // defines a width of Gauss distributed v_OUT values (-1 for fixed)
   compare(w1.OUT_persistence,w2.OUT_persistence,"OUT_persistence");   // average time gap in minutes between changes of the polarity

   // Antibodies
   // =============
   compare(w1.ic_k_on,w2.ic_k_on,"ic_k_on");
   compare(w1.ic_k_off,w2.ic_k_off,"ic_k_off");
   compare(w1.ag_threshold,w2.ag_threshold,"ag_threshold");
   compare(w1.antibodies_resolution,w2.antibodies_resolution,"antibodies_resolution");
   compare(w1.antibodies_production,w2.antibodies_production,"antibodies_production");
   compare(w1.antibodies_degradation,w2.antibodies_degradation,"antibodies_degradation");
   compare(w1.k_ic_exp_min,w2.k_ic_exp_min,"k_ic_exp_min");
   compare(w1.k_ic_exp_max,w2.k_ic_exp_max,"k_ic_exp_max");
   compare(w1.N_GC,w2.N_GC,"N_GC");
   compare(w1.V_blood,w2.V_blood,"V_blood");
   compare(w1.inject_antibody,w2.inject_antibody,"inject_antibody");
   compare(w1.injected_antibody_affinity,w2.injected_antibody_affinity,"injected_antibody_affinity");
   compare(w1.inject_antibody_time,w2.inject_antibody_time,"inject_antibody_time");
   compare(w1.inject_antibody_ASindex,w2.inject_antibody_ASindex,"inject_antibody_ASindex");

   // Photoactivation
   // ================
   compare(w1.photoactivation,w2.photoactivation,"photoactivation");
   compare(w1.photoactivation_t0,w2.photoactivation_t0,"photoactivation_t0");
   compare(w1.photoactivation_x0,w2.photoactivation_x0,"photoactivation_x0");
   compare(w1.photoactivation_y0,w2.photoactivation_y0,"photoactivation_y0");
   compare(w1.photoactivation_z0,w2.photoactivation_z0,"photoactivation_z0");
   compare(w1.photoactivation_delta_x,w2.photoactivation_delta_x,"photoactivation_delta_x");
   compare(w1.photoactivation_delta_y,w2.photoactivation_delta_y,"photoactivation_delta_y");
   compare(w1.photoactivation_delta_z,w2.photoactivation_delta_z,"photoactivation_delta_z");
   compare(w1.def_DEC205,w2.def_DEC205,"def_DEC205");
   compare(w1.inject_antiDEC205OVA,w2.inject_antiDEC205OVA,"inject_antiDEC205OVA");
   compare(w1.DEC205_induce_CBdifferentiation,w2.DEC205_induce_CBdifferentiation,"DEC205_induce_CBdifferentiation");
   compare(w1.retain_DEC205_ag,w2.retain_DEC205_ag,"retain_DEC205_ag");
   compare(w1.def_DEC205_t0,w2.def_DEC205_t0,"def_DEC205_t0");
   compare(w1.inject_antiDEC205OVA_t0,w2.inject_antiDEC205OVA_t0,"inject_antiDEC205OVA_t0");
   compare(w1.antiDEC205OVA_tend,w2.antiDEC205OVA_tend,"antiDEC205OVA_tend");
   compare(w1.p_DEC205,w2.p_DEC205,"p_DEC205");
   compare(w1.TC_dec205ova_time,w2.TC_dec205ova_time,"TC_dec205ova_time");
   compare(w1.TC_factor_dec205ova,w2.TC_factor_dec205ova,"TC_factor_dec205ova");
   compare(w1.DEC205_p_factor,w2.DEC205_p_factor,"DEC205_p_factor");
   compare(w1.DEC205_forces_output,w2.DEC205_forces_output,"DEC205_forces_output");

   // FDC:
   // ===============
   compare(w1.FDCnumber,w2.FDCnumber,"FDCnumber");
   compare(w1.FDCnetwork,w2.FDCnetwork,"FDCnetwork");
   compare(w1.posFDC,w2.posFDC,"posFDC");
   compare(w1.FDClength,w2.FDClength,"FDClength");
   compare(w1.FDCtransparent,w2.FDCtransparent,"FDCtransparent");
   compare(w1.FDCvesicle,w2.FDCvesicle,"FDCvesicle");
   compare(w1.mksignal,w2.mksignal,"mksignal");
   compare(w1.mkCXCL12,w2.mkCXCL12,"mkCXCL12");
   compare(w1.mkCXCL13,w2.mkCXCL13,"mkCXCL13");
   compare(w1.mk_SEMA4D,w2.mk_SEMA4D,"mk_SEMA4D");
   compare(w1.ag_per_FDC,w2.ag_per_FDC,"ag_per_FDC");
   compare(w1.ag_saturation_FDC,w2.ag_saturation_FDC,"ag_saturation_FDC");
   compare(w1.ag_distribution_mode,w2.ag_distribution_mode,"ag_distribution_mode");
   compare(w1.ag_detection_mode,w2.ag_detection_mode,"ag_detection_mode");

   // BETA-cells:
   // =======================
   compare(w1.BETA_Nini,w2.BETA_Nini,"BETA_Nini");                      // Total Number of initial betacells:
   compare(w1.BETA_pos,w2.BETA_pos,"BETA_pos");                         // position in space
   compare(w1.BETA_v_modi,w2.BETA_v_modi,"BETA_v_modi");                // modus of velocity state treatment
   compare(w1.BETA_n_v_states,w2.BETA_n_v_states,"BETA_n_v_states");    // # of velocity states
   compare(w1.BETA_radius,w2.BETA_radius,"BETA_radius");
   compare(w1.BETA_proliferate,w2.BETA_proliferate,"BETA_proliferate"); // Rate per hr
   compare(w1.BETA_max_pro,w2.BETA_max_pro,"BETA_max_pro");             // maximal distance for CB-proliferation from dividing cell
   compare(w1.BETA_grow,w2.BETA_grow,"BETA_grow");
   compare(w1.BETA_shrink,w2.BETA_shrink,"BETA_shrink");
   compare(w1.BETA_max_adhesion,w2.BETA_max_adhesion,"BETA_max_adhesion");  // maximum adhesion force in % of full stickness
   compare(w1.BETA_persistence,w2.BETA_persistence,"BETA_persistence");   // average time gap in minutes between changes of direction of the
   compare(w1.BETA_v,w2.BETA_v,"BETA_v");                               // velocity
   compare(w1.BETA_v_factor,w2.BETA_v_factor,"BETA_v_factor");      // for 2 velocities: the factor by which the velocity is reduced
   compare(w1.BETA_v_switch_deltat,w2.BETA_v_switch_deltat,"BETA_v_switch_deltat"); // Mean duration in a v-state in minutes
   compare(w1.BETA_v_cytosol,w2.BETA_v_cytosol,"BETA_v_cytosol");     // Strength of reshaping forces (cytosolic elements speed)
   compare(w1.BETA_elongation,w2.BETA_elongation,"BETA_elongation");    // Cell elongation by active movement
   compare(w1.BETA_K_elongation,w2.BETA_K_elongation,"BETA_K_elongation");  // Elongation in units of spherical cell radius, at which
   compare(w1.BETA_distance_tolerance,w2.BETA_distance_tolerance,"BETA_distance_tolerance");
   compare(w1.BETA_half_tolerance_deformation,w2.BETA_half_tolerance_deformation,"BETA_half_tolerance_deformation");
   compare(w1.BETA_smoothmove,w2.BETA_smoothmove,"BETA_smoothmove");   // Distributes a barycenter movement thought to overcome one
   compare(w1.v_resolution,w2.v_resolution,"v_resolution");
   compare(w1.s_resolution,w2.s_resolution,"s_resolution");
   compare(w1.alpha_resolution,w2.alpha_resolution,"alpha_resolution");
   compare(w1.delta_v,w2.delta_v,"delta_v");
   compare(w1.delta_s,w2.delta_s,"delta_s");
   compare(w1.delta_alpha,w2.delta_alpha,"delta_alpha");
   compare(w1.trackfrom,w2.trackfrom,"trackfrom");
   compare(w1.trackuntil,w2.trackuntil,"trackuntil");
   compare(w1.track_delta_t,w2.track_delta_t,"track_delta_t");
   compare(w1.tALL,w2.tALL,"tALL");
   compare(w1.tCB,w2.tCB,"tCB");
   compare(w1.tCC,w2.tCC,"tCC");
   compare(w1.tOUT,w2.tOUT,"tOUT");
   compare(w1.tTC,w2.tTC,"tTC");
   compare(w1.tBETA,w2.tBETA,"tBETA");
}
