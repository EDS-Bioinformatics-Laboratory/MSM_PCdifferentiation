#include "parameters.h"
#include <cmath>
#include <sstream>

using namespace std;

void parameters::init_param() {
  // Resizes the vector for paramteres values
  par.resize(N_par, 0.);
  // Resizes the vector for paramteres names
  names.resize(N_par);
}

parameters::parameters(int takeHyphasmaFile, string parfname) {
  // This constructor reads paramteres in Hyphasma format
  // Initialize parameters fileds
  init_param();
  cout << "Reading parameter file " << parfname
       << " in 'Hyphasma' format:" << endl;
  hyphasmaParameter hpar;

  // Read parameter in Htphasma format
  bool done = hpar.wahl(parfname.c_str(), true, true);
  if (!done) {
    cerr << "Error: Could not load parameters from " << parfname << endl;
  } else {
    cerr << "Parameter reading is done." << endl;
  }
  // Transfers the parameters from hyphasma format into Mafalda format
  matchFromHyphasma(hpar);
  parameter_file_name = parfname;
}

parameters::~parameters() {}

parameters::parameters() {
  // This constructor creates a verison of Built-in parameter values

  // Initialize parameters fileds "names" and "par"
  init_param();

  // A
  names[AgAmountperFDC] = "    Presented Antigen per FDC    ";
  par[AgAmountperFDC] = 3000;
  names[agSaturation] =
      "    Ag saturation per FDC fragment in units of threshold. 1:constant "
      "finding proability    ";
  par[agSaturation] = 20;
  names[Ag_threshold] =
      "    Threshold Ag-concentration for binding CC (in Mol):    ";
  par[Ag_threshold] = 1.e-8;

  // B
  names[Bcell_speed] = "    B-Cell Speed (um / hr)    ";
  par[Bcell_speed] = 7.5;
  names[Bcell_stddev] = "    deviation of B-cell speed (um/hr.)    ";
  par[Bcell_stddev] = -1;
  names[Bcell_tp] = "    B-Cell Persistent Time average (hr.)    ";
  par[Bcell_tp] = 1.5;
  names[Bcell_tp_stddev] = "    B-Cell Persistent Time stddev (hr.)    ";
  par[Bcell_tp_stddev] = 0;
  names[BCR_pool] = "    Size of initial B-cell receptor pool    ";
  par[BCR_pool] = 100;
  names[bcr]= "Bcr signal intensity";
  par[bcr]= 1;
  names[BLIMP1th]= "BLIMP1 threshlod for plasmacell differentiation";
  par[BLIMP1th]= 8;

  // C
  names[c_G1] = "    Phase g1 of cell cycle (hr.)    ";
  par[c_G1] = 2.5;
  names[c_S] = "    Phase S of cell cycle  (hr.)    ";
  par[c_S] = 1.5;
  names[c_G2] = "    Phase g2 of cell cycle  (hr.)    ";
  par[c_G2] = 2.5;
  names[c_M] = "    Phase M of cell cycle  (hr.)    ";
  par[c_M] = 0.5;
  names[c_G1_stddev] = "    Phase g1 of cell cycle stddev  (hr.)    ";
  par[c_G1_stddev] = 1;
  names[c_S_stddev] = "    Phase S of cell cycle stddev  (hr.)    ";
  par[c_S_stddev] = 1;
  names[c_G2_stddev] = "    Phase g2 of cell cycle stddev  (hr.)    ";
  par[c_G2_stddev] = 1;
  names[c_M_stddev] = "    Phase M of cell cycle stddev  (hr.)    ";
  par[c_M_stddev] = 1;
  names[Ccdif_delay_stddev] =
      "    Standard deviation for delay to differentiation.    ";
  par[Ccdif_delay_stddev] = 0;
  names[chemo_dx] = "    Lattice Chemokine Constant (um)    ";
  par[chemo_dx] = 5;
  names[CXCL12crit] =
      "    Critical CXCL12 concentration for desensitization (mol)    ";
  par[CXCL12crit] = 0.000000006;
  names[CXCL13crit] =
      "    Critical CXCL13 concentration for desensitization (mol) //(-1 for "
      "none)?????    ";
  par[CXCL13crit] = 8.e-11;
  names[CXCL12recrit] =
      "    Critical CXCL12 concentration for resensitization (mol)    ";
  par[CXCL12recrit] = 0.000000004;
  names[CXCL13recrit] =
      "    Critical CXCL13 concentration for resensitization (mol) //(-1 for "
      "none)????? ";
  par[CXCL13recrit] = 6.e-11;
  names[chemmax] = "    Maximum weigh of chemotaxis    ";
  par[chemmax] = 10;
  names[chemosteep] =
      "    Steepness of weight reduction with chemokine gradient (mol/l)    ";
  par[chemosteep] = 1.e+10;
  names[chemohalf] = "    Chemokine gradient of half weight (l/mol)    ";
  par[chemohalf] = 2.e-11;
  names[collectionFDCperiod] =
      "    Duration of CC collection of Antigen by serial encounters with FDC "
      "(hr.)    ";
  par[collectionFDCperiod] = 0.7;
  names[CB_radius] = "    Centroblast radius (um)    ";
  par[CB_radius] = 2.45;
  names[cd40]= "Cd40 signal intensity";
  par[cd40]= 50;

  // D
  names[DendriteLength] =
      "    Length FDC dendrites / dx (number of positions)    ";
  par[DendriteLength] = 40;
  names[dimension] = "    Lattice Dimensions    ";
  par[dimension] = 3;
  names[dt] = "    Time resolution (hr)    ";
  par[dt] = 0.002;
  names[dx] = "    Lattice Constant (um)    ";
  par[dx] = 5;
  names[DeleteAgInFreshCC] = "    Retained Ag is deleted in fresh CC    ";
  par[DeleteAgInFreshCC] = true;
  names[difDelay] =
      "    Delay cell differentiation after TC selection (hr.)    ";
  par[difDelay] = 6;

  // E
  names[expMin] = "Conversion of shape space affinity to (1/mol)";
  par[expMin] = 5.5;
  names[expMax] = "Conversion of shape space affinity to (1/mol)";
  par[expMax] = 9.5;
  names[eta] = "    Exponent of the hamming distance    ";
  par[eta] = 2;

  // G
  names[Gamma] = "    Width of gaussian affinity weight function    ";
  par[Gamma] = 2.8;

  // I
  names[InitialNumberSC] = "    Initial Number Stromal cells    ";
  par[InitialNumberSC] = 300;
  names[InitialNumberTC] = "    Initial Number T-cells    ";
  par[InitialNumberTC] = 250;
  names[InitialNumberCB] = "    Initial Number Centroblasts    ";
  par[InitialNumberCB] = 0;
  names[InitialNumberFDC] = "    Initial Number FDCs    ";
  par[InitialNumberFDC] = 200;

  // K
  names[kon] = "k_on for building immune complex (1/mol h)    ";
  par[kon] = 1.e6;

  names[koff] = "k_off for dissociation of immune complex (in /s):     ";
  par[koff] = 0.001;

  // L
  names[BCR_Length] = "    Length of BCRs    ";
  par[BCR_Length] = 4;

  // M
  names[macrophage] = "Rate of macrophage transport of dead cells (h):";
  par[macrophage] = 6.0;

  //    names[Memorycell_tp]     = "    Memory Cell speed (um / hr.)    ";
  //    par[Memorycell_tp]       =   0.0125      ;

  //    names[Memorycell_speed]  = "    Memory Cell polarity (degrees)    ";
  //    par[Memorycell_speed]    =   0.05        ;

  //    names[Memorycell_tp_stddev] = "    Memory Cell polarity (degrees)    ";
  //    par[Memorycell_tp_stddev] = 0           ;

  // N
  names[Nmax] = "    Maximum number of residues in one dimension    ";
  par[Nmax] = 9;
  names[NoMutFounderCells] = "    FounderCellsDoNotMutate    ";
  par[NoMutFounderCells] = false;
  names[nDiv] = "    Number of divisions of founder cells ";
  par[nDiv] = 12;
  names[nDiv_stddev] = "    stddev of Number of divisions  of founder cells ";
  par[nDiv_stddev] = 0;
  names[nDivinflow] = "    Number of divisions of influx Bcells    ";
  par[nDivinflow] = 6;
  Avogadro_constant = 6.02205e+23;  // mol^-1, Avogadro number

  // P
  names[Plasmacell_tp] = "    Plasma Cell persistence time (unit)   ";
  par[Plasmacell_tp] = 0.75;
  names[Plasmacell_speed] = "    Plasma Cell speed   (unit) ";
  par[Plasmacell_speed] = 3.0;
  names[Plasmacell_tp_stddev] = "    Plasma Cell polarity (degrees)    ";
  par[Plasmacell_tp_stddev] = -1;
  names[pMHCdepHill] =
      "    p-MHC dependent division number Hill (Hill coef. n_P)    ";
  par[pMHCdepHill] = 2;
  names[pMHCdepMin] =
      "    p-MHC dependent division number Hill (Hill coef. P_Min)    ";
  par[pMHCdepMin] = 1;
  names[pMHCdepMax] =
      "    p-MHC dependent division number Hill (Hill coef. P_Max)    ";
  par[pMHCdepMax] = 6;
  names[pMHCdepK] =
      "    p-MHC dependent division number Hill (Hill coef. K_P)    ";
  par[pMHCdepK] = 9;
  names[pmutB4StartMut] =
      "    Probability of mutation before first 24 hours     ";
  par[pmutB4StartMut] = 0;
  names[pmutAfterStartMut] =
      "    Probability of mutation after first 24 hours    ";
  par[pmutAfterStartMut] = 0.5;  // 0.5
  names[pmutAfterSelection] =
      "    Probability of mutation after selection (affinity dependant)   ";
  par[pmutAfterSelection] = 0;
  names[pmutAffinityExponent] =
      "    Affinity dependant mutation upon TC contact (affinity-exponent)    ";
  par[pmutAffinityExponent] = 1;
  names[pDivideAgAssymetric] =
      "    Probability to divide Ag assymetrically to daughter B-cell    ";
  par[pDivideAgAssymetric] = 0.72;
  names[polarityIndex] = "    Assymetric Distribution of Ag    ";
  par[polarityIndex] = 1.0;
  names[pSel] = "     Probability to be selected by FDC.    ";
  par[pSel] =
      20;  // -->original  (1/0.05)        ;  //#Rechcek, need to fix this
  names[pApoCC] = "    % Casp3+ LZ cells per hr. used as (apoptosis rate)     ";
  par[pApoCC] = 0;
  names[pApoCB] = "    % Casp3+ DZ cells per hr. used as (apoptosis rate)     ";
  par[pApoCB] = 0;
  names[p_dif] = "    Differentiation rate    ";
//Elena: Ask Danial Why noparameter value for p_dif? (see conversion)
  names[polarityBLIMP1]="Polarity of BLIMP1 in B-cells that divide Ag Asymmetrically [0-1]";
  par[polarityBLIMP1]= 0.5;
  names[polarityBCL6]="Polarity of BCL6 in B-cells that divide Ag Asymmetrically [0-1]";
  par[polarityBCL6]=0.5;
  names[polarityIRF4]="Polarity of IRF4 in B-cells that divide Ag Asymmetrically [0-1]";
  par[polarityIRF4]=0.5;

  // R
  names[radius] = "    Lattice Radius (um)    ";
  par[radius] = 160;
  names[rateCBinflow] = "    rate of inflow (cells/hr.)    ";
  par[rateCBinflow] = 2.0;

  // S
  names[smoothnessStopCBinflow] =
      "    smoothness of stop inflow CB (hr.) (-1 = no)    ";
  par[smoothnessStopCBinflow] = 6;
  names[StartMutation] = "    Start of mutation period  (hr.)    ";
  par[StartMutation] = 24;

  // T
  names[timeStopCBinflow] = "    time to stop inflow CB  (hr.)    ";
  par[timeStopCBinflow] = 96.0;
  names[Tcell_speed] = "    T-Cell Speed (um / hr.)    ";
  par[Tcell_speed] = 10;
  names[Tcell_tp] = "    T-Cell Persistent Time average (hr.)    ";
  par[Tcell_tp] = 1.7;  // hour
  names[Tcell_stddev] = "    deviation of T-cell speed (um/sec)    ";
  par[Tcell_stddev] = -1;
  names[Tcell_tp_stddev] = "    T-Cell Persistent Time stddev (hr.)    ";
  par[Tcell_tp_stddev] = 0;
  names[tcTime] = "    Duration of CC-Tc contact  (hr.)    ";
  par[tcTime] = 0.6;  // hour
  names[tcRescueTime] =
      "    Minimum duration of TC-CC-polarization for CC-rescue  (hr.)    ";
  par[tcRescueTime] = 0.5;  // hour 0.5
  names[testDelay] = "    Time gap between TFHC-CC binding tests  (hr.)    ";
  par[testDelay] = 0.02;  // Unit?
  names[tmax] = "    Maximum duration of GC simulation  (hr.)    ";
  par[tmax] = 504;
  names[tolight] = " Rate for differentiation of centroblasts to centrocytes ";
  par[tolight] = 0.1;  // in hours, #temporary, check hyphasma
  names[typeCD40signal]="Decisionn to run Affinity (1) or Fixed (0) CD40 signal models";
  par[typeCD40signal]=1;

  // W
  par[widthPI] = 0.04;
  names[widthPI] = "    Coefficient of variation arround Polarity Index    ";

  // Z
  par[zoneRatioGC] = 0.5;
  names[zoneRatioGC] =
      "    Ratio that determines the position of DZ in Germinal Center    ";
}

void parameters::writeparameters(string fname) {
  ofstream myfile;
  myfile.open(fname);
  if (!myfile) cerr << "ERROR! My parameter file empty";
  myfile << print();
  myfile.close();
}

//Elena: Reads a file that has first line with parameter name second with value. For future extension.
bool parameters::readparameters(string fname) {
  // danial: This needs to be fixed, in different versions there are different
  // types of reading parameters. This is not probabily compatible with current
  // code.
    //Elena: this is only necesary to read from a hyphasma file
  //  ifstream myfile;
  //  myfile.open(fname);
  //  if (myfile.is_open()) {
  //    string nameLine;  // parameter mame in imput file
  //    int count = 0;    // Counter for while loop
  //    while ((count < 1000000) && (getline(myfile, nameLine))) {
  //      double value = NAN;  // If value not found in names, parameter value
  //      will
  //                           // have NAN. Create an error mesage!
  //      myfile >> value;     // Open myfil and take value
  //      bool found = false;  // Check if nameLine (parameter of input file)
  //      was
  //                           // found twice or not found
  //      for (unsigned int i = 0; i < names.size(); i++) {
  //        if (nameLine.compare(names[i]) == 0)  // Compare name with nameline;
  //                                              // returns 0 if the strings
  //                                              are
  //                                              // identical.
  //        {
  //          // cerr << " reads " << nameLine << " value " << value << endl;
  //          par[i] = value;
  //          if (found) cerr << nameLine << "found 2 times or more" << endl;
  //          found = true;
  //        }
  //      }
  //      if (!found) cerr << nameLine << " not found" << endl;
  //      getline(myfile, nameLine);  // to finish line
  //    }
  //
  //    for (unsigned int i = 0; i < par.size(); i++) {
  //      if (isnan(par[i])) {
  //        cerr << "ERROR! Parameter missing in imput file" << i
  //             << " name:" << names[i] << endl;
  //      }
  //    }
  //    myfile.close();
  //    if (count >= 1000000) {
  //      cerr << "infinite loop" << endl;
  //    } else {
  //      return true;
  //    }
  //
  //  } else
  //    cerr << "ERROR! Unable to open file " << fname << endl;
  return false;
}

string parameters::print() {
  // This function returns a string that includes all parameters with their
  // values.
  stringstream res;
  for (int i = 0; i < (N_par - 1); i++) {
    res << names[i] << endl;
    res << par[i] << endl;
  }
  return res.str();
}

void parameters::matchFromHyphasma(hyphasmaParameter &hypar) {
  /**
   

      - Parameter : This function is written to put paramters from hyphasma
   parameter file to built-in parameter fields
   **/

  // A
  names[AgAmountperFDC] = "    Presented Antigen per FDC    ";
  par[AgAmountperFDC] = hypar.Value.ag_per_FDC;

  names[agSaturation] =
      "    Ag saturation per FDC fragment in units of threshold. 1:constant "
      "finding proability    ";
  par[agSaturation] = hypar.Value.ag_saturation_FDC;

  names[Ag_threshold] =
      "    Threshold Ag-concentration for binding CC (in Mol):    ";
  par[Ag_threshold] = hypar.Value.ag_threshold;

  // B
  names[Bcell_speed] = "    B-Cell Speed (um / hr)    ";
  par[Bcell_speed] = hypar.Value.v_CB;

  names[Bcell_stddev] = "    deviation of B-cell speed (um/hr.)    ";
  par[Bcell_stddev] = hypar.Value.v_CB_width;

  names[Bcell_tp] = "    B-Cell Persistent Time average (hr.)    ";
  par[Bcell_tp] = hypar.Value.CB_persistence;

  names[Bcell_tp_stddev] = "    B-Cell Persistent Time stddev (hr.)    ";
  par[Bcell_tp_stddev] = 0;

  names[BCR_pool] = "    Size of initial B-cell receptor pool    ";
  par[BCR_pool] = hypar.Value.totalBss;
  
  names[bcr]= "Bcr signal intensity";
    par[bcr]= hypar.Value.bcr;
    
  names[BLIMP1th]= "BLIMP1 threshlod for plasmacell differentiation";
    par[BLIMP1th]= hypar.Value.BLIMP1th;


  // C
  names[c_G1] = "    Phase g1 of cell cycle (hr.)    ";
  par[c_G1] = hypar.Value.CB_dt_G1;

  names[c_S] = "    Phase S of cell cycle  (hr.)    ";
  par[c_S] = hypar.Value.CB_dt_S;

  names[c_G2] = "    Phase g2 of cell cycle  (hr.)    ";
  par[c_G2] = hypar.Value.CB_dt_G2;

  names[c_M] = "    Phase M of cell cycle  (hr.)    ";
  par[c_M] = hypar.Value.CB_dt_M;

  names[c_G1_stddev] = "    Phase g1 of cell cycle stddev  (hr.)    ";
  par[c_G1_stddev] = hypar.Value.CB_dtphase_width;

  names[c_S_stddev] = "    Phase S of cell cycle stddev  (hr.)    ";
  par[c_S_stddev] = hypar.Value.CB_dtphase_width;

  names[c_G2_stddev] = "    Phase g2 of cell cycle stddev  (hr.)    ";
  par[c_G2_stddev] = hypar.Value.CB_dtphase_width;

  names[c_M_stddev] = "    Phase M of cell cycle stddev  (hr.)    ";
  par[c_M_stddev] = hypar.Value.CB_dtphase_width;

  names[Ccdif_delay_stddev] =
      "    Standard deviation for delay to differentiation.    ";
  par[Ccdif_delay_stddev] = 0;

  names[chemo_dx] = "    Lattice Chemokine Constant (um)    ";
  par[chemo_dx] = hypar.Value.dx_signal;

  names[CXCL12crit] =
      "    Critical CXCL12 concentration for desensitization (mol)    ";
  par[CXCL12crit] = hypar.Value.CXCL12crit;

  names[CXCL13crit] =
      "    Critical CXCL13 concentration for desensitization (mol) //(-1 for "
      "none)?????    ";
  par[CXCL13crit] = hypar.Value.CXCL13crit;

  names[CXCL12recrit] =
      "    Critical CXCL12 concentration for resensitization (mol)    ";
  par[CXCL12recrit] = hypar.Value.CXCL12recrit;

  names[CXCL13recrit] =
      "    Critical CXCL13 concentration for resensitization (mol) //(-1 for "
      "none)????? ";
  par[CXCL13recrit] = hypar.Value.CXCL13recrit;

  names[chemmax] = "    Maximum weigh of chemotaxis    ";
  par[chemmax] = hypar.Value.chemo_max;

  names[chemosteep] =
      "    Steepness of weight reduction with chemokine gradient (mol/l)    ";
  par[chemosteep] = hypar.Value.chemo_steep;

  names[chemohalf] = "    Chemokine gradient of half weight (l/mol)    ";
  par[chemohalf] = hypar.Value.chemo_half;

  names[collectionFDCperiod] =
      "    Duration of CC collection of Antigen by serial encounters with FDC  "
      "(hr.)    ";
  par[collectionFDCperiod] = hypar.Value.collectFDCperiod;

  names[CB_radius] = "    Centroblast radius (um)    ";
  par[CB_radius] = hypar.Value.CB_radius;
    
  names[cd40]= "Cd40 signal intensity";
    par[cd40]= hypar.Value.cd40;
  
  // D
  names[DendriteLength] =
      "    Length FDC dendrites / dx (number of positions)    ";
  par[DendriteLength] = hypar.Value.FDClength;

  names[dimension] = "    Lattice Dimensions    ";
  par[dimension] = hypar.Value.DimSpace;

  names[dt] = "    Time resolution (hr)    ";
  par[dt] = hypar.Value.deltat;

  names[dx] = "    Lattice Constant (um)    ";
  par[dx] = hypar.Value.dx;

  names[DeleteAgInFreshCC] = "    Retained Ag is deleted in fresh CC    ";
  par[DeleteAgInFreshCC] = hypar.Value.ag_loaded_CB_diff2output;

  names[difDelay] =
      "    Delay cell differentiation after TC selection (hr.)    ";
  par[difDelay] = hypar.Value.ccdiff_delay;

  // E
  names[expMin] = "Conversion of shape space affinity to (1/mol)";
  par[expMin] = hypar.Value.k_ic_exp_min;

  names[expMax] = "Conversion of shape space affinity to (1/mol)";
  par[expMax] = hypar.Value.k_ic_exp_max;

  names[eta] = "    Exponent of the hamming distance    ";
  par[eta] = hypar.Value.amplitudeGauss;

  // G
  names[Gamma] = "    Width of gaussian affinity weight function    ";
  par[Gamma] = hypar.Value.GammaGauss;

  // I
  names[InitialNumberSC] = "    Initial Number Stromal cells    ";
  par[InitialNumberSC] = 300;

  names[InitialNumberTC] = "    Initial Number T-cells    ";
  par[InitialNumberTC] = hypar.Value.totalTC;

  names[InitialNumberCB] = "    Initial Number Centroblasts    ";
  par[InitialNumberCB] =
      hypar.Value.totalB;  // Danial:changed from totalBss to total B

  names[InitialNumberFDC] = "    Initial Number FDCs    ";
  par[InitialNumberFDC] = hypar.Value.FDCnumber;

  // K
  names[kon] = "k_on for building immune complex (1/mol h)    ";
  par[kon] = hypar.Value.ic_k_on;

  names[koff] = "k_off for dissociation of immune complex (in /s):     ";
  par[koff] = hypar.Value.ic_k_off;

  // L
  names[BCR_Length] = "    Length of BCRs    ";
  par[BCR_Length] = hypar.Value.DimShapeSpace;

  // M
  names[macrophage] = "Rate of macrophage transport of dead cells (h):";
  par[macrophage] = hypar.Value.macrophage;

  //    names[Memorycell_tp]     = "    Memory Cell speed (um / hr.)    ";
  //    par[    Memorycell_tp     ] = hypar.Value.    OUT_persistence    ;

  //    names[Memorycell_speed]  = "    Memory Cell polarity (degrees)    ";
  //    par[    Memorycell_speed    ] = hypar.Value.    v_OUT    ;

  //    names[Memorycell_tp_stddev] = "    Memory Cell polarity (degrees)    ";
  //    par[    Memorycell_tp_stddev     ] = hypar.Value.    v_OUT_width    ;

  // N
  names[Nmax] = "    Maximum number of residues in one dimension    ";
  par[Nmax] = 9;

  names[NoMutFounderCells] = "    FounderCellsDoNotMutate    ";
  par[NoMutFounderCells] = false;

  names[nDiv] =
      "    Number of divisions of founder cells "
      "((StartDifferentiation-tmin)/cell cycle duration)    ";
  par[nDiv] = hypar.Value.CB_fixed_times_of_divisions_in_expansion;

  names[nDiv_stddev] =
      "    stddev of Number of divisions  of founder cells "
      "((StartDifferentiation-tmin)/cell cycle duration)     ";
  par[nDiv_stddev] = hypar.Value.stddev_initial_divisions;

  names[nDivinflow] = "    Number of divisions of influx Bcells    ";
  par[nDivinflow] =
      6;  // Elena: Is there any same parameter for hyphasma? CHECsK****

  Avogadro_constant = 6.02205e+23;  // mol^-1, Avogadro number

  // P
  names[Plasmacell_tp] = "    Plasma Cell persistence time (unit)   ";
  par[Plasmacell_tp] = hypar.Value.OUT_persistence;

  names[Plasmacell_speed] = "    Plasma Cell speed   (unit) ";
  par[Plasmacell_speed] = hypar.Value.v_OUT;

  names[Plasmacell_tp_stddev] = "    Plasma Cell polarity (degrees)    ";
  par[Plasmacell_tp_stddev] = hypar.Value.v_OUT_width;

  names[pMHCdepHill] =
      "    p-MHC dependent division number Hill (Hill coef. n_P)    ";
  par[pMHCdepHill] = hypar.Value.pMHC_dependent_nHill;

  names[pMHCdepMin] =
      "    p-MHC dependent division number Hill (Hill coef. P_Min)    ";
  par[pMHCdepMin] = hypar.Value.pMHC_dependent_P_min;

  names[pMHCdepMax] =
      "    p-MHC dependent division number Hill (Hill coef. P_Max)    ";
  par[pMHCdepMax] = hypar.Value.pMHC_dependent_P_max;

  names[pMHCdepK] =
      "    p-MHC dependent division number Hill (Hill coef. K_P)    ";
  par[pMHCdepK] = hypar.Value.pMHC_dependent_K;

  names[pmutB4StartMut] =
      "    Probability of mutation before first 24 hours     ";
  par[pmutB4StartMut] = 0;

  names[pmutAfterStartMut] =
      "    Probability of mutation after first 24 hours    ";
  par[pmutAfterStartMut] = hypar.Value.mutation;

  names[pmutAfterSelection] =
      "    Probability of mutation after selection (affinity dependant)   ";
  par[pmutAfterSelection] = hypar.Value.mutation_after_tc;

  names[pmutAffinityExponent] =
      "    Affinity dependant mutation upon TC contact (affinity-exponent)    ";
  par[pmutAffinityExponent] = hypar.Value.mutation_affinity_exponent;

  names[pDivideAgAssymetric] =
      "    Probability to divide Ag assymetrically to daughter B-cell    ";
  par[pDivideAgAssymetric] = hypar.Value.divide_ag_asymmetric;

  names[polarityIndex] = "    Assymetric Distribution of Ag    ";
  par[polarityIndex] = hypar.Value.asymmetric_polarity_index;

  names[pSel] = "     Probability to be selected by FDC.    ";
  par[pSel] = hypar.Value.selection;

  names[pApoCC] = "    % Casp3+ LZ cells per hr. used as (apoptosis rate)     ";
  par[pApoCC] = 0;

  names[pApoCB] = "    % Casp3+ DZ cells per hr. used as (apoptosis rate)     ";
  par[pApoCB] = 0;
    
  //Elena: Ask Danial Why noparameter value for p_dif? (see conversion)
  names[p_dif] = "    Differentiation rate    ";
    
    names[polarityBLIMP1]="Polarity of BLIMP1 in B-cells that divide Ag Asymmetrically [0-1]";
    par[polarityBLIMP1]= hypar.Value.pBLIMP1;
    
    names[polarityBCL6]="Polarity of BCL6 in B-cells that divide Ag Asymmetrically [0-1]";
    par[polarityBCL6]=hypar.Value.pBCL6;
    
    names[polarityIRF4]="Polarity of IRF4 in B-cells that divide Ag Asymmetrically [0-1]";
    par[polarityIRF4]=hypar.Value.pIRF4;

  // R
  names[radius] = "    Lattice Radius (um)    ";
  par[radius] = hypar.Value.GC_radius;

  names[rateCBinflow] = "    rate of inflow (cells/hr.)    ";
  par[rateCBinflow] = hypar.Value.newBCinflux_rate;

  // S
  names[smoothnessStopCBinflow] =
      "    smoothness of stop inflow CB (hr.) (-1 = no)    ";
  par[smoothnessStopCBinflow] = hypar.Value.smooth_stopBCinflux;

  names[StartMutation] = "    Start of mutation period  (hr.)    ";
  par[StartMutation] = hypar.Value.Start_Mutation;

  // T
  names[timeStopCBinflow] = "    time to stop inflow CB  (hr.)    ";
  par[timeStopCBinflow] = hypar.Value.newBCinflux_stop;

  names[Tcell_speed] = "    T-Cell Speed (um / hr.)    ";
  par[Tcell_speed] = hypar.Value.v_TC;

  names[Tcell_tp] = "    T-Cell Persistent Time average (hr.)    ";
  par[Tcell_tp] = hypar.Value.TC_persistence;

  names[Tcell_stddev] = "    deviation of T-cell speed (um/sec)    ";
  par[Tcell_stddev] = hypar.Value.v_TC_width;

  names[Tcell_tp_stddev] = "    T-Cell Persistent Time stddev (hr.)    ";
  par[Tcell_tp_stddev] = 0;

  names[tcTime] = "    Duration of CC-Tc contact  (hr.)    ";
  par[tcTime] = hypar.Value.TC_time;

  names[tcRescueTime] =
      "    Minimum duration of TC-CC-polarization for CC-rescue  (hr.)    ";
  par[tcRescueTime] = hypar.Value.TC_rescue_time;

  names[testDelay] = "    Time gap between TFHC-CC binding tests  (hr.)    ";
  par[testDelay] = hypar.Value.CC_test_delay;

  names[tmax] = "    Maximum duration of GC simulation  (hr.)    ";
  par[tmax] = hypar.Value.tmax;

  names[tolight] = " Rate for differentiation of centroblasts to centrocytes ";
  par[tolight] = hypar.Value.tolight;
    
  names[typeCD40signal]="Decision to run Affinity (1) or Fixed (0) CD40 signal models";
  par[typeCD40signal]=hypar.Value.type_CD40_signal;

  // W
  names[widthPI] = "    Coefficient of variation arround Polarity Index    ";
  par[widthPI] = hypar.Value.smooth_PI;
    
    

  // Z
  names[zoneRatioGC] =
      "    Ratio that determines the position of DZ in Germinal Center    ";
  par[zoneRatioGC] = 0.5;
}

void parameters::convert_parameters() {
  // It should be mentioned in output files that parameters are converted based
  // on time and size steps.  #temporary
    //#Recheck
  par[kon] = par[kon] * 3600. * par[Ag_threshold];  // /(Mol hour) ?
  par[koff] = par[koff] * 3600.;                    // /(Mol hour) ?
  par[macrophage] = (log(2) / par[macrophage]) * par[dt];
  par[testDelay] = int(par[testDelay] / par[dt] + 0.5);

  // Persistence time
  // CB
  if (par[Bcell_tp] > 60. * par[dt]) {
    par[Bcell_tp] = 60. * par[dt] / par[Bcell_tp];
  } else {
    par[Bcell_tp] = 1.;
  }

  // Out
  if (par[Plasmacell_tp] > 60. * par[dt]) {
    par[Plasmacell_tp] = 60. * par[dt] / par[Plasmacell_tp];
  } else {
    par[Plasmacell_tp] = 1.;
  }

  // TC
  if (par[Tcell_tp] > 60. * par[dt]) {
    par[Tcell_tp] = 60. * par[dt] / par[Tcell_tp];
  } else {
      par[Tcell_tp] = 1.;
  }

  par[p_dif] =//Elena: 1/par[tolight] should be in linkHyphasmaMafalda! this is only conversion!
    (1 / par[tolight]) * par[dt];  // #temporary make it specific for cc and cb

  par[chemosteep] = par[chemosteep] / (par[dx] * par[dx] * par[dx] * 1.e-15 *
                                       Avogadro_constant);  // in # molecules
  par[chemohalf] = par[chemohalf] * par[dx] * par[dx] * par[dx] * 1.e-15 * Avogadro_constant;

  // in # molecules
  par[CXCL12crit] =
      par[CXCL12crit] * par[dx] * par[dx] * par[dx] * 1.e-15 * Avogadro_constant;
  // in # molecules
  par[CXCL12recrit] =
      par[CXCL12recrit] * par[dx] * par[dx] * par[dx] * 1.e-15 * Avogadro_constant;
  // in # molecules
  par[CXCL13crit] =
      par[CXCL13crit] * par[dx] * par[dx] * par[dx] * 1.e-15 * Avogadro_constant;
  // in # molecules
  par[CXCL13recrit] =
      par[CXCL13recrit] * par[dx] * par[dx] * par[dx] * 1.e-15 * Avogadro_constant;

  //////////////////////////////////////////////////////////////////////////////////

  int smoothmove = 1;  // #temproary, move to parameters file
  par[Bcell_speed] = 60. * par[Bcell_speed] * par[dt] / par[dx];
  double smoothfactor = smoothmove * par[dx] * par[dx] /
                        (12.56 * par[CB_radius] * par[CB_radius]);
  double dfactor = 0.7;  // approximate value for use_D_correction==0
  if (par[Bcell_speed] * smoothfactor / (dfactor) > 1.0) {
    cout << "Centroblast-Diffusion (p="
         << par[Bcell_speed] * smoothfactor / (dfactor)
         << ") is too large for dx and dt !!!\n";
    exit(1);
  }

  // TC speed
  par[Tcell_speed] = 60. * par[Tcell_speed] * par[dt] / par[dx];
  if (par[Tcell_speed] > 0.5) {
    cout << "TC-motility (p=" << par[Tcell_speed]
         << ") is too large for dx and dt !!!\n";
    exit(1);
  }

  // Out speed
  par[Plasmacell_speed] = 60. * par[Plasmacell_speed] * par[dt] / par[dx];
  if (par[Plasmacell_speed] > 0.5) {
    cout << "TC-motility (p=" << par[Plasmacell_speed]
         << ") is too large for dx and dt !!!\n";
    exit(1);
  }

  par[DendriteLength] = int(par[DendriteLength] / par[dx]);
  par[pSel] = par[pSel] * par[dt];  //"Rate of positive
                                    //selection at FDCs" * dt , the rate of
                                    //position selection is Psel and we save the
                                    //new value that is the result of
                                    //multiplication by dt, in pSel itself.

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
}
