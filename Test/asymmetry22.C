// Macro to extract the aymmetry from histograms.
// Scale the AMO data to the PARA and PERP histograms using flux obtained for each energy bin using asymmetry16_full_fit.C.
// Use dilution factor obtained from text files (binned dilution factor)
// Normalise the PARA and PERP histograms by dividing them by the AMO histogram.

// This macro will fit a simplified function to the asymmetry in order to extract G
// ie. assuming that flux is the same for both PARA and PERP

// Code has been simplified from analysis_v8_G: just have one bin for energy and theta
// and not comparing the asymmetry to the MAID and SAID data sets

// !!!!!!!!!!!!!!!!!!!!!Notes on improvements/corrections required:!!!!!!!!!!!!!!!!!

// Only plotting with statistical errors at the moment.
// Need to check calculation of systematic errors.
// Need to fit the cos2phisin2phi function to the normalised PARA and PERP histograms to obtain phi offset
// Do I need to include the AMO flux as a term as divide through by AMO??????????

// ------------------------------------------------------------------------

// ----------------------- Constants to be edited by user -----------------

#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TChain.h"
#include <TTree.h>
#include <Riostream.h> 
#include <TStyle.h>

float Energy_at_peak;
float Energy_cms_min;
float Energy_cms_max;

int at_Ebin = 0;
int no_Ebin =4;  // number of energy bins
double Ebin_width =0.025; // width of energy bins
double Ebin_width_v[20][20]; // width of energy bins

double min_energy; // minimum cms energy            !!!!!!!!!Change for every coherent peak setting!!!!!!!!!!!!!!!!!!

double no_Thbin =10 ; // number of cos(theta) bins
double Thbin_width=0.2; // width of cos(theta) bins


int tot_bin_energy_butanol[20][100] ;
int tot_bin_energy_ch2[20][100];
int tot_bin2_energy_butanol[20][100] ;
int tot_bin2_energy_ch2[20][100];

int tot_angle_ch2_v[20] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
int tot_angle_butanol_v[20] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
int tot_energy_ch2_v[40] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
int tot_energy_butanol_v[40] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


// Beam polarisation -------------------------------------------------------------------------------------------------------------

double PARA_pol[6] = {0.0,0.0,0.0,0.0,0.0,0.0} ;  //  Initializing variables
double PARA_pol_sys_error[6] = {0.0,0.0,0.0,0.0,0.0,0.0} ;  
double PERP_pol[6] = {0.0,0.0,0.0,0.0,0.0,0.0} ;  // 
double PERP_pol_sys_error[6] = {0.0,0.0,0.0,0.0,0.0,0.0} ; 



// double PARA_pol[4] = {0.52,0.554419,0.769091,0.863043};  // 0.730GeV PARA beam polarisation
// double PARA_pol_sys_error[4] = {0.08,.0.08,0.08,0.08};  
// double PERP_pol[4] = {0.52,0.554419,0.769091,0.863043};  // 0.730GeV PERP beam polarisation
// double PERP_pol_sys_error[4] = {0.07,0.07,0.07,0.07}; 

// double PARA_pol[4] = {0.317273,0.551818,0.740435,0.831818};  // 0.930GeV PARA beam polarisation
// double PARA_pol_sys_error[4] = {};  
// double PERP_pol[4] = {0.518182,0.7,0.813043,0.859091};  // 0.930GeV PERP beam polarisation
// double PERP_pol_sys_error[4] = {}; 

// double PARA_pol[3] = {0.64,0.80381,0.854};  // 1.1GeV PARA beam polarisation
// double PARA_pol_sys_error[3] = {};  
// double PERP_pol[3] = {0.508,0.8,0.86};  // 1.1GeV PERP beam polarisation
// double PERP_pol_sys_error[3] = {}; 

// double PARA_pol[2] = {0.651,0.800667};  // 1.3GeV PARA beam polarisation
// double PARA_pol_sys_error[2] = {};  
// double PERP_pol[2] = {0.622333,0.818};  // 1.3GeV PERP beam polarisation
// double PERP_pol_sys_error[2] = {}; 

// double PARA_pol[3] = {0.484,0.6872,0.774615};  // 1.5GeV PARA beam polarisation
// double PARA_pol_sys_error[3] = {};  
// double PERP_pol[3] = {0.3952,0.6604,0.776154};  // 1.5GeV PERP beam polarisation
// double PERP_pol_sys_error[3] = {}; 

// double PARA_pol[2] = {0.424651,0.688};  // 1.7GeV PARA beam polarisation
// double PARA_pol_sys_error[2] = {};  
// double PERP_pol[2] = {0.290698,0.662889};  // 1.7GeV PERP beam polarisation
// double PERP_pol_sys_error[2] = {}; 

// double PARA_pol[5] = {0.421739,0.5825,0.683478,0.749565,0.7625};  // 1.9GeV PARA beam polarisation
// double PARA_pol_sys_error[5] = {};  
// double PERP_pol[5] = {0.307391,0.526667,0.667826,0.751304,0.7325};  // 1.9GeV PERP beam polarisation
// double PERP_pol_sys_error[5] = {}; 

// double PARA_pol[5] = {0.398148,0.526552,0.637143,0.714815,0.748667};  // 2.1GeV PARA beam polarisation
// double PARA_pol_sys_error[5] = {};  
// double PERP_pol[5] = {0.286667,0.458276,0.61,0.706296,0.748};  // 2.1GeV PERP beam polarisation
// double PERP_pol_sys_error[5] = {}; 

// double PARA_pol[6] = {0.250385,0.363571,0.474444,0.575,0.653333,0.695};  // 2.1GeV PARA beam polarisation
// double PARA_pol_sys_error[6] = {};  
// double PERP_pol[6] = {0.256154,0.362143,0.466667,0.567857,0.645926,0.697143};  // 2.1GeV PERP beam polarisation
// double PERP_pol_sys_error[6] = {}; 


// --------------------------------------------------------------------------------------------------------------

// Flux ratios-------------------------------------------------------------------------------------------------------------------

double val_flux_ratio_CH2[10][6];  //  positive flux ratio for scaling PARA [angle][energy]
double val_flux_ratio_butanol[10][6];
double val_flux_ratio_error_CH2[10][6];  //  error positive flux ratio for scaling PARA
double val_flux_ratio_error_butanol[10][6];

double flux_ratio_CH2[6] = {0.0,0.0,0.0,0.0,0.0,0.0};  //  positive flux ratio for scaling PARA
double flux_ratio_butanol[6] = {0.0,0.0,0.0,0.0,0.0,0.0};

// double flux_ratio_CH2[4] = {1.04936,1.01183,0.84698,0.914806};  // 0.730GeV positive flux ratio for scaling PARA
// double flux_ratio_butanol[4] = {1.09537,0.994094,0.893706,0.918806};

// double flux_ratio_CH2[4] = {2.39058,1.95599,1.67399,1.69843 };  // 0.730GeV negative flux ratio for scaling PARA
// double flux_ratio_butanol[4] = {2.13554,1.91095,1.60825,1.71249 };

// double flux_ratio_CH2[4] = { 1.53722,1.28628,1.41563,1.26966};  // 0.930GeV positive flux ratio for scaling PARA
// double flux_ratio_butanol[4] = { 1.53788,1.51439,1.36269,1.29688 };

// double flux_ratio_CH2[4] = { 0.979563,1.07009,1.01965,0.829674};  // 0.930GeV negative flux ratio for scaling PARA
// double flux_ratio_butanol[4] = {0.981861,1.05758,0.99833,0.841995 };

// double flux_ratio_CH2[3] = {0.630569,0.586224,0.543227 };  // 1.1GeV positive flux ratio for scaling PARA
// double flux_ratio_butanol[3] = { 0.629846,0.556939,0.535483 };

// double flux_ratio_CH2[3] = {0.221792,0.211066,0.225832 };  // 1.1GeV negative flux ratio for scaling PARA
// double flux_ratio_butanol[3] = {0.264991,0.242277,0.209651};

// double flux_ratio_CH2[2] = {1.17317,0.839595};  // 1.3GeV positive flux ratio for scaling PARA
// double flux_ratio_butanol[2] = {1.05876,0.867238 };

// double flux_ratio_CH2[2] = {1.29653,1.06023};  // 1.3GeV negative flux ratio for scaling PARA
// double flux_ratio_butanol[2] = { 1.26422,1.0542};

// double flux_ratio_CH2[3] = { 0.883502,0.948098,0.799394 };  // 1.5GeV positive flux ratio for scaling PARA
// double flux_ratio_butanol[3] = {0.920373,0.836235,0.74416 };

// double flux_ratio_CH2[3] = {1.35844,0.998401,0.807666 };  // 1.5GeV negative flux ratio for scaling PARA
// double flux_ratio_butanol[3] = {1.3026,0.962663,0.885399 };

// double flux_ratio_CH2[2] = {1.93535,1.53025 };  // 1.7GeV positive flux ratio for scaling PARA
// double flux_ratio_butanol[2] = { 1.87958,1.59807  };

// double flux_ratio_CH2[2] = { };  // 1.7GeV negative flux ratio for scaling PARA
// double flux_ratio_butanol[2] = { };

// double flux_ratio_CH2[5] = {1.20124,1.28662,1.18292,1.20703,0.839317 };  // 1.9GeV positive flux ratio for scaling PARA
// double flux_ratio_butanol[5] = {1.31659,1.16608,1.04352,0.928731,0.815808  };

// double flux_ratio_CH2[5] = {1.40767,1.37744,1.2224,0.98128,0.773147 };  // 1.9GeV negative flux ratio for scaling PARA
// double flux_ratio_butanol[5] = {1.45603,1.42382,1.23249,1.06424,0.917127 };

// double flux_ratio_CH2[5] = {1.44532,1.34613,1.26976,1.13251,0.776628 };  // 2.1GeV positive flux ratio for scaling PARA
// double flux_ratio_butanol[5] = {1.28098,1.38163,1.21543,1.06575,0.871218};

// double flux_ratio_CH2[5] = {1.2646,1.04638,0.929401,0.466963,0.925401 };  // 2.1GeV negative flux ratio for scaling PARA
// double flux_ratio_butanol[5] = { 1.18326,1.19632,1.02999,0.852153,0.852037};

// double flux_ratio_CH2[6] = {1.17625,0.91923,0.750992,0.900728,0.934838,0.719083 };  // 2.3GeV positive flux ratio for scaling PARA
// double flux_ratio_butanol[6] = {0.952848,1.01443,1.06914,0.965474,0.916638,0.740935};

// double flux_ratio_CH2[6] = { 0.956999,0.849629,0.871006,1.24237,0.700868,0.705015};  // 2.3GeV negative flux ratio for scaling PARA
// double flux_ratio_butanol[6] = { 0.892061,0.916926,0.912408,0.823409,0.807554,0.702889};

// ------------------------------------------------------------------------------------------------------------------------------

// Target polarisation-------------------------------------------------------------------------------------------

double TARG_pol = 0.0;            // target data
double TARG_pol_stat_error = 0.0; 

// double TARG_pol = 0.90928;            // 0.730GeV, positive target data
// double TARG_pol_stat_error = 0.00123; 

//  double TARG_pol = -0.793584;            // 0.730GeV, negative target data
// double TARG_pol_stat_error = 0.0010868; 

 // double TARG_pol = 0.92645;            // 0.930GeV, positive target data
 // double TARG_pol_stat_error = 0.0057478;

// double TARG_pol = -0.77699;            // 0.930GeV, negative target data
// double TARG_pol_stat_error = 0.0113;

// double TARG_pol = 0.89429;          //  1.1GeV, positive target data
// double TARG_pol_stat_error = 0.004861;

// double TARG_pol = -0.77484;            // 1.1GeV, negative target data
// double TARG_pol_stat_error = 0.008475; 

// double TARG_pol = 0.8672;            // 1.3GeV, positive target data
// double TARG_pol_stat_error = 0.0128; 

// double TARG_pol = -0.7788;            // 1.3GeV, negative target data
// double TARG_pol_stat_error = 0.004295; 

// double TARG_pol = 0.86202;            // 1.5GeV, positive target data
// double TARG_pol_stat_error = 0.0047519; 

// double TARG_pol = -0.784426;            // 1.5GeV, negative target data
// double TARG_pol_stat_error = 0.0030286;

// double TARG_pol = 0.84425;            // 1.7GeV, positive target data
// double TARG_pol_stat_error = 0.002324; 

// double TARG_pol = -0.77212;            // 1.7GeV, negative target data
// double TARG_pol_stat_error = 0.02078; 

// double TARG_pol = 0.84553;            // 1.9GeV, positive target data
// double TARG_pol_stat_error = 0.03837; 

// double TARG_pol = -0.8185;            // 1.9GeV, negative target data
// double TARG_pol_stat_error = 0.05658; 

// double TARG_pol = 0.84055;            // 2.1GeV, positive target data
// double TARG_pol_stat_error = 0.01768; 

// double TARG_pol = -0.76008;            // 2.1GeV, negative target data
// double TARG_pol_stat_error = 0.007724; 

// double TARG_pol = 0.88985;            // 2.3GeV, positive target data
// double TARG_pol_stat_error = 0.01258; 
			     
// double TARG_pol = -0.761904;            // 2.3GeV, negative target data
// double TARG_pol_stat_error = 0.251298; 
			     
// --------------------------------------------------------------------------------------
			     
// ---------------------- Dilution Factor ---------------------------------------------
			     
double dilution_factor_butanol = 0.5315;
double dilution_factor_stat_error_butanol = 0.009;
			     
double dilution_factor_array[10][10];
double dilution_factor_stat_error_array[10][10];

// -------------------------------------------------------------------------------------
			     
// ----------------------- Constants --------------------------------------

double BEAM_pol[5];  //  Average beam polarisation
double BEAM_pol_sys_error[5];  // Systematic error in average beam polarisations

// For the CH2 (unpolarised target)

double PARA_events_CH2=0;  // number of events in each histogram binned in theta and energy
double PERP_events_CH2=0;
double AMO_events_CH2=0;         // number of events in the AMO histogram
double scale_PARA_events_CH2_asym_flux=0; // ratio to scale PARA and PERP histograms
double scale_PARA_events_CH2_asym_flux_error=0;

double PARA_offset_CH2=0;  // The phi offset extracted from the fit of the normalised PARA histogram
double PERP_offset_CH2=0;  // The phi offset extracted from the fit of the normalised PERP histogram
double PARA_flux_ratio_CH2=0;  // F_PARA/F_AMO from the fit of the normalised PARA histogram
double PERP_flux_ratio_CH2=0;  // F_PERP/F_AMO extracted from the fit of the normalised PERP histogram
double PARA_flux_ratio_CH2_error=0;
double PERP_flux_ratio_CH2_error=0;
double offset_CH2=0;       // The average of these two phi offsets
double PARA_phi_offset_CH2[100][100]; // Array to store the offsets for each energy and theta bin
double PARA_phi_offset_CH2_error[100][100];
double PERP_phi_offset_CH2[100][100]; // Array to store the offsets for each energy and theta bin
double PERP_phi_offset_CH2_error[100][100];

int events_CH2[100][100][5];  // total number of events from CH2 histograms

// For the butanol (polarised) target

double PARA_events_butanol=0;  // number of events in each histogram binned in theta and energy
double PERP_events_butanol=0;
double AMO_events_butanol=0;         // number of events in the AMO histogram
double PARA_PERP_ratio_butanol=0;        // ratio of events in PARA and PERP hisotgrams

double scale_PARA_events_butanol_asym_flux=0; // ratio to scale PARA and PERP histograms
double scale_PARA_events_butanol_asym_flux_error=0;

double constrain_sigma[100][100];      // The value of sigma from the CH2 data fit multiplied by average polarisation

int events_butanol[100][100][5];  // total number of events from butanol histograms

float pi = TMath::Pi();

// Arrays to store fit parameters for CH2 (unpolarised target) data:

double Sigma[100][100];
double Sigma_error_stat[100][100];
double Sigma_butanol[100][100];
double Sigma_butanol_error_stat[100][100];
double Sigma_pol[2][100][100];
double Sigma_pol_error_stat[2][100][100];

double CH2_flux_ratio[100][100];
double CH2_flux_ratio_stat_error[100][100];
double CH2_pol_ratio[100][100];
double CH2_pol_ratio_stat_error[100][100];
double CH2_phi_offset[100][100];
double CH2_phi_offset_stat_error[100][100];

double x_angle_CH2[100][100];
double xe_angle_CH2[100][100];
double x_energy_CH2[100][100];
double xe_energy_CH2[100][100];

// Arrays to store fit parameters for butanol (polarised target) data:

double G[100][100];
double G_error_stat[100][100];
double G_pol[2][100][100];
double G_pol_error_stat[2][100][100];
double G_polfull[2][100][100];
double G_polfull_error_stat[2][100][100];
double butanol_flux_ratio[100][100];
double butanol_flux_ratio_stat_error[100][100];
double butanol_pol_ratio[100][100];
double butanol_pol_ratio_stat_error[100][100];
double Sigma_butanol_pol[2][100][100];
double Sigma_butanol_pol_stat_error[2][100][100];
double butanol_phi_offset[100][100];
double butanol_phi_offset_stat_error[100][100];

double x_angle_butanol[100][100];
double xe_angle_butanol[100][100];
double x_energy_butanol[100][100];
double xe_energy_butanol[100][100];

double PARA_flux_ratio_butanol=0;  // F_PARA/F_AMO from the fit of the normalised PARA histogram
double PERP_flux_ratio_butanol=0;  // F_PERP/F_AMO extracted from the fit of the normalised PERP histogram
double PARA_flux_ratio_butanol_error=0;
double PERP_flux_ratio_butanol_error=0;

// ------------------------------------------------------------------------

// ----------------------- Histograms ---------------------------------------

TH1 *PARA_phi_CH2[100][100];
TH1 *PERP_phi_CH2[100][100];
TH1 *PERP_phi_CH2_scaled_flux[100][100];
TH1 *AMO_phi_CH2[100][100];
TH1 *asym_CH2_h[100][100];
TH1 *PARA_asym_CH2[100][100];
TH1 *PERP_asym_CH2[100][100];
TH1 *AMO_asym2_CH2[100][100];

TH1 *PARA_phi_butanol[100][100];
TH1 *PERP_phi_butanol[100][100];
TH1 *AMO_phi_butanol[100][100];
TH1 *PERP_phi_butanol_scaled_flux[100][100];
TH1 *asym_butanol_h[100][100];
TH1 *PARA_asym_butanol[100][100];
TH1 *PERP_asym_butanol[100][100];

// ------------------------------------------------------------------------

// ------------------- Naming histograms ----------------------------------

char *PARA_CH2_name= new char[100];
char *PERP_CH2_name = new char[100];
char *AMO_CH2_name = new char[100];
char *PARA_asym_CH2_name = new char[100];
char *PERP_asym_CH2_name = new char[100];
char *asym_CH2_name = new char[100];
char *asym_CH2_title = new char[100];

char *PARA_butanol_name= new char[100];
char *PERP_butanol_name = new char[100];
char *AMO_butanol_name = new char[100];
char *asym_butanol_name = new char[100];
char *asym_butanol_title = new char[100];
char *PARA_asym_butanol_name = new char[100];
char *PERP_asym_butanol_name = new char[100];

// -----------------------------------------------------------------------

// ------------------------- Fitting functions for CH2 -------------------

// Set up cos2phi function to fit to the asymmetry:

//Old function if PARA and PERP beam polarisations are the same:

TF1 *cos2phi = new TF1("cos2phi","([0]+([1]*cos((2*(x-[2]))*TMath::DegToRad())))",-180,180);

//New function taking into account that PARA and PERP beam polarisations differ:

//TF1 *cos2phi = new TF1("cos2phi","(([0]*[1]*(cos((2*(x-[2]))*TMath::DegToRad())))/(2+([3]*[1]*(cos((2*(x-[2]))*TMath::DegToRad())))))",-180,180);

// Set up functions to fit to normalised phi distributions:

TF1 *PARA_fit_function_CH2 = new TF1("PARA_fit_function_CH2","([0]*(1+([1]*cos((2*(x-[2]))*TMath::DegToRad()))))");

TF1 *PERP_fit_function_CH2 = new TF1("PERP_fit_function_CH2","([0]*(1-([1]*cos((2*(x-[2]))*TMath::DegToRad()))))");

// -------------------------------------------------------------------------

// ------------------------- Fitting functions for butanol --------------
// Set up cos2phi+sin2phi  function to fit to the butanol asymmetry:

//Old function if PARA and PERP beam polarisations are the same:

TF1 *cos2phi_sin2phi = new TF1("cos2phi_sin2phi","([0]+([1]*cos((2*(x-[2]))*TMath::DegToRad()))+([3]*sin((2*(x-[2]))*TMath::DegToRad())))",-180,180);

//New function taking into account that PARA and PERP beam polarisations differ:

//TF1 *cos2phi_sin2phi = new TF1("cos2phi_sin2phi","((([0]*[1]*(cos((2*(x-[2]))*TMath::DegToRad())))+([0]*[3]*(sin((2*(x-[2]))*TMath::DegToRad()))))/(2+([4]*[1]*(cos((2*(x-[2]))*TMath::DegToRad())))+([4]*[3]*(sin((2*(x-[2]))*TMath::DegToRad())))))",-180,180);

// Set up functions to fit to normalised phi distributions:

TF1 *PARA_fit_function_butanol = new TF1("PARA_fit_function_butanol","([0]*(1+([1]*cos((2*(x-[2]))*TMath::DegToRad()))+([3]*sin((2*(x-[2]))*TMath::DegToRad()))))");

TF1 *PERP_fit_function_butanol = new TF1("PERP_fit_function_butanol","([0]*(1-([1]*cos((2*(x-[2]))*TMath::DegToRad()))-([3]*sin((2*(x-[2]))*TMath::DegToRad()))))");

// -------------------------------------------------------------------------

// -------------------- Asymmetry graphs ---------------------------------

TGraphErrors *asym_graph_CH2[100]; // Plot the asymmetry
TGraphErrors *asym_graph_CH2_single[20][100]; // Plot the asymmetry
char *graphname_CH2 = new char[100]; // name of the asymmetry graph
char *graphtitle_CH2 = new char[100]; // title of the asymmetry graph
double X_CH2[100]; // arrays to store points for the graph
double Y_CH2[100];
double Ye_CH2[100];
double Xe_CH2[100];
double X_CH2_single[20][100]; // arrays to store points for the graph
double Y_CH2_single[20][100];
double Ye_CH2_single[20][100];
double Xe_CH2_single[20][100];

TGraphErrors *asym_graph_butanol[100]; // Plot the asymmetry
TGraphErrors *asym_graph_butanol_pol[2][100]; // Plot the asymmetry
TGraphErrors *asym_graph_butanol_single[20][100]; // Plot the asymmetry
char *graphname_butanol = new char[100]; // name of the asymmetry graph
char *graphtitle_butanol = new char[100]; // title of the asymmetry graph
char *graphname_butanol_pol = new char[100]; // name of the asymmetry graph
char *graphtitle_butanol_pol = new char[100]; // title of the asymmetry graph

double X_butanol[100]; // arrays to store points for the graph
double Y_butanol[100];
double Ye_butanol[100];
double Xe_butanol[100];
double X_butanol_single[20][100]; // arrays to store points for the graph
double Y_butanol_single[20][100];
double Ye_butanol_single[20][100];
double Xe_butanol_single[20][100];

TGraphErrors *asym_graph_sigma_butanol[100]; // Plot the asymmetry
char *graphname_sigma_butanol = new char[100]; // name of the asymmetry graph
char *graphtitle_sigma_butanol = new char[100]; // title of the asymmetry graph
double X_sigma_butanol[100]; // arrays to store points for the graph
double Y_sigma_butanol[100];
double Ye_sigma_butanol[100];
double Xe_sigma_butanol[100];

float energy_label[100];  // array to srore energy values for naming asymmetry graph
float theta_label[100];   // array to store theta values for naming asymmetry graph
float theta_ax[100];   // array to store theta values for naming asymmetry graph
float etheta_ax[100];   // array to store theta values for naming asymmetry graph
int target_option;       // Is the target polarised?
int plotting_option;  // Do we want the asymmetry graph to have energy or theta on x axis?

TCanvas *c_CH2[100][100];    // Canvas on which to plot the asymmetries
TCanvas *c_butanol[100][100];    // Canvas on which to plot the asymmetries
TCanvas *c_G[100];    // Canvas on which to plot the asymmetries
TCanvas *c_Sigma[100];    // Canvas on which to plot the asymmetries
TCanvas *c_Sigma_butanol[100];    // Canvas on which to plot the asymmetries
      
// ------------------------------------------------------------------------

// ------------------------- Functions used in code ------------------------

void Asymmetry_plotting_theta(char const *outfile); // Plots the asymmetries as graphs
void Asymmetry_plotting_energy(char const *outfile); // Plots the asymmetries as graphs
void Labels(); // Creates labels for the asymmetries
void ReadParameters(int peak, int polar);
// ------------------------- Files to save flux ratios -----------------------



void ReadParameters(int peak, int polar){
  TChain input_chain("FROST_asym");
  input_chain.Add("/phys/linux/lzana/rootbeer2.2/Analysis/Data_structure/frost_asym.root");

  int n_Ebin,polarization;
  float Energy_peak, E_cms_min, E_cms_max, w_Ebin, width_Ebin;
  float v_para_pol[6] ; 
  float v_para_pol_sys[6] ; 
  float v_perp_pol[6] ; 
  float v_perp_pol_sys[6] ; 
  float v_flux_ratio_CH2[6] ; 
  float v_flux_ratio_butanol[6] ; 
  float v_TARG_pol ; 
  float v_TARG_pol_stat_error ;
  
  input_chain.SetBranchAddress("Energy_peak",&Energy_peak);
  input_chain.SetBranchAddress("Energy_cms_min",&E_cms_min);
  input_chain.SetBranchAddress("Energy_cms_max",&E_cms_max);
  input_chain.SetBranchAddress("no_Ebin",&n_Ebin);
  input_chain.SetBranchAddress("width_Ebin",&w_Ebin);
  input_chain.SetBranchAddress("polariz",&polarization);
  input_chain.SetBranchAddress("para_pol",v_para_pol);
  input_chain.SetBranchAddress("para_pol_sys",v_para_pol_sys);
  input_chain.SetBranchAddress("perp_pol",v_perp_pol);
  input_chain.SetBranchAddress("perp_pol_sys",v_perp_pol_sys);
  input_chain.SetBranchAddress("flux_ratio_CH2",v_flux_ratio_CH2);
  input_chain.SetBranchAddress("flux_ratio_butanol",v_flux_ratio_butanol);
  input_chain.SetBranchAddress("TARG_pol",&v_TARG_pol);
  input_chain.SetBranchAddress("TARG_pol_stat_error",&v_TARG_pol_stat_error);


  // Energy_min = (peak - 200.)/1000; // Energy in GeV
  // Energy_max = (peak + 100.)/1000; // Energy in GeV

  Int_t nentries = (Int_t)input_chain.GetEntries();
  for (int i=0; i<nentries ; i++) {
    input_chain.GetEntry(i);
    printf("Entry %i Energy_peak=%f (requested=%f ) polariz=%i (requested=%i )\n",i,Energy_peak,peak,polarization,polar);
    if (Energy_peak == peak && polarization == polar) {
      Energy_at_peak = Energy_peak;
      Energy_cms_min = E_cms_min /1000;
      Energy_cms_max = E_cms_max /1000; 
      no_Ebin = n_Ebin;
      Ebin_width = w_Ebin /1000;
      for (int j=0 ; j<6; j++) {
	PARA_pol[j] = v_para_pol[j];
	PARA_pol_sys_error[j] = v_para_pol_sys[j];
	PERP_pol[j] = v_perp_pol[j];
	PERP_pol_sys_error[j] = v_perp_pol_sys[j];
	flux_ratio_CH2[j] = v_flux_ratio_CH2[j];
	flux_ratio_butanol[j] = v_flux_ratio_butanol[j];
      }
      TARG_pol = v_TARG_pol;
      TARG_pol_stat_error = v_TARG_pol_stat_error;

      printf("Reading parameters from FROST data structure file /home/zana/rootbeer2.1/Analysis/Data_structure/frost_structure.root \n for Energy peak of %f and polarization of %d . Energy_cms_min = %f . width_Ebin =%f  \n ",peak,polar,Energy_cms_min,w_Ebin); 

      
    }

  }


}

// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
// Issues with polarization configuration:
// output_730_neg (perp)
// output_730_pos (perp)
void get_flux(char const *infile);

void asymmetry(char const *infile,char const *outfile,int peak, int polar){
  //void asymmetry(char *infile,TFile *dfile,char *outfile){
  
  TFile *input = new TFile(infile);
  char pol_file_730name[128];
  if (peak == 730 && polar == -1) sprintf(pol_file_730name,"output_pol05_730_neg.root") ;
  if (peak == 730 && polar == 1) sprintf(pol_file_730name,"output_pol05_730_pos.root") ;

  TFile *file730 = new TFile(pol_file_730name);

  TFile *output = new TFile(outfile,"RECREATE");

  // get parameters for this peak and polarization
  ReadParameters(peak,polar);

  // set the flux for each angle and energy
  get_flux(infile);

  int at_pol = int( polar + 1 ) /2;
  printf(" Polarization array = %d \n",at_pol);
  
  min_energy = Energy_cms_min;

  plotting_option =1;  // Plot theta on the x axis. 
                       // 1 = plotting energy on the x axis

  // Loop through all the PARA, PERP and AMO histograms produced by the analysis code

  // Open the file where the dilution factor is kept and put values into arrays
  double Energy_at_peak = double(peak)/1000.;
  char edgeSettingChar[128];
  char polarSettingChar[128];

  if (polar == 1) sprintf(polarSettingChar,"positive");
  if (polar == -1) sprintf(polarSettingChar,"negative");
  
  
  if (Energy_at_peak == 0.73) {
    sprintf(edgeSettingChar,"0_730GeV");
  }
  else if (Energy_at_peak == 0.93) {
    sprintf(edgeSettingChar,"0_930GeV");
  }
  else if (Energy_at_peak == 1.1) {
    sprintf(edgeSettingChar,"1_1GeV");
  }
  else if (Energy_at_peak == 1.3) {
    sprintf(edgeSettingChar,"1_3GeV");
  }
  else if (Energy_at_peak == 1.5) {
    sprintf(edgeSettingChar,"1_5GeV");
  }
  else if (Energy_at_peak == 1.7) {
    sprintf(edgeSettingChar,"1_7GeV");
  }
  else if (Energy_at_peak == 1.9) {
    sprintf(edgeSettingChar,"1_9GeV");
  }
  else if (Energy_at_peak == 2.1) {
    sprintf(edgeSettingChar,"2_1GeV");
  }
  else if (Energy_at_peak == 2.3) {
    sprintf(edgeSettingChar,"2_3GeV");
  }

  char name_dilutionfile[128];

  //  sprintf(name_dilutionfile,"v22_dilution_factor/%s/%s_%s_dilution_factor.txt",edgeSettingChar,edgeSettingChar,polarSettingChar);

  sprintf(name_dilutionfile,"DilutionFactor/dilution_factor_%d_%d.txt",peak,at_pol);
  ifstream dilutionfile;

  dilutionfile.open(name_dilutionfile);

  printf("Opening dilution file %s \n",name_dilutionfile);
  int j,k;
  double l,m;     // j = ebin, k = thbin, l=f, m=error

  while(!dilutionfile.eof()){

    dilutionfile >> j >> k >> l >> m;
    
    printf("energy=%d  theta=%d  l=%f  m=%f \n",j,k,l,m);
    dilution_factor_array[k][j]=l;
    dilution_factor_stat_error_array[k][j]=m;

  }
  TH1F * PARA_pol_histo[no_Ebin];
  TH1F * PERP_pol_histo[no_Ebin];

				  
  char histo_pol_beam[128];

  for(int e=0; e<no_Ebin; e++){
    
    for(int t=0; t<no_Thbin; t++){

      // To extract the asymmetry we need the following combinations of beam polarisation: for 730 GeV many failed with wrong fit at the time of data taking (for the PERP configuration). So I just get the values of the polarization where the fit was fine, and use it for the data where the beam polarization was not used as weight for the event. For the other values, the beam polarization has been included directly in the yields plots on an event by event base.
      
      if (peak == 730) {
	sprintf(histo_pol_beam,"para_pol_h%d",e+1);
	PARA_pol_histo[e] = (TH1F*)file730->Get(histo_pol_beam);
	sprintf(histo_pol_beam,"perp_pol_h%d",e+1);
	PERP_pol_histo[e] = (TH1F*)file730->Get(histo_pol_beam);
	PARA_pol[e] = PARA_pol_histo[e]->GetMean();
	PARA_pol_sys_error[e] = PARA_pol_histo[e]->GetRMS();
	PERP_pol[e] = PERP_pol_histo[e]->GetMean();
	PERP_pol_sys_error[e] = PERP_pol_histo[e]->GetRMS();
      }

      BEAM_pol[e] = ((PARA_pol[e]+PERP_pol[e])/2);
      BEAM_pol_sys_error[e] = TMath::Sqrt((PARA_pol_sys_error[e]*PARA_pol_sys_error[e])+(PERP_pol_sys_error[e]*PERP_pol_sys_error[e]));
   
      // Get name of histogram and save to PARA_name:
      
      sprintf(PARA_CH2_name,"phi_pip_PARA_Target3_E%d_Th%d",e+1,t+1);
      sprintf(PARA_butanol_name,"phi_pip_PARA_Target1_E%d_Th%d",e+1,t+1);

      // Find histogram:
      
      PARA_phi_CH2[t][e]=(TH1*)input->Get(PARA_CH2_name);
      PARA_phi_butanol[t][e]=(TH1*)input->Get(PARA_butanol_name);
   
      // Rebin for easy asymmetry viewing:
      
      PARA_phi_CH2[t][e]->Rebin(5);
      PARA_phi_butanol[t][e]->Rebin(5);

      // For error calculation later:
      
      PARA_phi_CH2[t][e]->Sumw2();
      PARA_phi_butanol[t][e]->Sumw2();
      
      // Find the number of events in the histogram:
      
      PARA_events_CH2 = PARA_phi_CH2[t][e]->Integral();
      PARA_events_butanol = PARA_phi_butanol[t][e]->Integral();
      
      cout << "Energy: " << e << ", Theta: " << t << ", Number of PARA events for CH2: " << PARA_events_CH2 << endl;
      cout << "Energy: " << e << ", Theta: " << t << ", Number of PARA events for butanol: " << PARA_events_butanol << endl;

      // Repeat the above for the PERP histogram:
      
      sprintf(PERP_CH2_name,"phi_pip_PERP_Target3_E%d_Th%d",e+1,t+1);
      sprintf(PERP_butanol_name,"phi_pip_PERP_Target1_E%d_Th%d",e+1,t+1);
      
      PERP_phi_CH2[t][e]=(TH1*)input->Get(PERP_CH2_name);
      PERP_phi_butanol[t][e]=(TH1*)input->Get(PERP_butanol_name);
      
      PERP_phi_CH2[t][e]->Rebin(5);
      PERP_phi_butanol[t][e]->Rebin(5);

      PERP_phi_CH2[t][e]->Sumw2();
      PERP_phi_butanol[t][e]->Sumw2();
      
      PERP_events_CH2 = PERP_phi_CH2[t][e]->Integral();
      PERP_events_butanol = PERP_phi_butanol[t][e]->Integral();
      
      cout << "Energy: " << e << ", Theta: " << t << ", Number of PERP events for CH2: " << PERP_events_CH2 << endl;
      cout << "Energy: " << e << ", Theta: " << t << ", Number of PERP events for butanol: " << PERP_events_butanol << endl;

      // Repeat the above with the AMO histogram:

      sprintf(AMO_CH2_name,"phi_pip_AMO_Target3_E%d_Th%d",e+1,t+1);
      sprintf(AMO_butanol_name,"phi_pip_AMO_Target1_E%d_Th%d",e+1,t+1);
      
      AMO_phi_CH2[t][e]=(TH1*)input->Get(AMO_CH2_name);
      AMO_phi_butanol[t][e]=(TH1*)input->Get(AMO_butanol_name);
      
      AMO_phi_CH2[t][e]->Rebin(5);
      AMO_phi_butanol[t][e]->Rebin(5);

      AMO_phi_CH2[t][e]->Sumw2();
      AMO_phi_butanol[t][e]->Sumw2();
      
      AMO_events_CH2 = AMO_phi_CH2[t][e]->Integral();
      AMO_events_butanol = AMO_phi_butanol[t][e]->Integral();
      
      cout << "Energy: " << e << ", Theta: " << t << ", Number of AMO events for CH2: " << AMO_events_CH2 << endl;
      cout << "Energy: " << e << ", Theta: " << t << ", Number of AMO events for butanol: " << AMO_events_butanol << endl;
      
      // Put the number of events into arrays
      
      events_CH2[t][e][0] = PARA_events_CH2;
      events_CH2[t][e][1] = PERP_events_CH2;
      events_CH2[t][e][2] = PARA_events_CH2 + PERP_events_CH2;
      events_CH2[t][e][3] = AMO_events_CH2;

      events_butanol[t][e][0] = PARA_events_butanol;
      events_butanol[t][e][1] = PERP_events_butanol;
      events_butanol[t][e][2] = PARA_events_butanol + PERP_events_butanol;
      events_butanol[t][e][3] = AMO_events_butanol;

      //---------------------------------------------------- Obtaining Sigma Asymmetry from CH2 data ------------------------------------------------------------

      // Obtain the asymmetry from the histograms if they contain points:
      
      if ((PARA_events_CH2 ==0) || (PERP_events_CH2 ==0)){
	tot_bin_energy_ch2[t][e] = 0;
	
	cout << "Nothing in either PARA or PERP or AMO for CH2 energy bin " << e << " theta bin " << t <<  " so no asymmetry will be created" << endl;
	
      }

      else {
	tot_bin_energy_ch2[t][e] = 1;
	cout << "energy bin: " << e << " theta bin: " << t << endl;

	// Now need to scale PARA and PERP histograms to eachother before extracting the asymmetry

	// Scale using the flux ratios:

	PERP_phi_CH2_scaled_flux[t][e] = (TH1*)PERP_phi_CH2[t][e]->Clone();
	PERP_phi_CH2_scaled_flux[t][e]->Sumw2();

	PERP_phi_CH2_scaled_flux[t][e]->Scale(val_flux_ratio_CH2[t][e]);
	
	// Check that the scaling worked:
	
	cout << "Number of PERP events in CH2 histogram after scaling with PARA (flux method): " << PERP_phi_CH2_scaled_flux[t][e]->Integral() << endl;

	// Now extract the asymmetry:
       
	// asym_CH2_h[t][e] = PERP_phi_CH2_scaled_flux[t][e]->GetAsymmetry(PARA_phi_CH2[t][e]);
	asym_CH2_h[t][e] = PARA_phi_CH2[t][e]->GetAsymmetry(PERP_phi_CH2_scaled_flux[t][e]);

	// Name the asymmetry

	sprintf(asym_CH2_name,"Asymmetry_CH2_Target_Ebin%d_Thbin%d_pol%d",e+1,t+1,at_pol);
	asym_CH2_h[t][e]->SetName(asym_CH2_name);

	// Create labels for energies and theta values

	Labels();

	// Create title for asymmetry

	sprintf(asym_CH2_title,"Asymmetry_%fGeV_to_%fGeV_%ddegrees_to_%ddegrees_for_CH2_target",energy_label[e],energy_label[e+1],theta_label[t+1],theta_label[t]);  

	asym_CH2_h[t][e]->SetTitle(asym_CH2_title);

	// Set options for statistics

	gStyle->SetOptStat("e");
	gStyle->SetOptFit(1111);

	// Set the phi-offset to zero

	cos2phi->FixParameter(2,0);

	// Now fit the cos2phi function to the asymmetry:
	
	asym_CH2_h[t][e]->Fit("cos2phi","QB");

	// Plot the asymmetry

	c_CH2[t][e]= new TCanvas();
	c_CH2[t][e]->cd();
	asym_CH2_h[t][e]->Draw();
	c_CH2[t][e]->Update();

	// Extract sigma from the fit parameters to plot later:
	if (peak == 730 && BEAM_pol[e]!=0.0 )  Sigma_pol[at_pol][t][e]=-((cos2phi->GetParameter(1))/BEAM_pol[e]);
	else Sigma_pol[at_pol][t][e]=-(cos2phi->GetParameter(1));
	// Calculate the statistical error in Sigma.

	if (peak == 730 && BEAM_pol[e]!=0.0)  Sigma_pol_error_stat[at_pol][t][e]=((cos2phi->GetParError(1))/BEAM_pol[e]);
	else Sigma_pol_error_stat[at_pol][t][e]=(cos2phi->GetParError(1)) ;
	// Get theta angle in degrees from bin number and the associated error
 
	x_angle_CH2[t][e]=theta_ax[t];
	xe_angle_CH2[t][e]=etheta_ax[t];              // Apparently the error in the xaxis is the bin width

	// Get energy in GeV from bin number and the associated error

	x_energy_CH2[t][e]= (e*Ebin_width)+(min_energy);

	xe_energy_CH2[t][e] = Ebin_width/2;                     // Apparently the error in the x axis is the bin width

	// Write to file
	output->cd();
	asym_CH2_h[t][e]->Write();

      } // End of loop selecting only histograms containing events

      //---------------------------------------------------------------------------------------------------------------------------------------------------------

      //---------------------------------------------------- Obtaining G Asymmetry from Butanol data ------------------------------------------------------------
      
      // Obtain the asymmetry from the histograms if they contain points:
      
      if ((PARA_events_butanol ==0) || (PERP_events_butanol ==0)){
	tot_bin_energy_butanol[t][e] = 0.0;
	
       	cout << "Nothing in either PARA or PERP for butanol energy bin " << e << " theta bin " << t <<  " so no asymmetry will be created" << endl;
	
      }
      
      else {
	tot_bin_energy_butanol[t][e] = 1.0;
	// Scale the PARA and PERP butanol histograms using the scale factors from the butanol flux ratios

	PERP_phi_butanol_scaled_flux[t][e] = (TH1*)PERP_phi_butanol[t][e]->Clone();
	PERP_phi_butanol_scaled_flux[t][e]->Sumw2();

	PERP_phi_butanol_scaled_flux[t][e]->Scale(val_flux_ratio_butanol[t][e]);
	
	// Check that the scaling worked:
	
	cout << "Number of PERP events in butanol histogram after scaling with PARA (flux method): " << PERP_phi_butanol_scaled_flux[t][e]->Integral() << endl;

	// Extract the asymmetry:
    
        // asym_butanol_h[t][e] = PERP_phi_butanol_scaled_flux[t][e]->GetAsymmetry(PARA_phi_butanol[t][e]);
        asym_butanol_h[t][e] = PARA_phi_butanol[t][e]->GetAsymmetry(PERP_phi_butanol_scaled_flux[t][e]);
	
	// Name the asymmetry

	sprintf(asym_butanol_name,"Asymmetry_Butanol_Target_Ebin%d_Thbin%d_pol%d",e+1,t+1,at_pol);
	asym_butanol_h[t][e]->SetName(asym_butanol_name);

	// Create labels for energies and theta values

	Labels();

	// Create title for asymmetry

	sprintf(asym_butanol_title,"Asymmetry_%fGeV_to_%fGeV_%ddegrees_to_%ddegrees_for_butanol_target",energy_label[e],energy_label[e+1],theta_label[t+1],theta_label[t]);  

	asym_butanol_h[t][e]->SetTitle(asym_butanol_title);

	// Set statistics and fitting options

	gStyle->SetOptStat("e");
	gStyle->SetOptFit(1111);

	// Set the phi-offset to zero

	cos2phi_sin2phi->FixParameter(2,0);

	// Now fit the cos2phi-sin2phi function to the asymmetry:
	
	asym_butanol_h[t][e] ->Fit("cos2phi_sin2phi","QB");

	// Plot the asymmetry

	c_butanol[t][e]= new TCanvas();
	c_butanol[t][e]->cd();
	asym_butanol_h[t][e]->Draw();
	c_butanol[t][e]->Update();

	// Extract G from the fit parameters to plot later:
	
	if (peak==730 && BEAM_pol[e]!=0.0) G_pol[at_pol][t][e]=((cos2phi_sin2phi->GetParameter(3))/(dilution_factor_array[t][e]*TARG_pol*BEAM_pol[e]));
	else G_pol[at_pol][t][e]= ((cos2phi_sin2phi->GetParameter(3))/(dilution_factor_array[t][e]*TARG_pol)) ;
	// Calculate the statistical error in G.
	// Don't forget need error for dilution factor!

	G_pol_error_stat[at_pol][t][e]= TMath::Abs(G_pol[at_pol][t][e])*(TMath::Sqrt(((cos2phi_sin2phi->GetParError(3)/cos2phi_sin2phi->GetParameter(3))*(cos2phi_sin2phi->GetParError(3)/cos2phi_sin2phi->GetParameter(3)))+((TARG_pol_stat_error/TARG_pol)*(TARG_pol_stat_error/TARG_pol))+((dilution_factor_stat_error_array[t][e]/dilution_factor_array[t][e])*(dilution_factor_stat_error_array[t][e]/dilution_factor_array[t][e]))));
	printf("!!!!!!!!!!!!   G for pol %d theta bin %d energy bin %d = %.3e  +/- %.3e \n",at_pol,t,e,G_pol[at_pol][t][e],G_pol_error_stat[at_pol][t][e]);
	// Extract Sigma from the butanol histogram to compare to the sigma value from CH2

	if (peak == 730 && BEAM_pol[e]!=0.0) {
	  Sigma_butanol_pol[at_pol][t][e]= ((cos2phi_sin2phi->GetParameter(1))/BEAM_pol[e]);
	  Sigma_butanol_pol_stat_error[at_pol][t][e]= ((cos2phi_sin2phi->GetParError(1))/BEAM_pol[e]);
	}
	else {
	  Sigma_butanol_pol[at_pol][t][e]= ((cos2phi_sin2phi->GetParameter(1)));
	  Sigma_butanol_pol_stat_error[at_pol][t][e]= ((cos2phi_sin2phi->GetParError(1)));
	}
	// Get theta angle in degrees from bin number and the associated error
 
	x_angle_butanol[t][e]=theta_ax[t];
	xe_angle_butanol[t][e]=etheta_ax[t];              // Apparently the error in the xaxis is the bin width

	// Get energy in GeV from bin number and the associated error

	x_energy_butanol[t][e+at_Ebin]= (e*Ebin_width)+(min_energy)+Ebin_width/2;
	printf("theta = %d ; energy = %d ; at energy = %d ; W= %.3f \n",t,e,e+at_Ebin,x_energy_butanol[t][e+at_Ebin]);
	xe_energy_butanol[t][e+at_Ebin] = Ebin_width/2;                     // Apparently the error in the x axis is the bin width

	// Write to file
	output->cd();
	asym_butanol_h[t][e]->Write();

      } // End of loop selecting only histograms containing events

      //---------------------------------------------------------------------------------------------------------------------------------------------------------

    } // End of loop over theta

  } // End of loop over energy

  // Now plot the asymmetries!

  // Asymmetry_plotting_theta();

  // Close all open files


  input->Close();
  output->Close();

} // End of asymmetry() function








// --------------------------------------------------------------------------------------



void asymmetry_pi0(char const *infile,char const *outfile,int peak, int polar){
  //void asymmetry(char *infile,TFile *dfile,char *outfile){
  
  TFile *input = new TFile(infile);
  char pol_file_730name[128];
  if (peak == 730 && polar == -1) sprintf(pol_file_730name,"output_pi0_730_neg.root") ;
  if (peak == 730 && polar == 1) sprintf(pol_file_730name,"output_pi0_730_pos.root") ;

  TFile *file730 = new TFile(pol_file_730name);

  TFile *output = new TFile(outfile,"RECREATE");

  // get parameters for this peak and polarization
  ReadParameters(peak,polar);

  // set the flux for each angle and energy
  get_flux(infile);

  int at_pol = int( polar + 1 ) /2;
  printf(" Polarization array = %d \n",at_pol);
  
  min_energy = Energy_cms_min;

  plotting_option =1;  // Plot theta on the x axis. 
                       // 1 = plotting energy on the x axis

  // Loop through all the PARA, PERP and AMO histograms produced by the analysis code

  // Open the file where the dilution factor is kept and put values into arrays
  double Energy_at_peak = double(peak)/1000.;
  char edgeSettingChar[128];
  char polarSettingChar[128];

  if (polar == 1) sprintf(polarSettingChar,"positive");
  if (polar == -1) sprintf(polarSettingChar,"negative");
  
  
  if (Energy_at_peak == 0.73) {
    sprintf(edgeSettingChar,"0_730GeV");
  }
  else if (Energy_at_peak == 0.93) {
    sprintf(edgeSettingChar,"0_930GeV");
  }
  else if (Energy_at_peak == 1.1) {
    sprintf(edgeSettingChar,"1_1GeV");
  }
  else if (Energy_at_peak == 1.3) {
    sprintf(edgeSettingChar,"1_3GeV");
  }
  else if (Energy_at_peak == 1.5) {
    sprintf(edgeSettingChar,"1_5GeV");
  }
  else if (Energy_at_peak == 1.7) {
    sprintf(edgeSettingChar,"1_7GeV");
  }
  else if (Energy_at_peak == 1.9) {
    sprintf(edgeSettingChar,"1_9GeV");
  }
  else if (Energy_at_peak == 2.1) {
    sprintf(edgeSettingChar,"2_1GeV");
  }
  else if (Energy_at_peak == 2.3) {
    sprintf(edgeSettingChar,"2_3GeV");
  }

  char name_dilutionfile[128];

  //  sprintf(name_dilutionfile,"v22_dilution_factor/%s/%s_%s_dilution_factor.txt",edgeSettingChar,edgeSettingChar,polarSettingChar);

  sprintf(name_dilutionfile,"DilutionFactor/dilution_factor_%d_%d.txt",peak,at_pol);
  ifstream dilutionfile;

  dilutionfile.open(name_dilutionfile);

  printf("Opening dilution file %s \n",name_dilutionfile);
  int j,k;
  double l,m;     // j = ebin, k = thbin, l=f, m=error

  while(!dilutionfile.eof()){

    dilutionfile >> j >> k >> l >> m;
    
    printf("energy=%d  theta=%d  l=%f  m=%f \n",j,k,l,m);
    dilution_factor_array[k][j]=l;
    dilution_factor_stat_error_array[k][j]=m;

  }
  TH1F * PARA_pol_histo[no_Ebin];
  TH1F * PERP_pol_histo[no_Ebin];

				  
  char histo_pol_beam[128];

  for(int e=0; e<no_Ebin; e++){
    
    for(int t=0; t<no_Thbin; t++){

      // To extract the asymmetry we need the following combinations of beam polarisation: for 730 GeV many failed with wrong fit at the time of data taking (for the PERP configuration). So I just get the values of the polarization where the fit was fine, and use it for the data where the beam polarization was not used as weight for the event. For the other values, the beam polarization has been included directly in the yields plots on an event by event base.
      
      if (peak == 730) {
	sprintf(histo_pol_beam,"para_pol_h%d",e+1);
	PARA_pol_histo[e] = (TH1F*)file730->Get(histo_pol_beam);
	sprintf(histo_pol_beam,"perp_pol_h%d",e+1);
	PERP_pol_histo[e] = (TH1F*)file730->Get(histo_pol_beam);
	PARA_pol[e] = PARA_pol_histo[e]->GetMean();
	PARA_pol_sys_error[e] = PARA_pol_histo[e]->GetRMS();
	PERP_pol[e] = PERP_pol_histo[e]->GetMean();
	PERP_pol_sys_error[e] = PERP_pol_histo[e]->GetRMS();
      }

      BEAM_pol[e] = ((PARA_pol[e]+PERP_pol[e])/2);
      BEAM_pol_sys_error[e] = TMath::Sqrt((PARA_pol_sys_error[e]*PARA_pol_sys_error[e])+(PERP_pol_sys_error[e]*PERP_pol_sys_error[e]));
   
      // Get name of histogram and save to PARA_name:
      
      sprintf(PARA_CH2_name,"phi_pip_PARA_Target3_E%d_Th%d",e+1,t+1);
      sprintf(PARA_butanol_name,"phi_pip_PARA_Target1_E%d_Th%d",e+1,t+1);

      // Find histogram:
      
      PARA_phi_CH2[t][e]=(TH1*)input->Get(PARA_CH2_name);
      PARA_phi_butanol[t][e]=(TH1*)input->Get(PARA_butanol_name);
   
      // Rebin for easy asymmetry viewing:
      
      PARA_phi_CH2[t][e]->Rebin(5);
      PARA_phi_butanol[t][e]->Rebin(5);

      // For error calculation later:
      
      PARA_phi_CH2[t][e]->Sumw2();
      PARA_phi_butanol[t][e]->Sumw2();
      
      // Find the number of events in the histogram:
      
      PARA_events_CH2 = PARA_phi_CH2[t][e]->Integral();
      PARA_events_butanol = PARA_phi_butanol[t][e]->Integral();
      
      cout << "Energy: " << e << ", Theta: " << t << ", Number of PARA events for CH2: " << PARA_events_CH2 << endl;
      cout << "Energy: " << e << ", Theta: " << t << ", Number of PARA events for butanol: " << PARA_events_butanol << endl;

      // Repeat the above for the PERP histogram:
      
      sprintf(PERP_CH2_name,"phi_pip_PERP_Target3_E%d_Th%d",e+1,t+1);
      sprintf(PERP_butanol_name,"phi_pip_PERP_Target1_E%d_Th%d",e+1,t+1);
      
      PERP_phi_CH2[t][e]=(TH1*)input->Get(PERP_CH2_name);
      PERP_phi_butanol[t][e]=(TH1*)input->Get(PERP_butanol_name);
      
      PERP_phi_CH2[t][e]->Rebin(5);
      PERP_phi_butanol[t][e]->Rebin(5);

      PERP_phi_CH2[t][e]->Sumw2();
      PERP_phi_butanol[t][e]->Sumw2();
      
      PERP_events_CH2 = PERP_phi_CH2[t][e]->Integral();
      PERP_events_butanol = PERP_phi_butanol[t][e]->Integral();
      
      cout << "Energy: " << e << ", Theta: " << t << ", Number of PERP events for CH2: " << PERP_events_CH2 << endl;
      cout << "Energy: " << e << ", Theta: " << t << ", Number of PERP events for butanol: " << PERP_events_butanol << endl;

      // Repeat the above with the AMO histogram:

      sprintf(AMO_CH2_name,"phi_pip_AMO_Target3_E%d_Th%d",e+1,t+1);
      sprintf(AMO_butanol_name,"phi_pip_AMO_Target1_E%d_Th%d",e+1,t+1);
      
      AMO_phi_CH2[t][e]=(TH1*)input->Get(AMO_CH2_name);
      AMO_phi_butanol[t][e]=(TH1*)input->Get(AMO_butanol_name);
      
      AMO_phi_CH2[t][e]->Rebin(5);
      AMO_phi_butanol[t][e]->Rebin(5);

      AMO_phi_CH2[t][e]->Sumw2();
      AMO_phi_butanol[t][e]->Sumw2();
      
      AMO_events_CH2 = AMO_phi_CH2[t][e]->Integral();
      AMO_events_butanol = AMO_phi_butanol[t][e]->Integral();
      
      cout << "Energy: " << e << ", Theta: " << t << ", Number of AMO events for CH2: " << AMO_events_CH2 << endl;
      cout << "Energy: " << e << ", Theta: " << t << ", Number of AMO events for butanol: " << AMO_events_butanol << endl;
      
      // Put the number of events into arrays
      
      events_CH2[t][e][0] = PARA_events_CH2;
      events_CH2[t][e][1] = PERP_events_CH2;
      events_CH2[t][e][2] = PARA_events_CH2 + PERP_events_CH2;
      events_CH2[t][e][3] = AMO_events_CH2;

      events_butanol[t][e][0] = PARA_events_butanol;
      events_butanol[t][e][1] = PERP_events_butanol;
      events_butanol[t][e][2] = PARA_events_butanol + PERP_events_butanol;
      events_butanol[t][e][3] = AMO_events_butanol;

      //---------------------------------------------------- Obtaining Sigma Asymmetry from CH2 data ------------------------------------------------------------

      // Obtain the asymmetry from the histograms if they contain points:
      
      if ((PARA_events_CH2 ==0) || (PERP_events_CH2 ==0)){
	tot_bin_energy_ch2[t][e] = 0;
	
	cout << "Nothing in either PARA or PERP or AMO for CH2 energy bin " << e << " theta bin " << t <<  " so no asymmetry will be created" << endl;
	
      }

      else {
	tot_bin_energy_ch2[t][e] = 1;
	cout << "energy bin: " << e << " theta bin: " << t << endl;

	// Now need to scale PARA and PERP histograms to eachother before extracting the asymmetry

	// Scale using the flux ratios:

	PERP_phi_CH2_scaled_flux[t][e] = (TH1*)PERP_phi_CH2[t][e]->Clone();
	PERP_phi_CH2_scaled_flux[t][e]->Sumw2();

	PERP_phi_CH2_scaled_flux[t][e]->Scale(val_flux_ratio_CH2[t][e]);
	
	// Check that the scaling worked:
	
	cout << "Number of PERP events in CH2 histogram after scaling with PARA (flux method): " << PERP_phi_CH2_scaled_flux[t][e]->Integral() << endl;

	// Now extract the asymmetry:
       
	// asym_CH2_h[t][e] = PERP_phi_CH2_scaled_flux[t][e]->GetAsymmetry(PARA_phi_CH2[t][e]);
	asym_CH2_h[t][e] = PARA_phi_CH2[t][e]->GetAsymmetry(PERP_phi_CH2_scaled_flux[t][e]);

	// Name the asymmetry

	sprintf(asym_CH2_name,"Asymmetry_CH2_Target_Ebin%d_Thbin%d_pol%d",e+1,t+1,at_pol);
	asym_CH2_h[t][e]->SetName(asym_CH2_name);

	// Create labels for energies and theta values

	Labels();

	// Create title for asymmetry

	sprintf(asym_CH2_title,"Asymmetry_%fGeV_to_%fGeV_%ddegrees_to_%ddegrees_for_CH2_target",energy_label[e],energy_label[e+1],theta_label[t+1],theta_label[t]);  

	asym_CH2_h[t][e]->SetTitle(asym_CH2_title);

	// Set options for statistics

	gStyle->SetOptStat("e");
	gStyle->SetOptFit(1111);

	// Set the phi-offset to zero

	cos2phi->FixParameter(2,0);

	// Now fit the cos2phi function to the asymmetry:
	
	asym_CH2_h[t][e]->Fit("cos2phi","QB");

	// Plot the asymmetry

	c_CH2[t][e]= new TCanvas();
	c_CH2[t][e]->cd();
	asym_CH2_h[t][e]->Draw();
	c_CH2[t][e]->Update();

	// Extract sigma from the fit parameters to plot later:
	if (peak == 730 && BEAM_pol[e]!=0.0 )  Sigma_pol[at_pol][t][e]=-((cos2phi->GetParameter(1))/BEAM_pol[e]);
	else Sigma_pol[at_pol][t][e]=-(cos2phi->GetParameter(1));
	// Calculate the statistical error in Sigma.

	if (peak == 730 && BEAM_pol[e]!=0.0)  Sigma_pol_error_stat[at_pol][t][e]=((cos2phi->GetParError(1))/BEAM_pol[e]);
	else Sigma_pol_error_stat[at_pol][t][e]=(cos2phi->GetParError(1)) ;
	// Get theta angle in degrees from bin number and the associated error
 
	x_angle_CH2[t][e]=theta_ax[t];
	xe_angle_CH2[t][e]=etheta_ax[t];              // Apparently the error in the xaxis is the bin width

	// Get energy in GeV from bin number and the associated error

	x_energy_CH2[t][e]= (e*Ebin_width)+(min_energy);

	xe_energy_CH2[t][e] = Ebin_width/2;                     // Apparently the error in the x axis is the bin width

	// Write to file
	output->cd();
	asym_CH2_h[t][e]->Write();

      } // End of loop selecting only histograms containing events

      //---------------------------------------------------------------------------------------------------------------------------------------------------------

      //---------------------------------------------------- Obtaining G Asymmetry from Butanol data ------------------------------------------------------------
      
      // Obtain the asymmetry from the histograms if they contain points:
      
      if ((PARA_events_butanol ==0) || (PERP_events_butanol ==0)){
	tot_bin_energy_butanol[t][e] = 0.0;
	
       	cout << "Nothing in either PARA or PERP for butanol energy bin " << e << " theta bin " << t <<  " so no asymmetry will be created" << endl;
	
      }
      
      else {
	tot_bin_energy_butanol[t][e] = 1.0;
	// Scale the PARA and PERP butanol histograms using the scale factors from the butanol flux ratios

	PERP_phi_butanol_scaled_flux[t][e] = (TH1*)PERP_phi_butanol[t][e]->Clone();
	PERP_phi_butanol_scaled_flux[t][e]->Sumw2();

	PERP_phi_butanol_scaled_flux[t][e]->Scale(val_flux_ratio_butanol[t][e]);
	
	// Check that the scaling worked:
	
	cout << "Number of PERP events in butanol histogram after scaling with PARA (flux method): " << PERP_phi_butanol_scaled_flux[t][e]->Integral() << endl;

	// Extract the asymmetry:
    
        // asym_butanol_h[t][e] = PERP_phi_butanol_scaled_flux[t][e]->GetAsymmetry(PARA_phi_butanol[t][e]);
        asym_butanol_h[t][e] = PARA_phi_butanol[t][e]->GetAsymmetry(PERP_phi_butanol_scaled_flux[t][e]);
	
	// Name the asymmetry

	sprintf(asym_butanol_name,"Asymmetry_Butanol_Target_Ebin%d_Thbin%d_pol%d",e+1,t+1,at_pol);
	asym_butanol_h[t][e]->SetName(asym_butanol_name);

	// Create labels for energies and theta values

	Labels();

	// Create title for asymmetry

	sprintf(asym_butanol_title,"Asymmetry_%fGeV_to_%fGeV_%ddegrees_to_%ddegrees_for_butanol_target",energy_label[e],energy_label[e+1],theta_label[t+1],theta_label[t]);  

	asym_butanol_h[t][e]->SetTitle(asym_butanol_title);

	// Set statistics and fitting options

	gStyle->SetOptStat("e");
	gStyle->SetOptFit(1111);

	// Set the phi-offset to zero

	cos2phi_sin2phi->FixParameter(2,0);

	// Now fit the cos2phi-sin2phi function to the asymmetry:
	
	asym_butanol_h[t][e] ->Fit("cos2phi_sin2phi","QB");

	// Plot the asymmetry

	c_butanol[t][e]= new TCanvas();
	c_butanol[t][e]->cd();
	asym_butanol_h[t][e]->Draw();
	c_butanol[t][e]->Update();

	// Extract G from the fit parameters to plot later:
	
	if (peak==730 && BEAM_pol[e]!=0.0) G_pol[at_pol][t][e]=((cos2phi_sin2phi->GetParameter(3))/(dilution_factor_array[t][e]*TARG_pol*BEAM_pol[e]));
	else G_pol[at_pol][t][e]= ((cos2phi_sin2phi->GetParameter(3))/(dilution_factor_array[t][e]*TARG_pol)) ;
	// Calculate the statistical error in G.
	// Don't forget need error for dilution factor!

	G_pol_error_stat[at_pol][t][e]= TMath::Abs(G_pol[at_pol][t][e])*(TMath::Sqrt(((cos2phi_sin2phi->GetParError(3)/cos2phi_sin2phi->GetParameter(3))*(cos2phi_sin2phi->GetParError(3)/cos2phi_sin2phi->GetParameter(3)))+((TARG_pol_stat_error/TARG_pol)*(TARG_pol_stat_error/TARG_pol))+((dilution_factor_stat_error_array[t][e]/dilution_factor_array[t][e])*(dilution_factor_stat_error_array[t][e]/dilution_factor_array[t][e]))));
	printf("!!!!!!!!!!!!   G for pol %d theta bin %d energy bin %d = %.3e  +/- %.3e \n",at_pol,t,e,G_pol[at_pol][t][e],G_pol_error_stat[at_pol][t][e]);
	// Extract Sigma from the butanol histogram to compare to the sigma value from CH2

	if (peak == 730 && BEAM_pol[e]!=0.0) {
	  Sigma_butanol_pol[at_pol][t][e]= ((cos2phi_sin2phi->GetParameter(1))/BEAM_pol[e]);
	  Sigma_butanol_pol_stat_error[at_pol][t][e]= ((cos2phi_sin2phi->GetParError(1))/BEAM_pol[e]);
	}
	else {
	  Sigma_butanol_pol[at_pol][t][e]= ((cos2phi_sin2phi->GetParameter(1)));
	  Sigma_butanol_pol_stat_error[at_pol][t][e]= ((cos2phi_sin2phi->GetParError(1)));
	}
	// Get theta angle in degrees from bin number and the associated error
 
	x_angle_butanol[t][e]=theta_ax[t];
	xe_angle_butanol[t][e]=etheta_ax[t];              // Apparently the error in the xaxis is the bin width

	// Get energy in GeV from bin number and the associated error

	x_energy_butanol[t][e+at_Ebin]= (e*Ebin_width)+(min_energy)+ Ebin_width/2;
	printf("theta = %d ; energy = %d ; at energy = %d ; W= %.3f \n",t,e,e+at_Ebin,x_energy_butanol[t][e+at_Ebin]);
	xe_energy_butanol[t][e+at_Ebin] = Ebin_width/2;                     // Apparently the error in the x axis is the bin width

	// Write to file
	output->cd();
	asym_butanol_h[t][e]->Write();

      } // End of loop selecting only histograms containing events

      //---------------------------------------------------------------------------------------------------------------------------------------------------------

    } // End of loop over theta

  } // End of loop over energy

  // Now plot the asymmetries!

  // Asymmetry_plotting_theta();

  // Close all open files


  input->Close();
  output->Close();

} // End of asymmetry() function




// ----------------------- Labels() ----------------------------------------------------

void Labels(){

  // Create a label for the energy title
  
  for (int e=0; e<(no_Ebin+1);e++){
    
    energy_label[e] = (e*Ebin_width)+(min_energy);
    
  }
  
  // Create an angle label for the title
  
  for(int t=0; t< (no_Thbin+1); t++){
    
    theta_label[t] = (TMath::ACos(((Thbin_width)*t)-1.))*(180/TMath::Pi());
    //theta_label[t] = t+1;
    if (t>0) {
      theta_ax[t-1] = (theta_label[t] + theta_label[t-1]) * 0.5;
      etheta_ax[t-1] = TMath::Abs(theta_ax[t-1]-theta_label[t-1]);
      
    }
    
  }
  
}

//--------------------------------------------------------------------------------------

// ----------------------- Asymmetry_plotting() -----------------------------------------

// Plot the cos2phi_sin2phi amplitude against either theta or energy

void Asymmetry_plotting_theta(char const *outfile){
  
  // Plot the asymmetry with theta on the x axis

  TFile *output = new TFile(outfile,"RECREATE");

  int pass_Ebin = 0;
    // Create canvases on which to draw the asymmetries
    
    for (int e=0; e<no_Ebin; e++){
      
      c_G[e] = new TCanvas();
      c_Sigma[e] = new TCanvas();
      c_Sigma_butanol[e] = new TCanvas();

    }
 
    // Now plot the asymmetry along with MAID and SAID solutions
      
      for (int e=0; e<no_Ebin; e++){
	
	// Put the values of theta and cos2phi_sin2phi amplitude along with their errors into arrays ready for plotting


	for (int t=0; t<no_Thbin; t++){

	  // Define the values of G, Sigma, Sigma_butanol as weighted mean with their error

	  if (G_pol_error_stat[0][t][e] > 0.0 && G_pol_error_stat[1][t][e] > 0.0 )  {
	    G_polfull[0][t][e+at_Ebin] = G_pol[0][t][e];
	    G_polfull[1][t][e+at_Ebin] = G_pol[1][t][e];
	    G_polfull_error_stat[0][t][e+at_Ebin] = G_pol_error_stat[0][t][e];
	    G_polfull_error_stat[1][t][e+at_Ebin] = G_pol_error_stat[1][t][e];

	    G_error_stat[t][e+at_Ebin] = pow( 1./ (pow( G_pol_error_stat[0][t][e],-2) +  pow( G_pol_error_stat[1][t][e],-2)) , 0.5);
	    G[t][e+at_Ebin] =   (G_pol[0][t][e]* pow( G_pol_error_stat[0][t][e],-2) + G_pol[1][t][e]*  pow( G_pol_error_stat[1][t][e],-2)) / (pow( G_pol_error_stat[0][t][e],-2) +  pow( G_pol_error_stat[1][t][e],-2));
	  }
	  else if  (G_pol_error_stat[0][t][e] > 0.0 && G_pol_error_stat[1][t][e] == 0.0 ) {
	    G_polfull[0][t][e+at_Ebin] = G_pol[0][t][e];
	    G_polfull[1][t][e+at_Ebin] = G_pol[1][t][e];
	    G_polfull_error_stat[0][t][e+at_Ebin] = G_pol_error_stat[0][t][e];
	    G_polfull_error_stat[1][t][e+at_Ebin] = G_pol_error_stat[1][t][e];

	    G_error_stat[t][e+at_Ebin] = G_pol_error_stat[0][t][e] ;
	    G[t][e+at_Ebin] = G_pol[0][t][e];
	  }
	  else if  (G_pol_error_stat[0][t][e] == 0.0 && G_pol_error_stat[1][t][e] > 0.0 ) {
	    G_polfull[0][t][e+at_Ebin] = G_pol[0][t][e];
	    G_polfull[1][t][e+at_Ebin] = G_pol[1][t][e];
	    G_polfull_error_stat[0][t][e+at_Ebin] = G_pol_error_stat[0][t][e];
	    G_polfull_error_stat[1][t][e+at_Ebin] = G_pol_error_stat[1][t][e];

	    G_error_stat[t][e+at_Ebin] = G_pol_error_stat[1][t][e] ;
	    G[t][e+at_Ebin] = G_pol[1][t][e];
	  }
	  else {
	    G_polfull[0][t][e+at_Ebin] = G_pol[0][t][e];
	    G_polfull[1][t][e+at_Ebin] = G_pol[1][t][e];
	    G_polfull_error_stat[0][t][e+at_Ebin] = G_pol_error_stat[0][t][e];
	    G_polfull_error_stat[1][t][e+at_Ebin] = G_pol_error_stat[1][t][e];

	    G_error_stat[t][e+at_Ebin] = 0.0;
	    G[t][e+at_Ebin] = 0.0;
	    printf("No values for G for theta bin=%d and energy bin=%d \n",t,e);
	  }

	  if (Sigma_pol_error_stat[0][t][e] > 0.0 && Sigma_pol_error_stat[1][t][e] > 0.0 )  {
	    Sigma_error_stat[t][e+at_Ebin] = pow( 1./ (pow( Sigma_pol_error_stat[0][t][e],-2) +  pow( Sigma_pol_error_stat[1][t][e],-2)) , 0.5);
	    Sigma[t][e+at_Ebin] =   (Sigma_pol[0][t][e]* pow( Sigma_pol_error_stat[0][t][e],-2) + Sigma_pol[1][t][e]*  pow( Sigma_pol_error_stat[1][t][e],-2)) / (pow( Sigma_pol_error_stat[0][t][e],-2) +  pow( Sigma_pol_error_stat[1][t][e],-2));
	  }
	  else if  (Sigma_pol_error_stat[0][t][e] > 0.0 && Sigma_pol_error_stat[1][t][e] == 0.0 ) {
	    Sigma_error_stat[t][e+at_Ebin] = Sigma_pol_error_stat[0][t][e] ;
	    Sigma[t][e+at_Ebin] = Sigma_pol[0][t][e];
	  }
	  else if  (Sigma_pol_error_stat[0][t][e] == 0.0 && Sigma_pol_error_stat[1][t][e] > 0.0 ) {
	    Sigma_error_stat[t][e+at_Ebin] = Sigma_pol_error_stat[1][t][e] ;
	    Sigma[t][e+at_Ebin] = Sigma_pol[1][t][e];
	  }
	  else {
	    Sigma_error_stat[t][e+at_Ebin] = 0.0;
	    Sigma[t][e+at_Ebin] = 0.0;
	    printf("No values for Sigma for theta bin=%d and energy bin=%d \n",t,e);
	  }

	  if (Sigma_butanol_pol_stat_error[0][t][e] > 0.0 && Sigma_butanol_pol_stat_error[1][t][e] > 0.0 )  {
	    Sigma_butanol_error_stat[t][e+at_Ebin] = pow( 1./ (pow( Sigma_butanol_pol_stat_error[0][t][e],-2) +  pow( Sigma_butanol_pol_stat_error[1][t][e],-2)) , 0.5);
	    Sigma_butanol[t][e+at_Ebin] =   (Sigma_butanol_pol[0][t][e]* pow( Sigma_butanol_pol_stat_error[0][t][e],-2) + Sigma_butanol_pol[1][t][e]*  pow( Sigma_butanol_pol_stat_error[1][t][e],-2)) / (pow( Sigma_butanol_pol_stat_error[0][t][e],-2) +  pow( Sigma_butanol_pol_stat_error[1][t][e],-2));
	  }
	  else if  (Sigma_butanol_pol_stat_error[0][t][e] > 0.0 && Sigma_butanol_pol_stat_error[1][t][e] == 0.0 ) {
	    Sigma_butanol_error_stat[t][e+at_Ebin] = Sigma_butanol_pol_stat_error[0][t][e] ;
	    Sigma_butanol[t][e+at_Ebin] = Sigma_butanol_pol[0][t][e];
	  }
	  else if  (Sigma_butanol_pol_stat_error[0][t][e] == 0.0 && Sigma_butanol_pol_stat_error[1][t][e] > 0.0 ) {
	    Sigma_butanol_error_stat[t][e+at_Ebin] = Sigma_butanol_pol_stat_error[1][t][e] ;
	    Sigma_butanol[t][e+at_Ebin] = Sigma_butanol_pol[1][t][e];
	  }
	  else {
	    Sigma_butanol_error_stat[t][e+at_Ebin] = 0.0;
	    Sigma_butanol[t][e+at_Ebin] = 0.0;
	    printf("No values for Sigma_butanol for theta bin=%d and energy bin=%d \n",t,e);
	  }

	  if (tot_bin_energy_butanol[t][e] == 1) {
	    X_butanol[t] = x_angle_butanol[t][e];
	    Xe_butanol[t] = xe_angle_butanol[t][e];
	  
	    Y_butanol[t] = G[t][e+at_Ebin];
	    Ye_butanol[t] = G_error_stat[t][e+at_Ebin];

	    X_sigma_butanol[t] = x_angle_butanol[t][e];
	    Xe_sigma_butanol[t] = xe_angle_butanol[t][e];
	  
	    Y_sigma_butanol[t] = Sigma_butanol[t][e+at_Ebin];
	    Ye_sigma_butanol[t] = Sigma_butanol_error_stat[t][e+at_Ebin];
	    tot_energy_butanol_v[e+at_Ebin]++;
	    tot_bin2_energy_butanol[t][e+at_Ebin] = 1;
	  }

	  if (tot_bin_energy_ch2[t][e] == 1) {
	    X_CH2[t] = x_angle_CH2[t][e];
	    Xe_CH2[t] = xe_angle_CH2[t][e];
	  
	    Y_CH2[t] = Sigma[t][e+at_Ebin];
	    Ye_CH2[t] = Sigma_error_stat[t][e+at_Ebin];
	    tot_energy_ch2_v[e+at_Ebin]++;
	    tot_bin2_energy_ch2[t][e+at_Ebin] = 1;
	  }
	  
	}

      
      // Create a name for the graph (based on energy bins) and a title for the graph
      
      sprintf(graphname_butanol, "G_vs_theta_Ebin%d_Butanol_Target",e+at_Ebin+1);
      sprintf(graphtitle_butanol, "G vs Theta, CMS energy range %fGeV to %fGeV, Butanol Target",energy_label[e],energy_label[e+1]);

      sprintf(graphname_CH2, "Sigma_vs_theta_Ebin%d_CH2_Target",e+at_Ebin+1);
      sprintf(graphtitle_CH2, "Sigma vs Theta, CMS energy range %fGeV to %fGeV, CH2 Target",energy_label[e],energy_label[e+1]);

      sprintf(graphname_sigma_butanol, "Sigma_vs_theta_Ebin%d_Butanol_Target",e+at_Ebin+1);
      sprintf(graphtitle_sigma_butanol, "Sigma vs Theta, CMS energy range %fGeV to %fGeV,Butanol Target",energy_label[e],energy_label[e+1]);
      
      // Plot the asymmetry graph, set name, title and label x axis

      asym_graph_butanol[e+at_Ebin] = new TGraphErrors(tot_energy_butanol_v[e+at_Ebin],X_butanol,Y_butanol,Xe_butanol,Ye_butanol);
      asym_graph_butanol[e+at_Ebin]->SetName(graphname_butanol);
      asym_graph_butanol[e+at_Ebin]->SetTitle(graphtitle_butanol);
      asym_graph_butanol[e+at_Ebin]->GetXaxis()->SetTitle("Theta in CMS /degrees");
      //asym_graph_butanol[e]->GetYaxis()->SetMaximum(1);
      //asym_graph_butanol[e]->GetYaxis()->SetMinimum(-1);
      //asym_graph_butanol[e]->GetYaxis()->SetRangeUser(-1,1);
      gStyle->SetOptStat("e");
      gStyle->SetOptFit(1111);

      c_G[e]->cd();
      asym_graph_butanol[e+at_Ebin]->Draw("AP*");
      asym_graph_butanol[e+at_Ebin]->Write();
      c_G[e]->Update();

      asym_graph_CH2[e+at_Ebin] = new TGraphErrors(tot_energy_ch2_v[e+at_Ebin],X_CH2,Y_CH2,Xe_CH2,Ye_CH2);
      asym_graph_CH2[e+at_Ebin]->SetName(graphname_CH2);
      asym_graph_CH2[e+at_Ebin]->SetTitle(graphtitle_CH2);
      asym_graph_CH2[e+at_Ebin]->GetXaxis()->SetTitle("Theta in CMS /degrees");
      //asym_graph_CH2[e]->GetYaxis()->SetMaximum(1);
      //asym_graph_CH2[e]->GetYaxis()->SetMinimum(-1);
      gStyle->SetOptStat("e");
      gStyle->SetOptFit(1111);

      c_Sigma[e]->cd();
      asym_graph_CH2[e+at_Ebin]->Draw("AP*");
      asym_graph_CH2[e+at_Ebin]->Write();
      c_Sigma[e]->Update();

      asym_graph_sigma_butanol[e+at_Ebin] = new TGraphErrors(tot_energy_butanol_v[e+at_Ebin],X_sigma_butanol,Y_sigma_butanol,Xe_sigma_butanol,Ye_sigma_butanol);
      asym_graph_sigma_butanol[e+at_Ebin]->SetName(graphname_sigma_butanol);
      asym_graph_sigma_butanol[e+at_Ebin]->SetTitle(graphtitle_sigma_butanol);
      asym_graph_sigma_butanol[e+at_Ebin]->GetXaxis()->SetTitle("Theta in CMS /degrees");
      //asym_graph_sigma_butanol[e]->GetYaxis()->SetMaximum(1);
      //asym_graph_sigma_butanol[e]->GetYaxis()->SetMinimum(-1);
      gStyle->SetOptStat("e");
      gStyle->SetOptFit(1111);

      c_Sigma_butanol[e]->cd();
      asym_graph_sigma_butanol[e+at_Ebin]->Draw("AP*");
      asym_graph_sigma_butanol[e+at_Ebin]->Write();
      c_Sigma_butanol[e]->Update();

      
      }
      at_Ebin = at_Ebin + no_Ebin; 
      output->Close();
}
    

void Asymmetry_plotting_energy(char const *outfile) {    
  
  // Plot the asymmetry with energy on the x axis

  TFile *output = new TFile(outfile,"RECREATE");

    // Create an angle label for the title

    for(int t=0; t< (no_Thbin+1); t++){
      
      //theta_label[t] = (TMath::ACos(((Thbin_width)*t)-1))*(180/pi);
      theta_label[t] = t+1;
      
    }

    // Create canvases on which to draw the asymmetries
    
    for (int t=0; t<no_Thbin; t++){
      
      c_G[t] = new TCanvas();
      c_Sigma[t] = new TCanvas();
      c_Sigma_butanol[t] = new TCanvas();
      
    }
 
    // Fill arrays with energy and cos2phi_sin2phi amplitude along with their errors ready for plotting
    
    for (int t=0; t<no_Thbin; t++){
      for (int e=0; e<at_Ebin; e++){

	if (tot_bin2_energy_butanol[t][e] == 1) {
	  X_butanol[tot_angle_butanol_v[t]] = x_energy_butanol[t][e];
	  Xe_butanol[tot_angle_butanol_v[t]]= xe_energy_butanol[t][e];
	
	  Y_butanol[tot_angle_butanol_v[t]] = G[t][e];
	  Ye_butanol[tot_angle_butanol_v[t]] = G_error_stat[t][e];

	  X_sigma_butanol[tot_angle_butanol_v[t]] = x_energy_butanol[t][e];
	  Xe_sigma_butanol[tot_angle_butanol_v[t]] = xe_energy_butanol[t][e];
	
	  Y_sigma_butanol[tot_angle_butanol_v[t]] = Sigma_butanol[t][e];
	  Ye_sigma_butanol[tot_angle_butanol_v[t]] = Sigma_butanol_error_stat[t][e];
	  tot_angle_butanol_v[t]++;
	  printf("Butanol data for theta=%d ; energy=%d \n",t,e);
	}
	

	if (tot_bin2_energy_ch2[t][e] == 1) {
	  X_CH2[tot_angle_ch2_v[t]] = x_energy_CH2[t][e]+(Ebin_width/2);
	  Xe_CH2[tot_angle_ch2_v[t]]= xe_energy_CH2[t][e];
	
	  Y_CH2[tot_angle_ch2_v[t]] = Sigma[t][e];
	  Ye_CH2[tot_angle_ch2_v[t]] = TMath::Abs(Sigma_error_stat[t][e]);
	  tot_angle_ch2_v[t]++;
	  printf("CH2 data for theta=%d ; energy=%d \n",t,e);

	}
	printf("theta = %d ; energy= %d ; passed = %d ; tot butanol=%d ; tot CH2=%d ; X = %.3f  ; G=%.3f  ; S = %.3f \n",t,e,tot_bin2_energy_butanol[t][e],tot_angle_butanol_v[t],tot_angle_ch2_v[t],x_energy_butanol[t][e],G[t][e],Sigma[t][e]);
	  
      }
      
      // Create graph name and title
      
      sprintf(graphname_butanol, "G_vs_energy_Thbin%d_Butanol_Target",t+1);
      sprintf(graphtitle_butanol, "G vs Energy, Theta Bin %f, Butanol Target",theta_label[t]);


      sprintf(graphname_CH2, "Sigma_vs_energy_Thbin%d_CH2_Target",t+1);
      sprintf(graphtitle_CH2, "Sigma vs Energy, Theta Bin %f, CH2 Target",theta_label[t]);

      sprintf(graphname_sigma_butanol, "Sigma_vs_energy_Thbin%d_Butanol_Target",t+1);
      sprintf(graphtitle_sigma_butanol, "Sigma vs Energy, Theta Bin %f, Butanol Target",theta_label[t]);
      
      // Plot the graph, set title, name and label x axis

      asym_graph_butanol[t] = new TGraphErrors(tot_angle_butanol_v[t],X_butanol,Y_butanol,Xe_butanol,Ye_butanol);
      asym_graph_butanol[t]->SetName(graphname_butanol);
      asym_graph_butanol[t]->SetTitle(graphtitle_butanol);
      asym_graph_butanol[t]->GetXaxis()->SetTitle("Energy in CMS /GeV");
      //asym_graph_butanol[e]->GetYaxis()->SetMaximum(1);
      //asym_graph_butanol[e]->GetYaxis()->SetMinimum(-1);
      gStyle->SetOptStat("n");
      gStyle->SetOptFit(0011);

      c_G[t]->cd();
      asym_graph_butanol[t]->Draw("AP*");
      asym_graph_butanol[t]->Write();
      c_G[t]->Update();

      sprintf(graphname_butanol_pol, "G_vs_energy_Thbin%d_Butanol_Target_pol0",t+1);
      sprintf(graphtitle_butanol_pol, "G vs Energy, Theta Bin %f, Butanol Target POS(red) vs NEG(black)",theta_label[t]);

      asym_graph_butanol_pol[0][t] = new TGraphErrors(tot_angle_butanol_v[t],X_butanol,G_polfull[0][t],Xe_butanol,G_polfull_error_stat[0][t]);
      asym_graph_butanol_pol[0][t]->SetName(graphname_butanol_pol);
      asym_graph_butanol_pol[0][t]->SetTitle(graphtitle_butanol_pol);
      sprintf(graphname_butanol_pol, "G_vs_energy_Thbin%d_Butanol_Target_pol1",t+1);
      asym_graph_butanol_pol[1][t] = new TGraphErrors(tot_angle_butanol_v[t],X_butanol,G_polfull[1][t],Xe_butanol,G_polfull_error_stat[1][t]);
      asym_graph_butanol_pol[1][t]->SetName(graphname_butanol_pol); 
      asym_graph_butanol_pol[1][t]->SetTitle(graphtitle_butanol_pol);     
      asym_graph_butanol_pol[1][t]->SetLineColor(2);
      // Plot the graph, set title, name and label x axis
      c_G[t]->cd();
      asym_graph_butanol_pol[0][t]->Draw("AP");
      asym_graph_butanol_pol[1][t]->Draw("SAMEP");
      asym_graph_butanol_pol[0][t]->Write();
      asym_graph_butanol_pol[1][t]->Write();
      c_G[t]->Update();

      asym_graph_CH2[t] = new TGraphErrors(tot_angle_ch2_v[t],X_CH2,Y_CH2,Xe_CH2,Ye_CH2);
      asym_graph_CH2[t]->SetName(graphname_CH2);
      asym_graph_CH2[t]->SetTitle(graphtitle_CH2);
      asym_graph_CH2[t]->GetXaxis()->SetTitle("Energy in CMS /GeV");
      //asym_graph_CH2[e]->GetYaxis()->SetMaximum(1);
      //asym_graph_CH2[e]->GetYaxis()->SetMinimum(-1);
      gStyle->SetOptStat("n");
      gStyle->SetOptFit(0011);

      c_Sigma[t]->cd();
      asym_graph_CH2[t]->Draw("AP*");
      asym_graph_CH2[t]->Write();
      c_Sigma[t]->Update();

      // Plot the graph, set title, name and label x axis

      asym_graph_sigma_butanol[t] = new TGraphErrors(tot_angle_butanol_v[t],X_sigma_butanol,Y_sigma_butanol,Xe_sigma_butanol,Ye_sigma_butanol);
      asym_graph_sigma_butanol[t]->SetName(graphname_sigma_butanol);
      asym_graph_sigma_butanol[t]->SetTitle(graphtitle_sigma_butanol);
      asym_graph_sigma_butanol[t]->GetXaxis()->SetTitle("Energy in CMS /GeV");
      // asym_graph_sigma_butanol[e]->GetYaxis()->SetRangeUser(-1,1);
      gStyle->SetOptStat("n");
      gStyle->SetOptFit(0011);

      c_Sigma_butanol[t]->cd();
      asym_graph_sigma_butanol[t]->Draw("AP*");
      asym_graph_sigma_butanol[t]->Write();
      c_Sigma_butanol[t]->Update();

      //asym_graph[e]->Write();
      
    }

    output->Close();
    
}

// -----------------------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------



void Run_all() {
  asymmetry("output_730_pos.root","histo_730_pos.root",730,1);
  asymmetry("output_730_neg.root","histo_730_neg.root",730,-1);
  Asymmetry_plotting_theta("asym_theta_730.root");
  asymmetry("output_930_pos.root","histo_930_pos.root",930,1);
  asymmetry("output_930_neg.root","histo_930_neg.root",930,-1);
  Asymmetry_plotting_theta("asym_theta_930.root");
  asymmetry("output_1100_pos.root","histo_1100_pos.root",1100,1);
  //  asymmetry("output_1100_neg.root","histo_1100_neg.root",1100,-1);
  Asymmetry_plotting_theta("asym_theta_1100.root");
  asymmetry("output_1300_pos.root","histo_1300_pos.root",1300,1);
  asymmetry("output_1300_neg.root","histo_1300_neg.root",1300,-1);
  Asymmetry_plotting_theta("asym_theta_1300.root");
  asymmetry("output_1500_pos.root","histo_1500_pos.root",1500,1);
  asymmetry("output_1500_neg.root","histo_1500_neg.root",1500,-1);
  Asymmetry_plotting_theta("asym_theta_1500.root");
  asymmetry("output_1700_pos.root","histo_1700_pos.root",1700,1);
  asymmetry("output_1700_neg.root","histo_1700_neg.root",1700,-1);
  Asymmetry_plotting_theta("asym_theta_1700.root");
  asymmetry("output_1900_pos.root","histo_1900_pos.root",1900,1);
  asymmetry("output_1900_neg.root","histo_1900_neg.root",1900,-1);
  Asymmetry_plotting_theta("asym_theta_1900.root");
  asymmetry("output_2100_pos.root","histo_2100_pos.root",2100,1);
  asymmetry("output_2100_neg.root","histo_2100_neg.root",2100,-1);
  Asymmetry_plotting_theta("asym_theta_2100.root");
  asymmetry("output_2300_pos.root","histo_2300_pos.root",2300,1);
  asymmetry("output_2300_neg.root","histo_2300_neg.root",2300,-1);
  Asymmetry_plotting_theta("asym_theta_2300.root");
  Asymmetry_plotting_energy("asym_energy_all.root");

}


void Run_all_pol05() {
  asymmetry("output_pol05_730_pos.root","histo_pol05_730_pos.root",730,1);
  asymmetry("output_pol05_730_neg.root","histo_pol05_730_neg.root",730,-1);
  Asymmetry_plotting_theta("asym_pol05_theta_730.root");
  asymmetry("output_pol05_930_pos.root","histo_pol05_930_pos.root",930,1);
  asymmetry("output_pol05_930_neg.root","histo_pol05_930_neg.root",930,-1);
  Asymmetry_plotting_theta("asym_pol05_theta_930.root");
  asymmetry("output_pol05_1100_pos.root","histo_pol05_1100_pos.root",1100,1);
  asymmetry("output_pol05_1100_neg.root","histo_pol05_1100_neg.root",1100,-1);
  Asymmetry_plotting_theta("asym_pol05_theta_1100.root");
  asymmetry("output_pol05_1300_pos.root","histo_pol05_1300_pos.root",1300,1);
  asymmetry("output_pol05_1300_neg.root","histo_pol05_1300_neg.root",1300,-1);
  Asymmetry_plotting_theta("asym_pol05_theta_1300.root");
  asymmetry("output_pol05_1500_pos.root","histo_pol05_1500_pos.root",1500,1);
  asymmetry("output_pol05_1500_neg.root","histo_pol05_1500_neg.root",1500,-1);
  Asymmetry_plotting_theta("asym_pol05_theta_1500.root");
  asymmetry("output_pol05_1700_pos.root","histo_pol05_1700_pos.root",1700,1);
  asymmetry("output_pol05_1700_neg.root","histo_pol05_1700_neg.root",1700,-1);
  Asymmetry_plotting_theta("asym_pol05_theta_1700.root");
  asymmetry("output_pol05_1900_pos.root","histo_pol05_1900_pos.root",1900,1);
  asymmetry("output_pol05_1900_neg.root","histo_pol05_1900_neg.root",1900,-1);
  Asymmetry_plotting_theta("asym_pol05_theta_1900.root");
  asymmetry("output_pol05_2100_pos.root","histo_pol05_2100_pos.root",2100,1);
  asymmetry("output_pol05_2100_neg.root","histo_pol05_2100_neg.root",2100,-1);
  Asymmetry_plotting_theta("asym_pol05_theta_2100.root");
  asymmetry("output_pol05_2300_pos.root","histo_pol05_2300_pos.root",2300,1);
  asymmetry("output_pol05_2300_neg.root","histo_pol05_2300_neg.root",2300,-1);
  Asymmetry_plotting_theta("asym_pol05_theta_2300.root");
  Asymmetry_plotting_energy("asym_pol05_energy_all.root");

}

void Run_all_pol055() {
  asymmetry("output_pol055_730_pos.root","histo_pol055_730_pos.root",730,1);
  asymmetry("output_pol055_730_neg.root","histo_pol055_730_neg.root",730,-1);
  Asymmetry_plotting_theta("asym_pol055_theta_730.root");
  asymmetry("output_pol055_930_pos.root","histo_pol055_930_pos.root",930,1);
  asymmetry("output_pol055_930_neg.root","histo_pol055_930_neg.root",930,-1);
  Asymmetry_plotting_theta("asym_pol055_theta_930.root");
  asymmetry("output_pol055_1100_pos.root","histo_pol055_1100_pos.root",1100,1);
  asymmetry("output_pol055_1100_neg.root","histo_pol055_1100_neg.root",1100,-1);
  Asymmetry_plotting_theta("asym_pol055_theta_1100.root");
  asymmetry("output_pol055_1300_pos.root","histo_pol055_1300_pos.root",1300,1);
  asymmetry("output_pol055_1300_neg.root","histo_pol055_1300_neg.root",1300,-1);
  Asymmetry_plotting_theta("asym_pol055_theta_1300.root");
  asymmetry("output_pol055_1500_pos.root","histo_pol055_1500_pos.root",1500,1);
  asymmetry("output_pol055_1500_neg.root","histo_pol055_1500_neg.root",1500,-1);
  Asymmetry_plotting_theta("asym_pol055_theta_1500.root");
  asymmetry("output_pol055_1700_pos.root","histo_pol055_1700_pos.root",1700,1);
  asymmetry("output_pol055_1700_neg.root","histo_pol055_1700_neg.root",1700,-1);
  Asymmetry_plotting_theta("asym_pol055_theta_1700.root");
  asymmetry("output_pol055_1900_pos.root","histo_pol055_1900_pos.root",1900,1);
  asymmetry("output_pol055_1900_neg.root","histo_pol055_1900_neg.root",1900,-1);
  Asymmetry_plotting_theta("asym_pol055_theta_1900.root");
  asymmetry("output_pol055_2100_pos.root","histo_pol055_2100_pos.root",2100,1);
  asymmetry("output_pol055_2100_neg.root","histo_pol055_2100_neg.root",2100,-1);
  Asymmetry_plotting_theta("asym_pol055_theta_2100.root");
  asymmetry("output_pol055_2300_pos.root","histo_pol055_2300_pos.root",2300,1);
  asymmetry("output_pol055_2300_neg.root","histo_pol055_2300_neg.root",2300,-1);
  Asymmetry_plotting_theta("asym_pol055_theta_2300.root");
  Asymmetry_plotting_energy("asym_pol055_energy_all.root");

}

void Run_all_pol045() {
  asymmetry("output_pol045_730_pos.root","histo_pol045_730_pos.root",730,1);
  asymmetry("output_pol045_730_neg.root","histo_pol045_730_neg.root",730,-1);
  Asymmetry_plotting_theta("asym_pol045_theta_730.root");
  asymmetry("output_pol045_930_pos.root","histo_pol045_930_pos.root",930,1);
  asymmetry("output_pol045_930_neg.root","histo_pol045_930_neg.root",930,-1);
  Asymmetry_plotting_theta("asym_pol045_theta_930.root");
  asymmetry("output_pol045_1100_pos.root","histo_pol045_1100_pos.root",1100,1);
  asymmetry("output_pol045_1100_neg.root","histo_pol045_1100_neg.root",1100,-1);
  Asymmetry_plotting_theta("asym_pol045_theta_1100.root");
  asymmetry("output_pol045_1300_pos.root","histo_pol045_1300_pos.root",1300,1);
  asymmetry("output_pol045_1300_neg.root","histo_pol045_1300_neg.root",1300,-1);
  Asymmetry_plotting_theta("asym_pol045_theta_1300.root");
  asymmetry("output_pol045_1500_pos.root","histo_pol045_1500_pos.root",1500,1);
  asymmetry("output_pol045_1500_neg.root","histo_pol045_1500_neg.root",1500,-1);
  Asymmetry_plotting_theta("asym_pol045_theta_1500.root");
  asymmetry("output_pol045_1700_pos.root","histo_pol045_1700_pos.root",1700,1);
  asymmetry("output_pol045_1700_neg.root","histo_pol045_1700_neg.root",1700,-1);
  Asymmetry_plotting_theta("asym_pol045_theta_1700.root");
  asymmetry("output_pol045_1900_pos.root","histo_pol045_1900_pos.root",1900,1);
  asymmetry("output_pol045_1900_neg.root","histo_pol045_1900_neg.root",1900,-1);
  Asymmetry_plotting_theta("asym_pol045_theta_1900.root");
  asymmetry("output_pol045_2100_pos.root","histo_pol045_2100_pos.root",2100,1);
  asymmetry("output_pol045_2100_neg.root","histo_pol045_2100_neg.root",2100,-1);
  Asymmetry_plotting_theta("asym_pol045_theta_2100.root");
  asymmetry("output_pol045_2300_pos.root","histo_pol045_2300_pos.root",2300,1);
  asymmetry("output_pol045_2300_neg.root","histo_pol045_2300_neg.root",2300,-1);
  Asymmetry_plotting_theta("asym_pol045_theta_2300.root");
  Asymmetry_plotting_energy("asym_pol045_energy_all.root");

}


void Run_all_pi0() {
  asymmetry_pi0("output_pi0_730_pos.root","histo_pi0_730_pos.root",730,1);
  asymmetry_pi0("output_pi0_730_neg.root","histo_pi0_730_neg.root",730,-1);
  Asymmetry_plotting_theta("asym_pi0_theta_730.root");
  asymmetry_pi0("output_pi0_930_pos.root","histo_pi0_930_pos.root",930,1);
  asymmetry_pi0("output_pi0_930_neg.root","histo_pi0_930_neg.root",930,-1);
  Asymmetry_plotting_theta("asym_pi0_theta_930.root");
  asymmetry_pi0("output_pi0_1100_pos.root","histo_pi0_1100_pos.root",1100,1);
  asymmetry_pi0("output_pi0_1100_neg.root","histo_pi0_1100_neg.root",1100,-1);
  Asymmetry_plotting_theta("asym_pi0_theta_1100.root");
  asymmetry_pi0("output_pi0_1300_pos.root","histo_pi0_1300_pos.root",1300,1);
  asymmetry_pi0("output_pi0_1300_neg.root","histo_pi0_1300_neg.root",1300,-1);
  Asymmetry_plotting_theta("asym_pi0_theta_1300.root");
  asymmetry_pi0("output_pi0_1500_pos.root","histo_pi0_1500_pos.root",1500,1);
  asymmetry_pi0("output_pi0_1500_neg.root","histo_pi0_1500_neg.root",1500,-1);
  Asymmetry_plotting_theta("asym_pi0_theta_1500.root");
  asymmetry_pi0("output_pi0_1700_pos.root","histo_pi0_1700_pos.root",1700,1);
  asymmetry_pi0("output_pi0_1700_neg.root","histo_pi0_1700_neg.root",1700,-1);
  Asymmetry_plotting_theta("asym_pi0_theta_1700.root");
  asymmetry_pi0("output_pi0_1900_pos.root","histo_pi0_1900_pos.root",1900,1);
  asymmetry_pi0("output_pi0_1900_neg.root","histo_pi0_1900_neg.root",1900,-1);
  Asymmetry_plotting_theta("asym_pi0_theta_1900.root");
  asymmetry_pi0("output_pi0_2100_pos.root","histo_pi0_2100_pos.root",2100,1);
  asymmetry_pi0("output_pi0_2100_neg.root","histo_pi0_2100_neg.root",2100,-1);
  Asymmetry_plotting_theta("asym_pi0_theta_2100.root");
  asymmetry_pi0("output_pi0_2300_pos.root","histo_pi0_2300_pos.root",2300,1);
  asymmetry_pi0("output_pi0_2300_neg.root","histo_pi0_2300_neg.root",2300,-1);
  Asymmetry_plotting_theta("asym_pi0_theta_2300.root");
  Asymmetry_plotting_energy("asym_pi0_energy_all.root");

}



void get_flux(char const *infile){
  
  TFile *input_flux = new TFile(infile);
  


  // Loop through all the PARA, PERP and AMO histograms produced by the analysis code
  
  for(int e=0; e<no_Ebin; e++){
    
    for(int t=0; t<no_Thbin; t++){
      
      // Get name of histogram and save to PARA_name:
      
      sprintf(PARA_CH2_name,"phi_pip_PARA_Target3_E%d_Th%d",e+1,t+1);
      sprintf(PARA_butanol_name,"phi_pip_PARA_Target1_E%d_Th%d",e+1,t+1);
      
      // Find histogram:
      
      PARA_phi_CH2[t][e]=(TH1*)input_flux->Get(PARA_CH2_name);
      PARA_phi_butanol[t][e]=(TH1*)input_flux->Get(PARA_butanol_name);
      
      // Rebin for easy asymmetry viewing:
      
      PARA_phi_CH2[t][e]->Rebin(5);
      PARA_phi_butanol[t][e]->Rebin(5);
      
      // For error calculation later:
      
      PARA_phi_CH2[t][e]->Sumw2();
      PARA_phi_butanol[t][e]->Sumw2();
      
      // Find the number of events in the histogram:
      
      PARA_events_CH2 = PARA_phi_CH2[t][e]->Integral();
      PARA_events_butanol = PARA_phi_butanol[t][e]->Integral();
      
      cout << "Energy: " << e << ", Theta: " << t << ", Number of PARA events for CH2: " << PARA_events_CH2 << endl;
      cout << "Energy: " << e << ", Theta: " << t << ", Number of PARA events for butanol: " << PARA_events_butanol << endl;
      
      // Repeat the above for the PERP histogram:
      
      sprintf(PERP_CH2_name,"phi_pip_PERP_Target3_E%d_Th%d",e+1,t+1);
      sprintf(PERP_butanol_name,"phi_pip_PERP_Target1_E%d_Th%d",e+1,t+1);
      
      PERP_phi_CH2[t][e]=(TH1*)input_flux->Get(PERP_CH2_name);
      PERP_phi_butanol[t][e]=(TH1*)input_flux->Get(PERP_butanol_name);
      
      PERP_phi_CH2[t][e]->Rebin(5);
      PERP_phi_butanol[t][e]->Rebin(5);
      
      PERP_phi_CH2[t][e]->Sumw2();
      PERP_phi_butanol[t][e]->Sumw2();
      
      PERP_events_CH2 = PERP_phi_CH2[t][e]->Integral();
      PERP_events_butanol = PERP_phi_butanol[t][e]->Integral();
      
      cout << "Energy: " << e << ", Theta: " << t << ", Number of PERP events for CH2: " << PERP_events_CH2 << endl;
      cout << "Energy: " << e << ", Theta: " << t << ", Number of PERP events for butanol: " << PERP_events_butanol << endl;
      
      // Repeat the above with the AMO histogram:
      
      sprintf(AMO_CH2_name,"phi_pip_AMO_Target3_E%d_Th%d",e+1,t+1);
      sprintf(AMO_butanol_name,"phi_pip_AMO_Target1_E%d_Th%d",e+1,t+1);
      
      AMO_phi_CH2[t][e]=(TH1*)input_flux->Get(AMO_CH2_name);
      AMO_phi_butanol[t][e]=(TH1*)input_flux->Get(AMO_butanol_name);
      
      AMO_phi_CH2[t][e]->Rebin(5);
      AMO_phi_butanol[t][e]->Rebin(5);
      
      AMO_phi_CH2[t][e]->Sumw2();
      AMO_phi_butanol[t][e]->Sumw2();
      
      AMO_events_CH2 = AMO_phi_CH2[t][e]->Integral();
      AMO_events_butanol = AMO_phi_butanol[t][e]->Integral();
      
      cout << "Energy: " << e << ", Theta: " << t << ", Number of AMO events for CH2: " << AMO_events_CH2 << endl;
      cout << "Energy: " << e << ", Theta: " << t << ", Number of AMO events for butanol: " << AMO_events_butanol << endl;
      
      // Put the number of events into arrays
      
      events_CH2[t][e][0] = PARA_events_CH2;
      events_CH2[t][e][1] = PERP_events_CH2;
      events_CH2[t][e][2] = PARA_events_CH2 + PERP_events_CH2;
      events_CH2[t][e][3] = AMO_events_CH2;
      
      events_butanol[t][e][0] = PARA_events_butanol;
      events_butanol[t][e][1] = PERP_events_butanol;
      events_butanol[t][e][2] = PARA_events_butanol + PERP_events_butanol;
      events_butanol[t][e][3] = AMO_events_butanol;
      
      //---------------------------------------------------- Obtaining CH2 flux ratio ------------------------------------------------------------
      
      // Obtain the asymmetry from the histograms if they contain points:
      
      if ((PARA_events_CH2 ==0) || (PERP_events_CH2 ==0) || (AMO_events_CH2 ==0)){
	
	cout << "Nothing in either PARA or PERP for CH2 energy bin " << e << " theta bin " << t <<  " so no asymmetry will be created" << endl;
	
      }
      
      else {
	
	cout << "energy bin: " << e << " theta bin: " << t << endl;
	
	// Divide the PARA histogram by the AMO histogram to remove acceptance effects
	
	PARA_asym_CH2[t][e] = (TH1*)PARA_phi_CH2[t][e]->Clone();
	PARA_asym_CH2[t][e]->Sumw2();
	PARA_asym_CH2[t][e]->Divide(AMO_phi_CH2[t][e]);
	
	sprintf(PARA_asym_CH2_name,"Normalised_PARA_histogram_CH2_Ebin%d_Thbin%d",e+1,t+1);
	
	PARA_asym_CH2[t][e]->SetName(PARA_asym_CH2_name);

	// Set the phi-offset to zero

	PARA_fit_function_CH2->FixParameter(2,0);
	
	//Fit the normalised asymmetry to obtain phi offset and the flux ratio
	
	PARA_asym_CH2[t][e]->Fit("PARA_fit_function_CH2","B");
	
	// Now write this histogram to file
	
	
	// Extract the flux ratio F_PARA/F_AMO
	
	PARA_flux_ratio_CH2 = PARA_fit_function_CH2->GetParameter(0);
	PARA_flux_ratio_CH2_error = PARA_fit_function_CH2->GetParError(0);
	
	cout << "CH2 PARA Flux Ratio: " << PARA_flux_ratio_CH2 << " +/- " << PARA_flux_ratio_CH2_error << endl;
	
	// Divide the PERP histogram by the AMO histogram to remove acceptance effects
	
	PERP_asym_CH2[t][e] = (TH1*)PERP_phi_CH2[t][e]->Clone();
	PERP_asym_CH2[t][e]->Sumw2();
	PERP_asym_CH2[t][e]->Divide(AMO_phi_CH2[t][e]);
	
	sprintf(PERP_asym_CH2_name,"Normalised_PERP_histogram_CH2_Ebin%d_Thbin%d",e+1,t+1);

	PERP_asym_CH2[t][e]->SetName(PERP_asym_CH2_name);

	// Set the phi-offset to zero

	PERP_fit_function_CH2->FixParameter(2,0);
	
	// Fit the normalised asymmetry to obtain phi offset

	PERP_asym_CH2[t][e]->Fit("PERP_fit_function_CH2","B");

	// Now write this histogram to file

	// Extract the flux ratio F_PERP/F_AMO
	
	PERP_flux_ratio_CH2 = PERP_fit_function_CH2->GetParameter(0);
	PERP_flux_ratio_CH2_error = PERP_fit_function_CH2->GetParError(0);
	
	cout << "CH2 PERP Flux Ratio: " << PERP_flux_ratio_CH2 << " +/- " << PERP_flux_ratio_CH2_error << endl;

	// Find the flux ratio F_PARA/F_PERP

	val_flux_ratio_CH2[t][e]=0.;
	val_flux_ratio_error_CH2[t][e]=0.;

	if (  PERP_flux_ratio_CH2 != 0.0 ) {
	  val_flux_ratio_CH2[t][e]=PARA_flux_ratio_CH2/PERP_flux_ratio_CH2;
	  val_flux_ratio_error_CH2[t][e] = TMath::Sqrt(((PARA_flux_ratio_CH2_error/PARA_flux_ratio_CH2)*(PARA_flux_ratio_CH2_error/PARA_flux_ratio_CH2))+((PERP_flux_ratio_CH2_error/PERP_flux_ratio_CH2)*(PERP_flux_ratio_CH2_error/PERP_flux_ratio_CH2)));
	}
	if (val_flux_ratio_CH2[t][e] == 0.0) val_flux_ratio_CH2[t][e] =1; // put to 1 the flux normalization in case there are no counts in the PARA_flux_ratio_CH2
	cout << "Final CH2 flux ratio: " << val_flux_ratio_CH2[t][e] << "   +/-   " << val_flux_ratio_error_CH2[t][e] << endl; 
	
	
      } // End of loop selecting only histograms containing events
      
      //---------------------------------------------------------------------------------------------------------------------------------------------------------
      
      //---------------------------------------------------- Obtaining flux ratio from Butanol data ------------------------------------------------------------
      
      // Obtain the asymmetry from the histograms if they contain points:
      
      if ((PARA_events_butanol ==0) || (PERP_events_butanol ==0) || (AMO_events_butanol ==0)){
	
       	cout << "Nothing in either PARA or PERP for butanol energy bin " << e << " theta bin " << t <<  " so no asymmetry will be created" << endl;
	
      }
      
      else {
	
	//--------------------------------------------------------------------------
	
	// Divide the PARA histogram by the AMO histogram to remove acceptance effects

	PARA_asym_butanol[t][e] = (TH1*)PARA_phi_butanol[t][e]->Clone();
	PARA_asym_butanol[t][e]->Sumw2();
	PARA_asym_butanol[t][e]->Divide(AMO_phi_butanol[t][e]);

	sprintf(PARA_asym_butanol_name,"Normalised_PARA_histogram_butanol_Ebin%d_Thbin%d",e+1,t+1);

	PARA_asym_butanol[t][e]->SetName(PARA_asym_butanol_name);

	// Set the phi-offset to zero

	PARA_fit_function_butanol->FixParameter(2,0);
	
	//Fit the normalised asymmetry to obtain phi offset and the flux ratio

	PARA_asym_butanol[t][e]->Fit("PARA_fit_function_butanol","B");

	// Now write this histogram to file


	// Extract the flux ratio F_PARA/F_AMO

	PARA_flux_ratio_butanol = PARA_fit_function_butanol->GetParameter(0);
	PARA_flux_ratio_butanol_error = PARA_fit_function_butanol->GetParError(0);
	
	cout << "Butanol PARA Flux Ratio: " << PARA_flux_ratio_butanol << " +/- " << PARA_flux_ratio_butanol_error << endl;
	
	// Divide the PERP histogram by the AMO histogram to remove acceptance effects
	
	PERP_asym_butanol[t][e] = (TH1*)PERP_phi_butanol[t][e]->Clone();
	PERP_asym_butanol[t][e]->Sumw2();
	PERP_asym_butanol[t][e]->Divide(AMO_phi_butanol[t][e]);
	
	sprintf(PERP_asym_butanol_name,"Normalised_PERP_histogram_butanol_Ebin%d_Thbin%d",e+1,t+1);
	
	PERP_asym_butanol[t][e]->SetName(PERP_asym_butanol_name);

	// Set the phi-offset to zero

	PERP_fit_function_butanol->FixParameter(2,0);
	
	// Fit the normalised asymmetry to obtain phi offset
	
	PERP_asym_butanol[t][e]->Fit("PERP_fit_function_butanol","B");
	
	// Now write this histogram to file
	
	
	// Extract the flux ratio F_PERP/F_AMO
	
	PERP_flux_ratio_butanol = PERP_fit_function_butanol->GetParameter(0);
	PERP_flux_ratio_butanol_error = PERP_fit_function_butanol->GetParError(0);
	
	cout << "Butanol PERP Flux Ratio: " << PERP_flux_ratio_butanol << " +/- " << PERP_flux_ratio_butanol_error << endl;

	// Extract the flux ratio F_PERP/F_AMO
	
	PERP_flux_ratio_butanol = PERP_fit_function_butanol->GetParameter(0);
	PERP_flux_ratio_butanol_error = PERP_fit_function_butanol->GetParError(0);
	
	cout << "butanol PERP Flux Ratio: " << PERP_flux_ratio_butanol << " +/- " << PERP_flux_ratio_butanol_error << endl;

	// Find the flux ratio F_PARA/F_PERP

	val_flux_ratio_butanol[t][e]=0.;
	val_flux_ratio_error_butanol[t][e]=0.;
	if ( PERP_flux_ratio_butanol != 0.0) {
	  val_flux_ratio_butanol[t][e]=PARA_flux_ratio_butanol/PERP_flux_ratio_butanol;
	  val_flux_ratio_error_butanol[t][e] = TMath::Sqrt(((PARA_flux_ratio_butanol_error/PARA_flux_ratio_butanol)*(PARA_flux_ratio_butanol_error/PARA_flux_ratio_butanol))+((PERP_flux_ratio_butanol_error/PERP_flux_ratio_butanol)*(PERP_flux_ratio_butanol_error/PERP_flux_ratio_butanol)));
	}
	cout << "Final butanol flux ratio: " << val_flux_ratio_butanol[t][e] << "   +/-   " << val_flux_ratio_error_butanol[t][e] << endl; 

		
      } // End of loop selecting only histograms containing events
      
      //     counter++;

      //---------------------------------------------------------------------------------------------------------------------------------------------------------
      
    } // End of loop over theta
    
  } // End of loop over energy
  
  input_flux->Close();
  
} // End of get_flux() function
