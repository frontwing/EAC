#ifndef PARAMETER_H_
#define PARAMETER_H_


//particle number
const int N_p_initial   = 20000;
const int N_p_diagnosis = 1;

//tokamak parameter
const double a_b   = 1.3;
const double R0    = 1.0;
const double a     = 0.33333333;
const double B0    = 1.0;
/*************** for n=1 TAE case ******************************/
const double q_0   = 1.1;
/*************** for n=3 TAE case ******************************/
//const double q_0   = 1.25;
const double q_a   = q_0 + 1;
double psi_0 = -B0*a*a/2/(q_0+0.5);
double psi_a = 0.0;
double Delta_psi = psi_a - psi_0;

//grid parameter

const int mR    = 128;
const int mZ    = 128;
const int mphi  = 36;
const double dR = 2*(a*a_b)/(mR-1);
const double dZ = 2*(a*a_b)/(mZ-1);


//alpha = omega_A/omega_c
const double alpha = 1;

//how many modes in simulation
const int NUM_MODE = 1;
/*************** for n=1 TAE case ******************************/
double n[NUM_MODE];
double omega_A[NUM_MODE];

//alpha = omega_c^2/omega_A^2
const double alpha1      = 3600;
const double v_A         = 1.0/sqrt(alpha1);
const double omega_Alfen = 1.0/sqrt(alpha1);
//slowing down distrubution parameter
double v_0=1.7*v_A,v_c = 0.58*v_0,c_0=1.0,c_1=0.37,v_th=0.1;
//double omega_A =  -0.3266*omega_Alfen;
//double omega_A =  -0.3117*omega_Alfen;

/****************************************************************/


/*************** for n=3 TAE case ******************************
double n = 3;
//alpha = omega_c^2/omega_A^2
//const double alpha1      = 32859.56872978810701169948432369;
const double alpha1      = 32860;
const double v_A         = 1.0/sqrt(alpha1);
const double omega_Alfen = 1.0/sqrt(alpha1);
//slowing down distrubution parameter
double v_0=1.71*v_A,v_c = 0.5848*v_0,c_0=1.0,c_1=0.25,v_th=0.1;
//double omega_A = -0.32062*omega_Alfen;
double omega_A = -0.3132*omega_Alfen;
/****************************************************************/


//Calculate particle Orbit Frequency (true) or not(false)
bool COF = false;
//Orbit in Phase Spase : 
//1 : (E,P_phi)
//2 : (Lambda,E)
//3 : (Lambda,P_phi)
double OPS =2;
/*****parameter for calculation of particle orbit frequency******/
int    NUM_PPHI             = -1;
int    NUM_E                = 100;
int    NUM_LAMBDA           = 100;
double LAMBDA               = 0.80;
double P_phi_0              = -0.045;
int    NUM_TRACK_PER_PPHI   = -1;
int    NUM_TRACK_PER_E      = 10;
int    NUM_TRACK_PER_LAMBDA = 10;
bool   OUTPUT_TRACK;
bool   LOSS;
bool   CO_PASSING = false;//true for co-passing while false for counter-passing particles
//parameter for (Lambda,E)
double LAMBDA_min  = 0.0;
double LAMBDA_max  = 1.2;
/****************************************************************/

/***************Calculate Phase including Nonlinear terms (true) or not(false)****************/
//Phase : n*phi + p*theta + omega
//for resonances, dot{Phase} = 0 e.g. n*omega_phi + p*omega_theta + t = 0
//In classical, omega_phi and omega_theta are calculated on unperturbed particle orbit
//If including nonlinear terms, dPhi/dr and dA_||/dr terms will be added
//and frequnecy of phase is also be calulated in E-P_phi phase space : dot{alpha}
bool CPN = false;
double t_start_CPN = 600;//unit: Alfen time, time started to calculate phase Delta_alpha
int index_t_start_CPN;
int NUM_LINE_AMP;//total lines of file amplitude.dat

/***********************************************************/

/***************output perturbed_field_vs_time at (R_field,Z_field,phi_field)*****************/
double R_field    = R0 + 0.5*a;
double Z_field    = 0.0*a;
double phi_field  = 0.0;



//Evolution of Distribution function as function of Time (true) or not(false)
bool EDT = false;
//output time evolution of Distribution function as function of P_phi and Energy (true) or not(false)
bool DPE = false;
//output time evolution of Distribution function as function of P_phi and K(K = ~P_phi - n/omega_A*~E) (true) or not(false)
bool DPK = false;
//output f and deltaf in E-P_phi phase space by Intergrating All Lambda!!
bool IAL = true;
/*****parameter for calculation of the evolution of distribution function******/
const int GRID_P_PHI      = 180;
double Lambda_0           = 0.80;
double E_0                = 0.25*v_A*v_A;
double windows_width      = 2e-2;
double P_phi_min;
double P_phi_max;
double dP_phi;
int    TIME_STEP_INTERVAL = 50;
double myf_vs_P_phi[GRID_P_PHI],mydeltaf_vs_P_phi[GRID_P_PHI],mynum[GRID_P_PHI];
double f_vs_P_phi[GRID_P_PHI],deltaf_vs_P_phi[GRID_P_PHI],num[GRID_P_PHI];
//second group variable of different pitch angle variable Lambda_0
double P_phi_min1;
double P_phi_max1;
double dP_phi1;
double myf_vs_P_phi1[GRID_P_PHI],mydeltaf_vs_P_phi1[GRID_P_PHI],mynum1[GRID_P_PHI];
double f_vs_P_phi1[GRID_P_PHI],deltaf_vs_P_phi1[GRID_P_PHI],num1[GRID_P_PHI];

//-----------------------------------------------------------------------------------
//extra parameter for f vs. P_phi and E 2d-plot
const int GRID_E      = 180;
double E_min;
double E_max;
double dE;
double myf_vs_P_phi_and_E[GRID_P_PHI][GRID_E],mydeltaf_vs_P_phi_and_E[GRID_P_PHI][GRID_E],mynum_2d[GRID_P_PHI][GRID_E];
double f_vs_P_phi_and_E[GRID_P_PHI][GRID_E],deltaf_vs_P_phi_and_E[GRID_P_PHI][GRID_E],num_2d[GRID_P_PHI][GRID_E];




/****************************************************************/


double v_par_max =  1.5*v_0;
double v_par_min = -1.5*v_0;
double v_per_max =  1.5*v_0;
double v_per_min =  0.0*v_0;
double v_max     =  1.5*v_0;
double v_cut     =  0.1*v_0;

//volume average beta vaule
const double beta = 0.020;

//total time
double total_time = 500/omega_Alfen;

//time step
double dt = 0.035*sqrt(alpha1);

//time initialization
double t=0;

//background damping rate
const double gamma_d[NUM_MODE]  = {0.000*omega_Alfen};
//const double gamma_d[NUM_MODE]  = {0.014*omega_Alfen,0.008*omega_Alfen};

//drag rate
const double nu                 = 0.0*omega_Alfen;

//pitch angle scattering effect can be added (true) or not (false)
bool SCATTERING                 = false;

//slowing down effect can be added (true) or not (false)
bool SLOWING_DOWN               = false;

//source and sink effect can be considered (true) or not (false)
bool SOURCE_AND_SINK            = false;

//particle loading in velocity space is uniform(true) or non-uniform(false) 
bool UNIFORM_LOADING            = true;

//parameter for TAE input
/* ndat   :number of radial grid
   mdat   :number of poloidal grid
   m_min  :minimum poloidal number
   m_max  :maximum poloidal number
   omega_A:mode frequency
   A      :mode amplitude
   n      :toridal mode number
*/
const int ndat = 201;
const int mdat = 21;
int ndat_NOVA,mode_number,m_min_NOVA;
const int m_min=-3,m_max=20;
double xi[m_max-m_min+1][ndat],data[26][ndat];

double X[NUM_MODE],Y[NUM_MODE];
//double X = B0*0.0003101,Y = B0*0.00000000;//beta = 0.0126
//double X = B0*0.0002879,Y = B0*0.00000000;//beta = 0.0026


//double X = B0*0.00003281,Y = B0*0.00000000;//beta = 0.00007
//double X = B0*0.00006685,Y = B0*0.00000000;//beta = 0.00012
//double X = B0*0.00011880,Y = B0*0.00000000;//beta = 0.00017
//double X = B0*0.00019230,Y = B0*0.00000000;//beta = 0.00022
//double X = B0*0.00028510,Y = B0*0.00000000;//beta = 0.00027
//double X = B0*0.00039830,Y = B0*0.00000000;//beta = 0.00032
//double X = B0*0.00064660,Y = B0*0.00000000;//beta = 0.00042
//double X = B0*0.00099160,Y = B0*0.00000000;//beta = 0.00052
//double X = B0*0.004825,Y = B0*0.00000000;//beta = 0.00027,time average case
//X[0] = B0*0.0005112,Y[0] = B0*0.00000000;
//X[1] = B0*0.0003800,Y[1] = B0*0.00000000;//omega_n_2 = 0.300,multiple mode case

//self-consistent wave evolution?
bool AMP = true;

//Output saving resonant particles orbit information
bool RES = false;

//saving Single Particle Information
bool SPI = false;

//How many steps record resonant particles
int RES_STEP = 1;

//full-f option
bool FULL_F = false;

//time average option
bool AVG = false;

//Method 1 : calc current on grid
//Method 2 : sum all particles
int METHOD = 2;

//weather using error function
bool ERF = true;

//deltaf-method:
//1. dw/dt = - (f/g - w)dlnf_0/dt
//2. dw/dt = - 1/g*df_0/dt
int DELTAF_METHOD = 2;

//how many wave period to average
const double WAVE_PERIOD_NUM=1;

//linear condition for particle orbit
bool LINEAR = false;

//weather include deltaB velocity ?
bool DELTAB = true;

//simple drift velocity.v_d = c*v_z
bool SIMPLE_V_D = false;

//assumption for weight equation
//condition 1 : no any assumption
//condition 2 : dp_phidt = n/omega_A*dEdt
//condition 3 : dEdt     = omega_A/n*dp_phidt
int WEIGHT_CONDITION = 1;

//weather phase is fixed?
bool PHASE_FIXED = false;

//weather A_par is derived from deltaPhi by integration or not?
bool A_PAR_INT = false;

//  eliminate adiabatic term using
//  0.retain adiabatic term : f_a = 0.0;
//  1.candy's method        : f_a = e*R*deltaA_||*B_phi/B df_0dP_phi + e*deltaPhi*df_0dE - mu*deltaB*B0/E*df_0dLambda
//  2.fu's method           : f_a = -xi \cdot \grad f = e*R*deltaA_||*B_phi/B df_0dP_phi
int ADIABATIC = 1;

//check dp_phidt
bool CHECK_DPDT = false;

//check dEdt
bool CHECK_DEDT = false;

//check dKdt
bool CHECK_DKDT = false;

//electrostatic case?
bool ELECTROSTATIC = false;

//simplst case:
//B*_// = B curl_B = 0,so curl_b = (b times grad_B)/B , E_// != 0
bool SIMPLE_CASE = false;

//including deltaPhi and deltaA_par
bool INCLUDE_PHI_A_PAR = true;

//unperturbed P_phi?
bool unperturbed_P_PHI = false;

//xi is uniform in r(true) or sqrt(psi)(false)
bool UNIFORM_RADIAL = false;

//diagnosis for single particle:
//v_per,v_par
double SINGLE_R     = R0 + 0.21*a;
double DELTA_R	    = 0.18*a;
double SINGLE_Z     = 0.00;
double SINGLE_phi   = 0.00;
double SINGLE_v_par = -0.5*v_A;
double SINGLE_v_per = 0.0*v_A;
double SINGLE_Lambda= 0.75;//add by 2013.10.26,use SINGLE_Lambda to define SINGLE_v_per;
double SINGLE_mu;
//add by 2014-09-18, fix mu and K0
bool   TEK_mu      = true;
double E_0_TEK     = 0.79*v_A*v_A;
double P_phi_0_TEK = -0.040;
double mu_0_TEK    = 0.75*v_A*v_A/B0;
double sigma_TEK   = -1;

//set same intial condition for all particles
bool SAME_CONDITION = false;

//dpdt and dEdt use hamiltonian descreption(true) or not(false)
bool HAMILTONIAN = false;

//magnetic field only has perpendicular part(true) or not (false)
bool B_PER = false;

//recycle method : 1: z = -z 2: remove particles
int RECYCLE_METHOD = 1;

//amplitude equation should be modify(true) or not(false)
bool AMP_MOD_TERM = false;

//derivative of phi is analytical soluntion (true) or not(false)
bool DERIVATIVE_OF_PHI = false;

//use Shafranov flux coordinates (true) or not(false)
bool SHAFRANOV = true;

//  eliminate adbatic term using
//  1.transform Poloidal angle	:	theta	--> theta_f
//  2.transform Toroidal angle	:	phi		-->	zeta
int FLUX_COORDINATE = 1;

//use two different condition for trapped particles and passing particles (true) or not(false)
bool AVERAGE_PSI = false;

//Jacobian has approxmation to O(e^2)(false) or not(true)
bool JACOBIAN_ACCURACY = true;

//Shafranov shift has q-profile(true) or constant(false)
bool SHIFT_ACCURACY = true;

//multiple mode (true) or not(false)
bool MULTIPLE_MODE = false;

//isotropy case(ture) or anisotropic case(false),here anisotropic is correspond to Lambda = 0!
bool ISOTROPY = true;

//numerical equilibrium
bool NUMERICAL_EQUILIBRIUM = false;

//Output K at t=0 for every Particle
bool OKP = true;


/***************Specific Pitch angle Loading****************/
//Specific Pitch angle Loading : Lambda = mu*B0/E
bool SPL = false;
double Lambda_0_load;
//Output Number of Specific pitch angle
bool ONS = true;
double Lambda_0_load_1 = 0.75;
double Lambda_0_load_2 = 0.99;
int N_Lambada_1 = N_p_initial;
int N_Lambada_2 = 0;
bool INPUT_C_F = false;
double c_f_specific = 3.55083e-12;
/***********************************************************/

/***************Finite width of Gaussian distrbution of Pitch angle****************/
//Specific Pitch angle Loading                         : Lambda = mu*B0/E
//Finite width of Gaussian distrbution of Pitch angle  : f(Lambda) ~ exp[(Lambda-Lambda_0)^2/Delta_Lambda^2]
bool FGP = false;
double Lambda_0_FGP = 0.50;
double Delta_Lambda = 0.10;

//----Double pitch angle distrbution with FGP          : f(Lambda) ~ c1_FGP*exp[(Lambda-Lambda_1)^2/Delta_Lambda_1^2]+c2_FGP*exp[(Lambda-Lambda_2)^2/Delta_Lambda_2^2]
bool DFGP = false;
double Lambda_1_FGP   = 0.25;
double Lambda_2_FGP   = 0.99;
double Delta_Lambda_1 = 0.05;
double Delta_Lambda_2 = 0.05;
double c1_FGP         = 0.3635;//coefficient of proportion for Lambda_1_FGP particle
double c2_FGP         = 0.6365;//coefficient of proportion for Lambda_2_FGP particle
/***********************************************************/

/***************Time Evolution of KAM surfaces****************/
bool TEK = false;


//output delta-f vs. phase space diagnostics (true) or not(false)
bool DIAGNOSTICS = true;

//how many time interval region of total simulation for diagnostics
int NUM_TIME_DIAGNOSTICS = 20;

/************output deltaf in E-P_phi phase space**************/
bool   dEP     = true;
//use mu Instead of Lambda to fix
bool   mIL  = true;
double mu_0_1 = 0.25*v_A*v_A/B0;
double mu_0_2 = 0.50*v_A*v_A / B0;
double mu_0_3 = 0.75*v_A*v_A / B0;
double mu_0_4 = 1.00*v_A*v_A / B0;


/**********output Deltaf in Lambda-P_phi phase space**********/
//only valid for isotropy case
bool   dLP               = false;
double E_0_dLP           = 0.5*v_A*v_A; 
double windows_width_dLP = 5e-2;

/************output Deltaf in Lambda-E phase space************/
//only valid for isotropy case
bool   dLE               = false;
double P_phi_0_dLE       = -0.040;
double windows_width_dLE = 5e-2;


/*********output Particle_Region in XX-Lambda space***********/
//XX is psi_average,R,P_phi,E
bool PRL = false;

//add by 2014.05.20
//EP current J_h Including Magnetization Current (true) or not(false)
bool IMC = true;


#endif
