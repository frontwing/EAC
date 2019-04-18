



#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <ctime>

#ifndef TEST_H_
#define TEST_H_
#include "tokamak.h"
#include "parameter.h"
using namespace std;

/*
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
*/



#include <omp.h>
#include <mpi.h>


clock_t start,finish;
double totaltime;

int mymphi=mphi;

double dphi;

//equilibrium field grid variable ~(R,Z)
scalar2d<double> psi(mR,mZ),B_value(mR,mZ),g_eq(mR,mZ);
vector2d<double> B(mR,mZ),b(mR,mZ),curl_b(mR,mZ),grad_B(mR,mZ),grad_g_eq(mR,mZ);

//pertubated field grid varable
scalar<double> deltaPhi_sin_3D[NUM_MODE];
scalar<double> deltaPhi_cos_3D[NUM_MODE];
vector<double>   deltaE_sin_3D[NUM_MODE];
vector<double>   deltaE_cos_3D[NUM_MODE];

//for out_put deltaB at axis
scalar<double> deltaA_par_cos_3D[NUM_MODE];
vector<double> deltaB_initial_3D[NUM_MODE];

scalar2d<double> Xi_initia_cos[NUM_MODE],Xi_initia_sin[NUM_MODE];
scalar2d<double> deltaA_par_sin[NUM_MODE],deltaA_par_cos[NUM_MODE],deltaPhi_sin[NUM_MODE],deltaPhi_cos[NUM_MODE];
vector2d<double> grad_deltaA_par_sin[NUM_MODE],grad_deltaA_par_cos[NUM_MODE],grad_deltaPhi_sin[NUM_MODE],grad_deltaPhi_cos[NUM_MODE];

//local pertubated field varable
double deltaE_local[3][2][2],deltaE_X_local[3][2][2],deltaE_Y_local[3][2][2],deltaB_local[3][2][2],deltaPhi_local[2][2],deltaA_par_local[2][2],grad_deltaA_par_local[3][2][2],grad_deltaPhi_local[3][2][2],grad_deltaPhi_X_local[3][2][2],grad_deltaPhi_Y_local[3][2][2];


double myJ_dot_E_X[NUM_MODE],myJ_dot_E_Y[NUM_MODE],J_dot_E_X[NUM_MODE],J_dot_E_Y[NUM_MODE];

//coordinate system convert form (r,theta,phi) -->  (R,Z,phi)
double r[mR][mZ],theta[mR][mZ],THETA_S[mR][mZ],THETA_F[mR][mZ];

//deltaN is change of per process
//particle number defition
int N_p = N_p_initial;
const int deltaN = (int)(0.1*N_p_initial) + N_p_diagnosis*8;

//transfer between two process variavle definition
Particle buffer_out[2*deltaN],buffer_in[2*deltaN];
Particle buffer_out_left[deltaN],buffer_in_left[deltaN],buffer_out_right[deltaN],buffer_in_right[deltaN];
Particle marker[N_p_initial + deltaN];
int N_lost,N_total;

//definition for particle species 
Species species[]={Species("ion")};

//input files list;
string input_list[NUM_MODE];

//per process boundary
double phi_left,phi_right;

//f ~ PDF(Physical  Distribution Function)
//g ~ NDF(Numerical Distribution Function) 
double f[N_p_initial+deltaN],g[N_p_initial+deltaN];

//Wave energy
double W[NUM_MODE];

//nonorthogonality of field variables
double mod_X_RHS[NUM_MODE],mod_Y_RHS[NUM_MODE];


//variables related to calculate particle Orbit Frequency
int    TIMESTEP;
double dZ_of_step1;
bool   Z_reverse_sign;

//variables related to time average
int T_RECORD=0;
int T_INITIA=0;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
const int T_AVG_NUM = (int)(WAVE_PERIOD_NUM*2*PI/abs(0.3132*omega_Alfen)/dt);
double *X_RECORD = new double[T_AVG_NUM];
double *Y_RECORD = new double[T_AVG_NUM];
double X_AVG = 0.0,Y_AVG = 0.0;

//index of resonant particle 
int *res_index = new int[N_p_diagnosis];

//normalization factor of distrubution funtion
double c_f;

//Time base ~ cos(n*phi+omega*t) & sin(n*phi+omega*t)
double COS_NPHI_PLUS_OMEGA_T,SIN_NPHI_PLUS_OMEGA_T;

//MPI variable
int myid,numprocs;

//input amplitude file X_input,Y_input,t_input
double *X_input,*Y_input,*t_input;

//calc phase changing
double alpha_current,alpha_previous,alpha_initial,Delta_alpha;


void equilibrium_field_2d();
void initia_TAE();
void initia_TAE_NOVA();
void r2R_accuracy();
void initia_output();
void input_TAE_NOVA(const int s);
void perturbed_field_cos(const int s);
void perturbed_field_sin(const int s);
void calc_wave_energy();
void loading_particle();
void amplitude();
void clean_create_files();
void initia_single_particle();
//add by 2013-10-27
void initia_single_particle_MPI();
void initia_parameter_output();
void stepon_RK4();
void kick_particle();
void diagnostics_movie_use_MPI_Gather(int index_of_time_record,double Lambda_0);
double get_B_value(const particle_vector<double> &X);
//add by 2012-3-23
bool calc_oribt_frequency(bool LOSS,bool OUTPUT_TRACK);
//add by 2011-6-20
void scattering();
//add by 2011-7-30
void source_and_sink();
void slowing_down();
//add by 2012-6-9
void numerical_equilibrium();
//add by 2012-7-22
void diagnostics_output_f_vs_P_phi(double Lambda_0,double E_0,double windows_width);
void diagnostics_output_f_vs_P_phi1(double Lambda_0,double E_0,double windows_width);
//add by 2012-9-19
void diagnostics_output_f_vs_P_phi_and_E(double Lambda_0,double windows_width);
//add by 2013-3-7
void diagnostics_on_number_of_particle();
//add by 2013-10-25
int load_amp_file();
//add by 2013-10-27
void output_Time_Evolution_of_KAM(int index_of_time_record);
//add by 2013-10-29
void output_perturbed_field_vs_time(double R_field,double Z_field,double phi_field);
//add by 2013-11-21
void diagnostics_output_deltaf_vs_P_phi_and_Lambda(int index_of_time_record,double E_0,double windows_width);
//add by 2013-11-21
void diagnostics_output_deltaf_vs_E_and_Lambda(int index_of_time_record,double P_phi_0,double windows_width);
//add by 2013-12-06
void output_particle_region();
//add by 2014-07-29
void diagnostics_movie_use_MPI_Gather_fixed_mu(int index_of_time_record, double mu_0);
//add by 2014-09-18
void initia_single_particle_MPI_fix_mu(double E_0, double P_phi_0, double mu_0, int sigma);


void equilibrium_field_2d()
{
	scalar2d<double> R(mR,mZ);
	for(int i=0;i<mR;i++)
		for(int j=0;j<mZ;j++)
			R[i][j] = R0-a*a_b+i*dR;


	//initia condition

	//q(psi) = q_0 + (psi - psi(0))/(psi(a)-psi(0)) , psi(a) = 0 , psi(0) = 1/(2*(q_0+0.5))*B_0*a*a; 
	//calculate poloidal flux psi
    if(NUMERICAL_EQUILIBRIUM)
        numerical_equilibrium();
    else
    {
	    if(SHAFRANOV)
	    {
		    for(int i=0;i<mR;i++)
			    for(int j=0;j<mZ;j++)
			    {
				    double r2 = r[i][j]*r[i][j];
				    psi[i][j] = psi_0*(q_0+1) + 0.5*sqrt(4*pow(psi_0*(q_0+1),2.0)-4*psi_0*B0*(r2-a*a));
			    }
	    }
	    else
	    {
            for(int i=0;i<mR;i++)
                for(int j=0;j<mZ;j++)
                {
                    double Z =   -a*a_b+j*dZ;
				    double R = R0-a*a_b+i*dR;
				    double r2 = (R-R0)*(R-R0)+Z*Z;
	    		    psi[i][j] = psi_0*(q_0+1) + 0.5*sqrt(4*pow(psi_0*(q_0+1),2.0)-4*psi_0*B0*(r2-a*a));
                }
        }
    }

    if(!NUMERICAL_EQUILIBRIUM)
        g_eq = B0*R0;

    grad_g_eq = grad(g_eq);
    


	//output psi(R,Z)
	if(myid == 0)
	{
		ofstream psi_RZ_output("psi_RZ.dat",ios::app);
        ofstream g_RZ_output("g_RZ.dat",ios::app);
		ofstream R_output("R.dat",ios::app);
		ofstream Z_output("Z.dat",ios::app);
		for(int i=0;i<mR;i++)
		{
			for(int j=0;j<mZ;j++)
			{
				double Z =   -a*a_b+j*dZ;
				double R = R0-a*a_b+i*dR;
                g_RZ_output   << g_eq[i][j] << " ";
				psi_RZ_output <<  psi[i][j] << " ";
				R_output << R << " ";
				Z_output << Z << " ";
			}
            g_RZ_output << endl;
			psi_RZ_output << endl;
			R_output << endl;
			Z_output << endl;
		}
        psi_RZ_output.close();
        g_RZ_output.close();
        R_output.close();
        Z_output.close();
	}

	//calculate three component of magnetic field
	  vector2d<double> BB(1.0/R*ddZ(psi),-1.0/R*ddR(psi),g_eq/R);
	  B = BB;


	//calculate scalar value of magnetic field
	  B_value = abs(B);


	//calculate unit vector of magnetic field 
      b=B/B_value;    
      
    //calculate curl b
      curl_b = curl(b);

    //calculate gradient B
      grad_B = grad(B_value);

}

void initia_TAE_NOVA()
{
	//input TAE mode structre

	if(myid ==0)
		cout << "Start reading TAE mode structure!!!!!!!!!!!!!" << endl;

	for(int s=0;s<NUM_MODE;s++)
    {
        char *input_file = (char*)input_list[s].c_str();
        ifstream fin(input_file);
	    if(!fin)
	    {
		    cout << "!!!!!!!!!!!outzpsi file does not exsit!!!!!!!!!!" << endl;
		    exit(1);
	    }


        string str;
	    int m;
	    getline(fin,str);// NTOR=    1, Om/OmA=  0.6546E+00
	    getline(fin,str);// NOSURF,  MTOT
	    fin >> ndat_NOVA >> mode_number;
	    double **xi_NOVA = new double* [mode_number]; 
	    for(int i=0;i<mode_number;i++) 
            xi_NOVA[i] = new double[ndat_NOVA];
	    getline(fin,str);//empty line
	    for(int i=0;i<mode_number;i++)
	    {
		    getline(fin,str);// ** m, n[0], Om/OmA =
		    fin >> m >> n[s] >> omega_A[s];
		    omega_A[s] = omega_A[s]/q_a/sqrt(alpha1);
    //add by 12.7
		    omega_A[s] = - omega_A[s];
		    if(i == 0)
			    m_min_NOVA = m;
		    for(int j=0;j<ndat_NOVA;j++)
			    fin  >> xi_NOVA[m - m_min_NOVA][j];
		    getline(fin,str);//empty line
        }
        //identify data is correct
	    if(myid ==0)
	    {
            string fh;
            if(NUM_MODE == 1)
                fh = "TAE_out.dat";
            else
            {
                stringstream ss;
                ss << s+1;
                fh = "TAE_out_mode" + ss.str() + ".dat";
            }
            char *filename = (char*)fh.c_str();
		    ofstream TAE_out(filename,ios::app);
		    for(int j=0;j<ndat_NOVA;j++)
		    {
			    for(int i=0;i<mode_number;i++)
				    TAE_out << setiosflags(ios::fixed) << setprecision(8) << xi_NOVA[i][j] << " ";
			    TAE_out << endl;
		    }
            fin.close();
            TAE_out.close();
        }
        if(myid ==0)
            cout << "Inputing is OK!!!!!!!!!!!!!" << endl;

        for(int j=0;j<ndat_NOVA;j++)
            for(int i=0;i<mode_number;i++)
                xi[i][j] = xi_NOVA[i][j];
        input_TAE_NOVA(s);
	    perturbed_field_cos(s);
	    perturbed_field_sin(s);
        for(int i=0;i<mode_number;i++)
			delete [] xi_NOVA[i];
		delete [] xi_NOVA;
    }
	

	

/**************************************************************************************************
			    output Elctric field and magnetic field on phi = 0 cross section
***************************************************************************************************/
	if(myid ==0)
	{
		ofstream E_field_out("E_field_initia.dat",ios::app);
		for(int j=0;j<mZ;j++)
		{
			for(int i=0;i<mR;i++)
				for(int k=1;k<2;k++)
					E_field_out << setiosflags(ios::scientific) << setprecision(12) << deltaE_cos_3D[0][0][i][j][k] << " " ;
			E_field_out << endl;
		}
        for(int s=0;s<NUM_MODE;s++)
        {
            string fh;
            if(NUM_MODE == 1)
                fh = "Xi_initia.dat";
            else
            {
                stringstream ss;
                ss << s+1;
                fh = "Xi_initia_mode" + ss.str() + ".dat";
            }
            char *filename = (char*)fh.c_str();
		    ofstream Xi_out(filename,ios::app);
		    for(int j=0;j<mZ;j++)
		    {
			    for(int i=0;i<mR;i++)
				    Xi_out << setiosflags(ios::fixed)  << setprecision(12) << Xi_initia_cos[s][i][j] << " " ;
			    Xi_out << endl;
		    }
            Xi_out.close();
        }

		ofstream Phi_field_out("phi_field_initia.dat",ios::app);
		for(int j=0;j<mZ;j++)
		{
			for(int i=0;i<mR;i++)
				Phi_field_out << setiosflags(ios::fixed)  << setprecision(12) << deltaPhi_cos[0][i][j] << " " ;
			Phi_field_out << endl;
		}
		ofstream A_par_field_out("A_par_field_initia.dat",ios::app);
		for(int j=0;j<mZ;j++)
		{
			for(int i=0;i<mR;i++)
			    if(A_PAR_INT)
				    A_par_field_out << setiosflags(ios::fixed)  << setprecision(12) << deltaA_par_sin[0][i][j] << " " ;
			    else
				    A_par_field_out << setiosflags(ios::fixed)  << setprecision(12) << deltaA_par_cos[0][i][j] << " " ;
			A_par_field_out << endl;
		}
        E_field_out.close();
        Phi_field_out.close();
        A_par_field_out.close();
    }
}


/********************************************************************************************
*                                                                                           *
*                           convert xi --> deltaPhi & deltaA_par                            *
*                    deltaPhi_mn = C_1*xi_mn deltaA_par_mn = C_2*deltaPhi_mn                *
*                                                                                           *
*                   where   C_1 = omega*B^2/(m*R*B_phi*J^{-1}+(n*B_0/q/R)^2)                *
*                           C_2 = 1/(omega*B)*(-m*J^{-1}+n*R*B_phi*|grad(phi)|^2)           *
*                                                                                           *
*     p.s. here we neglect m=0 component for deltaPhi due to singlerality of conversion     *
*                                                                                           *
********************************************************************************************/
void input_TAE_NOVA(const int s)
{
	double d_sqrt_psi = 1.0/(ndat_NOVA-1);
	double dr = a/(ndat_NOVA-1);
	double weight;

	deltaPhi_cos[s]   = 0.0;
    deltaPhi_sin[s]   = 0.0;
	deltaA_par_cos[s] = 0.0;
	deltaA_par_sin[s] = 0.0;
	Xi_initia_cos[s]  = 0.0;
	Xi_initia_sin[s]  = 0.0;

    deltaA_par_cos_3D[s] = 0.0;

    deltaPhi_cos_3D[s]   = 0.0;
	deltaPhi_sin_3D[s]   = 0.0;

	double aa=1.0;//if r^2 = (R-R0)^2 + Z^2 > 1 , aa = 0

	ofstream file_deltaPhi("deltaPhi.dat",ios::app);
	if(myid == 0)
		for(int i=0;i<mR;i++)
			for(int j=0;j<mZ;j++)
				if(r[i][j] < a)
					file_deltaPhi << setiosflags(ios::fixed) << setprecision(12) << r[i][j] << " " ;

	if(myid == 0)
		file_deltaPhi << endl;

	for(int m=m_min_NOVA;m<m_min_NOVA+mode_number;m++)
	{
		for(int i=0;i<mR;i++)
			for(int j=0;j<mZ;j++)
				for(int k=0;k<mymphi;k++)
				{
					double psi_par,q,r2;
					double Jacobian_1,coeff1,coeff2;
					double grad_parallel;
					double zeta_f,theta_f;
					double delta,ddeltadtheta,Delta_prime;
					double R = R0-a*a_b+i*dR;
					double Z =   -a*a_b+j*dZ;
					double phi = phi_left + k*dphi;
					if(SHAFRANOV)
					{
                        if(NUMERICAL_EQUILIBRIUM)
                        {
                            psi_par = psi[i][j];
                            if(psi[i][j] < 1e-7)
                            {
                                r2 = r[i][j]*r[i][j];
                                psi_par = psi_0*(q_0+1) + 0.5*sqrt(4*pow(psi_0*(q_0+1),2.0)-4*psi_0*B0*(r2-a*a));
                            }
                        }
                        else
                        {
                            r2 = r[i][j]*r[i][j];
                            psi_par = psi_0*(q_0+1) + 0.5*sqrt(4*pow(psi_0*(q_0+1),2.0)-4*psi_0*B0*(r2-a*a));
                        }
						q = q_0 + (psi_0 - psi_par)/(psi_0-psi_a);
						if(FLUX_COORDINATE == 1)
						{
							zeta_f	= phi;
							theta_f	= THETA_F[i][j];
                            //g_eq = B_phi*R
							if(JACOBIAN_ACCURACY)
								Jacobian_1 = g_eq[i][j]/q/(R*R);
							else
								Jacobian_1 = g_eq[i][j]/q/(R0*R0*(1+2*r[i][j]/R0*cos(THETA_F[i][j])));
						}
						else if(FLUX_COORDINATE == 2)
						{
							ddeltadtheta	= 5.0*r[i][j]/4/R0*cos(THETA_S[i][j]);
							delta			= 5.0*r[i][j]/4/R0*sin(THETA_S[i][j]);
							Delta_prime		= r[i][j]/4/R0;
                            zeta_f = phi - q*delta;
							theta_f = THETA_S[i][j];
							Jacobian_1 = B0/q/R/(1-Delta_prime*cos(THETA_S[i][j]));
						}
					}
					else
					{
						r2 = (R-R0)*(R-R0)+Z*Z;
						psi_par = psi_0*(q_0+1) + 0.5*sqrt(4*pow(psi_0*(q_0+1),2.0)-4*psi_0*B0*(r2-a*a));
						q = q_0 + (psi_0 - psi_par)/(psi_0);
//						double Jacobian_1 = B[2][i][j][k]*R0/q/R;
						Jacobian_1 = B0/q/R;
					}
					if(SHAFRANOV)
					{
						if(FLUX_COORDINATE == 1)
						{
							coeff1  =	omega_A[s]*pow(B_value[i][j],2.0)/(n[s]*pow(B0/q/R,2.0)*r2 + m*B[2][i][j]*R*Jacobian_1);
							coeff2 = 1.0/(omega_A[s]*B_value[i][j])*(n[s]*B[2][i][j]/R - m*Jacobian_1);
						}
						else if(FLUX_COORDINATE == 2)
						{
							coeff1  = omega_A[s]*pow(B_value[i][j],2.0)/(n[s]*pow(B0/q,2.0)*r2/(R*R) + (m - n[s]*q*ddeltadtheta)*B[2][i][j]*R*Jacobian_1);
							coeff2 = 1.0/(omega_A[s]*B_value[i][j])*(n[s]*B0*R0/(R*R) - (m - n[s]*q*ddeltadtheta)*Jacobian_1);
							//avoid singularity!!
							if((m - n[s]*q*ddeltadtheta) < 0)
							{
								coeff1 = 0.00;
								coeff2 = 0.00;
							}
						}
						grad_parallel = -(-m*Jacobian_1+n[s]*R0*B0/(R*R))/B_value[i][j];
                    }
					else
					{
						coeff1  =		   omega_A[s]*pow(B_value[i][j],2.0)/(n[s]*pow(B0/q/R,2.0)*r2 + m*B[2][i][j]*R*Jacobian_1);
						coeff2 = 1.0/(omega_A[s]*B_value[i][j])*(n[s]*B[2][i][j]/R - m*Jacobian_1);
					}
					int r_index;
                    bool OUT_OF_BOUDARY;
                    if(NUMERICAL_EQUILIBRIUM)
                    {
                        if(abs(psi[i][j]) > 1e-7)
                            OUT_OF_BOUDARY = false;
                        else
                            OUT_OF_BOUDARY = true;
                    }
                    else
                    {
                        if(r[i][j] < a)
                            OUT_OF_BOUDARY = false;
                        else
                            OUT_OF_BOUDARY = true;
                    }

                        
					if(!OUT_OF_BOUDARY)
					{
						if(UNIFORM_RADIAL)
						{
							r_index = (int)((r[i][j])/dr);
							weight  = ((r[i][j]) - r_index*dr)/dr;
							aa =1.0;
						}
						else
						{
							double norm_psi;
                            if(NUMERICAL_EQUILIBRIUM)
                                norm_psi = (psi[i][j] - psi_0)/(psi_a - psi_0);
                            else
                                norm_psi = (psi_par - psi_0)/(psi_a - psi_0);
							r_index = (int)(sqrt(norm_psi)/d_sqrt_psi);
							weight  = (sqrt(norm_psi) - r_index*d_sqrt_psi)/d_sqrt_psi;
							aa =1.0;
						}
					}
					else
					{
						r_index = 0;
						aa = 0;
					}
                    //neglect m=0 mode when calculate deltaPhi!!
					if(m == 0)
						aa = 0;
					if(SHAFRANOV)
					{
                        if(k == 0)
                        {
						    deltaPhi_cos[s][i][j]   += aa*coeff1*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1])*cos(m*theta_f);
						    deltaPhi_sin[s][i][j]   += aa*coeff1*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1])*sin(m*theta_f);
				
							
						    Xi_initia_cos[s][i][j]	+=        aa*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1])*cos(m*theta_f);
						    Xi_initia_sin[s][i][j]	+=        aa*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1])*sin(m*theta_f);
						    if(A_PAR_INT)
						    {
							    deltaA_par_cos[s][i][j] += aa*coeff1*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1])*cos(m*theta_f);
							    deltaA_par_sin[s][i][j] += aa*coeff1*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1])*sin(m*theta_f);
						    }
						    else
						    {
							    deltaA_par_cos[s][i][j] += aa*coeff1*coeff2*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1])*cos(m*theta_f);
							    deltaA_par_sin[s][i][j] += aa*coeff1*coeff2*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1])*sin(m*theta_f);
						    }
                        }

                        //output deltaB_max
                        if(A_PAR_INT)
                            deltaA_par_cos_3D[s][i][j][k] += aa*coeff1*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1])*cos(n[s]*zeta_f + m*theta_f);
                        else
                            deltaA_par_cos_3D[s][i][j][k] += aa*coeff1*coeff2*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1])*cos(n[s]*zeta_f + m*theta_f);

                        deltaPhi_cos_3D[s][i][j][k]  += aa*coeff1*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1])*cos(n[s]*zeta_f + m*theta_f);
						deltaPhi_sin_3D[s][i][j][k]  += aa*coeff1*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1])*sin(n[s]*zeta_f + m*theta_f);
					}
					else
					{
						deltaPhi_cos[s][i][j]   += aa*coeff1*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1])*cos(m*theta[i][j]);
						deltaPhi_sin[s][i][j]   += aa*coeff1*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1])*sin(m*theta[i][j]);
						Xi_initia_cos[s][i][j]	+=        aa*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1])*cos(m*theta[i][j]);
						Xi_initia_sin[s][i][j]	+=        aa*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1])*sin(m*theta[i][j]);
						if(A_PAR_INT)
						{
							deltaA_par_cos[s][i][j] += aa*coeff1*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1])*cos(m*theta[i][j]);
							deltaA_par_sin[s][i][j] += aa*coeff1*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1])*sin(m*theta[i][j]);
						}
						else
						{
							deltaA_par_cos[s][i][j] += aa*coeff1*coeff2*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1])*cos(m*theta[i][j]);
							deltaA_par_sin[s][i][j] += aa*coeff1*coeff2*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1])*sin(m*theta[i][j]);
						}
					}
					if(myid == 0 && r[i][j] < a && k == 0)
						file_deltaPhi << setiosflags(ios::fixed) << setprecision(12) << aa*coeff1*((1-weight)*xi[m-m_min_NOVA][r_index]+weight*xi[m-m_min_NOVA][r_index+1]) << " " ;
                }
				if(myid == 0)
					file_deltaPhi << endl;

                vector<double> b_3D(mR,mZ,mymphi);
                    for(int i=0;i<mR;i++)
                        for(int j=0;j<mZ;j++)
                            for(int k=0;k<mymphi;k++)
                                for(int d=0;d<3;d++)
                                    b_3D[d][i][j][k] = b[d][i][j];
                                    
                //output deltaB_max 
                deltaB_initial_3D[s] = curl(deltaA_par_cos_3D[s]*b_3D);
    }
    if(myid == 0)
        file_deltaPhi.close();
}

//compute cosine part of perturbed field
void perturbed_field_cos(const int s)
{
    vector<double> b_3D(mR,mZ,mymphi);
    for(int i=0;i<mR;i++)
        for(int j=0;j<mZ;j++)
            for(int k=0;k<mymphi;k++)
                for(int d=0;d<3;d++)
                    b_3D[d][i][j][k] = b[d][i][j];
    deltaE_cos_3D[s] = - grad(deltaPhi_cos_3D[s]) + (grad(deltaPhi_cos_3D[s])*b_3D)*b_3D;
	grad_deltaA_par_cos[s]  =   grad2d(deltaA_par_cos[s]);
	grad_deltaPhi_cos[s]    =   grad2d(deltaPhi_cos[s]);
}

//compute sine part of perturbed field
void perturbed_field_sin(const int s)
{
    vector<double> b_3D(mR,mZ,mymphi);
    for(int i=0;i<mR;i++)
        for(int j=0;j<mZ;j++)
            for(int k=0;k<mymphi;k++)
                for(int d=0;d<3;d++)
                    b_3D[d][i][j][k] = b[d][i][j];
    deltaE_sin_3D[s] = - grad(deltaPhi_sin_3D[s]) + (grad(deltaPhi_sin_3D[s])*b_3D)*b_3D;
	grad_deltaA_par_sin[s]  =   grad2d(deltaA_par_sin[s]);
	grad_deltaPhi_sin[s]	=   grad2d(deltaPhi_sin[s]);
}


//
/************************************************************************
*                                                                       *
*               convert (r,theta,phi) --> (R,Z,phi)                     *
*                                                                       *
*	S	:	Shafranov coordinates                                       *
*	f	:	flux coordinates                                            *
*	where,                                                              *
*	r_f		= r_S,                                                      *
*	theta_f	= theta_S - 2*eta*sin(theta_S),                             *
*	zeta_f	= zeta_S.                                                   *
*                                                                       *
*	R	= R0 - Dleta + r_S*cos(theta_S),                                *
*	Z	= r_S*sin(theta_S),                                             *
*	phi	= zeta_S.                                                       *
*                                                                       *
*	R	= R0 - Dleta + r_f*cos(theta_f) + eta*r_f*(cos(2*theta_f)-1),   *
*	Z	= r_f*sin(theta_f) + eta*r_f*sin(2*theta_f),                    *
*	phi	= zeta_f.                                                       *
*                                                                       *
*	x	=	r_f*cos(theta_S)                                            *
*	y	=	r_f*sin(theta_S)                                            *
*                                                                       *
*	where Delta(a) = 0                                                  *
*	                                                                    *
*	use Newton's method to find (r_S,theta_S) vs. (R,Z)                 *
*                                                                       *
************************************************************************/
void r2R_accuracy()
{
	if(SHAFRANOV)
	{

		double x,x_old,y;
		
		for(int i=0;i<mR;i++)
		{
			for(int j=0;j<mZ;j++)
			{
    			double Z =   -a*a_b+j*dZ;
				double R = R0-a*a_b+i*dR;
				double theta_initia;
				if(Z < 0)
					theta_initia = atan2(Z,R-R0)+ 2*PI;
				else
					theta_initia = atan2(Z,R-R0);
				//initial guess for x = r*cos(theta)
				x = sqrt(pow(R-R0,2.0) + pow(Z,2.0))*cos(theta_initia);
                if(sqrt(pow(R-R0,2.0)+pow(Z,2.0))>a)
                {
                    r[i][j] = sqrt(pow(R-R0,2.0) + pow(Z,2.0));
                    THETA_S[i][j] = 0.0;
                    THETA_F[i][j] = 0.0;
                }
                else
                {
				do
				    {
					    x_old = x;
					    double f,f_prime;
					    if(SHIFT_ACCURACY)
					    {
						    f		= -  (2*q_0+1)/(48*q_0*q_0*R0*a*a)*pow(x*x+Z*Z,2.0)   - 1.0/(8*R0)*(x*x+Z*Z) + x + R0 + (a*a)/8/R0 + (2*q_0+1)/(48*q_0*q_0*R0)*a*a - R;
						    f_prime = -4*(2*q_0+1)/(48*q_0*q_0*R0*a*a)*(x*x+Z*Z)*x        - 1.0/(4*R0)*x		 + 1;
					    }
					    else
					    {
						    f		= -1.0/(8*R0)*x*x + x + R0 - 1.0/8/R0*(Z*Z-a*a) - R;
						    f_prime = -1.0/(4*R0)*x   + 1;
					    }
					    x = x - f/f_prime;
				    }while(abs((x-x_old)/x_old) > 1e-5);
				    y = Z;
				    r[i][j]		= sqrt(x*x+y*y);
				    if(Z < 0)
					    THETA_S[i][j] = atan2(y,x)+ 2*PI;
				    else
					    THETA_S[i][j] = atan2(y,x);
				    //Dleta : Shafranov shfit
				    double Delta_prime,eta;
				    if(SHIFT_ACCURACY)
					    Delta_prime	= r[i][j]/(4*R0) + (2*q_0+1)/(12*q_0*q_0*R0*a*a)*pow(r[i][j],3.0);
				    else
					    Delta_prime	= r[i][j]/(4*R0);

				    eta      = 0.5*(Delta_prime+r[i][j]/R0);
				
				    THETA_F[i][j] = THETA_S[i][j] - 2*eta*sin(THETA_S[i][j]);
                }
			}
		}
	}
	else
	{
		for(int i=0;i<mR;i++)
			for(int j=0;j<mZ;j++)
			{
				double R = R0-a*a_b+i*dR;
				double Z =   -a*a_b+j*dZ;
				r[i][j] = sqrt(pow(R-R0,2.0) + pow(Z,2.0));
				if(Z < 0)
					theta[i][j] = atan2(Z,R-R0)+ 2*PI;
				else
					theta[i][j] = atan2(Z,R-R0);
			}
	}
}

void initia_output()
{
	if(myid==0)
	{
        ofstream q_output("q.dat",ios::app);
	    ofstream psi_output("psi.dat",ios::app);
		for(int i=0;i<mR;i++)
			for(int j=0;j<mZ;j++)
			{
				if(SHAFRANOV)
				{
					double r2 = r[i][j]*r[i][j];
					if(sqrt(r2) <= a)
					{
						double psi_par = psi_0*(q_0+1) + 0.5*sqrt(4*pow(psi_0*(q_0+1),2.0)-4*psi_0*B0*(r2-a*a));
						double q = q_0 + (psi_0 - psi_par)/(psi_0);
						q_output << setiosflags(ios::fixed) << setprecision(8) << q << " " << sqrt(r2) << " " << myid << endl;
						psi_output << setiosflags(ios::fixed) << setprecision(8) << psi_par << " " << sqrt(r2) << endl;
					}

				}
				else
				{
					double R = R0-a*a_b+i*dR;
					double Z =   -a*a_b+j*dZ;
					double r2 = (R-R0)*(R-R0)+Z*Z;
					if(sqrt(r2) <= a)
					{
						double psi_par = psi_0*(q_0+1) + 0.5*sqrt(4*pow(psi_0*(q_0+1),2.0)-4*psi_0*B0*(r2-a*a));
						double q = q_0 + (psi_0 - psi_par)/(psi_0);
						q_output << setiosflags(ios::fixed) << setprecision(8) << q << " " << sqrt(r2) << " " << myid << endl;
						psi_output << setiosflags(ios::fixed) << setprecision(8) << psi_par << " " << sqrt(r2) << endl;
					}
				}
			}
        q_output.close();
        psi_output.close();
    }

}


/************************************************************************
*                                                                       *
*                   Description for wave energy                         *
*   W_s :   <deltaE_X ^2/(mu_0*v_A^2)>_V                                *                                                                    *
*	W_c :   <deltaE_Y ^2/(mu_0*v_A^2)>_V                                *
*	                                                                    *
*	nonorthogonality :                                                  *
*	Delta_X = <2*deltaE_sin*times(deltaE_cos,B)/omega/(B*B)>_V,         *
*	Delta_Y	= Delta_X,                                                  *
*                                                                       *
************************************************************************/
void calc_wave_energy()
{
	
    for(int s=0;s<NUM_MODE;s++)
    {
        double myW = 0.0;

	    scalar<double> W_s(mR,mZ,mymphi),W_c(mR,mZ,mymphi);

        vector<double> B_3D(mR,mZ,mymphi);
        scalar<double> B_value_3D(mR,mZ,mymphi);
        for(int i=0;i<mR;i++)
            for(int j=0;j<mZ;j++)
                for(int k=0;k<mymphi;k++)
                {
                    B_value_3D[i][j][k] = B_value[i][j];
                    for(int d=0;d<3;d++)
                        B_3D[d][i][j][k] = B[d][i][j];
                }
        W_s = times(deltaE_sin_3D[s],B_3D)*times(deltaE_sin_3D[s],B_3D)/(B_value_3D*B_value_3D*B_value_3D*B_value_3D);
	    W_c = times(deltaE_cos_3D[s],B_3D)*times(deltaE_cos_3D[s],B_3D)/(B_value_3D*B_value_3D*B_value_3D*B_value_3D);

	    bool OUT_OF_BOUDARY;
        for(int i=0;i<mR;i++)
		    for(int j=0;j<mZ;j++)
			    for(int k=0;k<mymphi-1;k++)
                {
                    if(NUMERICAL_EQUILIBRIUM)
                    {
                        if(abs(psi[i][j]) > 1e-7)
                            OUT_OF_BOUDARY = false;
                        else
                            OUT_OF_BOUDARY = true;
                    }
                    else
                    {
                        if(r[i][j] < a)
                            OUT_OF_BOUDARY = false;
                        else
                            OUT_OF_BOUDARY = true;
                    }
				    if(!OUT_OF_BOUDARY)
				    {
					    double R  = R0-a*a_b+i*dR;
					    double dV = R*dR*dZ*dphi;
					    myW += 0.5*(W_s[i][j][k] + W_c[i][j][k])*alpha1*dV;
				    }
                }

    //add by 2011-4-26
	    scalar<double> mod_X(mR,mZ,mymphi),mod_Y(mR,mZ,mymphi);
	    double mymod_X,mymod_Y;
	    mod_X = 1.0/omega_A[s]*deltaE_sin_3D[s]*(times(deltaE_cos_3D[s],B_3D))/(B_3D*B_3D) - 1.0/omega_A[s]*deltaE_cos_3D[s]*(times(deltaE_sin_3D[s],B_3D))/(B_3D*B_3D);
	    mod_Y = 1.0/omega_A[s]*deltaE_cos_3D[s]*(times(deltaE_sin_3D[s],B_3D))/(B_3D*B_3D) - 1.0/omega_A[s]*deltaE_sin_3D[s]*(times(deltaE_cos_3D[s],B_3D))/(B_3D*B_3D);
	    mymod_X = 0.0;
	    mymod_Y = 0.0;
	    for(int i=0;i<mR;i++)
		    for(int j=0;j<mZ;j++)
			    for(int k=0;k<mymphi-1;k++)
                {
                    if(NUMERICAL_EQUILIBRIUM)
                    {
                        if(abs(psi[i][j]) > 1e-7)
                            OUT_OF_BOUDARY = false;
                        else
                            OUT_OF_BOUDARY = true;
                    }
                    else
                    {
                        if(r[i][j] < a)
                            OUT_OF_BOUDARY = false;
                        else
                            OUT_OF_BOUDARY = true;
                    }
				    if(!OUT_OF_BOUDARY)
				    {
					    double R  = R0-a*a_b+i*dR;
					    double dV = R*dR*dZ*dphi;
					    mymod_X += mod_X[i][j][k]*dV;
					    mymod_Y += mod_Y[i][j][k]*dV;
				    }
                }
	    MPI::COMM_WORLD.Reduce(&mymod_X,&mod_X_RHS[s],1,MPI::DOUBLE,MPI_SUM,0);
	    MPI::COMM_WORLD.Reduce(&mymod_Y,&mod_Y_RHS[s],1,MPI::DOUBLE,MPI_SUM,0);
	    
        
        mod_X_RHS[s] /= (2*PI*PI*a*a*R0);
	    mod_Y_RHS[s] /= (2*PI*PI*a*a*R0);

        mod_X_RHS[s] /= numprocs;
	    mod_Y_RHS[s] /= numprocs;

	    MPI::COMM_WORLD.Reduce(&myW,&W[s],1,MPI::DOUBLE,MPI_SUM,0);

	    //volume average -- <W> = (\int WdV)/V , V = (2*PI*R0)*(PI*a^2)
	    W[s] /= (2*PI*PI*a*a*R0);
        W[s] /= numprocs;



	    if(myid == 0)
		    cout << setiosflags(ios::fixed) << setprecision(15) << "Wave energy[" << s << "] : " << W[s] << " mod_X_RHS : " << mod_X_RHS[s] << " mod_Y_RHS :" << mod_Y_RHS[s] << endl;
    }
}


/************************************************************************************************
*                                                                                               *
*                               loading particles                                               *
*                                                                                               *                   
*                   f ~ PDF(Physical  Distribution Function)                                    *
*                   g ~ NDF(Numerical Distribution Function)                                    *
*                                                                                               *
* two type loading method :                                                                     *
*                           #1      uniform loading (e.g g is constant)                         *
*                           #2  non-uniform loading (e.g g=f)                                   *
* two type distribution   :                                                                     *
*                           #1  Slowing-down    distrubution function                           *
*                           #2  Maxwellian      distrubution function                           *
*                                                                                               *
*                                                                                               *
************************************************************************************************/
void loading_particle()
{
	double beta_f_avg,beta_g_avg;
	double beta_f_center,beta_g_center;
	double mybeta_f_center,mybeta_g_center;
	double dr;
	double weight_line_auxiliary[2][2],weight_square[2][2];
    double myN_not_use,N_not_use;

	//calc beta value
	scalar<double> beta_f(mR,mZ,mymphi),beta_g(mR,mZ,mymphi);
	double mybeta_f,mybeta_g;
	beta_f = 0.0;
	beta_g = 0.0;
	mybeta_f = 0.0;
	mybeta_g = 0.0;
    mybeta_f_center = 0.0;
    mybeta_g_center = 0.0;
	ofstream ff("f.dat",ios::app);
	ofstream f_vs_v("f_vs_v.dat",ios::app);
	ofstream f_vs_p_phi("f_vs_p_phi.dat",ios::app);

		
//  particle's position loading
/*------------------------------------------------------------------------------|
|                                                                               |
|                Random loading in real space (R^2,Z,phi)                       |
|                                                                               |
|------------------------------------------------------------------------------*/

	srand((unsigned)time(NULL)*(myid+1.0));
	for(int p=0;p<N_p;p++)
	{
		double R2,R,Z,phi;
        bool OUT_OF_BOUDARY=true;
		do
		{
			R2  =  pow(R0-a,2.0) + (pow(R0+a,2.0)-pow(R0-a,2.0))*(rand()*1.0/RAND_MAX);
			Z   = -a + 2*a*(rand()*1.0/RAND_MAX);
			phi =  2*PI*(rand()*1.0/(RAND_MAX+1.0));
			R   =  sqrt(R2);
			marker[p].X[0] = R;
			marker[p].X[1] = Z;
			marker[p].X[2] = phi;
			marker[p].id   = myid*N_p + p;
            if(NUMERICAL_EQUILIBRIUM)
            {
                double weight_line_auxiliary[2][2],weight_square[2][2];

		        int i = (int)((R-R0+a*a_b)/dR);
		        int j = (int)((Z+a*a_b)/dZ);
		
		        double R_grid   = R0-a*a_b+i*dR;
		        double Z_grid   =   -a*a_b+j*dZ;

		
		        weight_line_auxiliary[0][0] = (R - R_grid)/dR;
		        weight_line_auxiliary[1][0] = (Z - Z_grid)/dZ;

		        weight_line_auxiliary[0][1] = 1.0 - weight_line_auxiliary[0][0];
		        weight_line_auxiliary[1][1] = 1.0 - weight_line_auxiliary[1][0];
				

		        for(int ii=0;ii<2;ii++)
			        for(int jj=0;jj<2;jj++)
				        weight_square[ii][jj] = weight_line_auxiliary[0][1-ii]*weight_line_auxiliary[1][1-jj];

		        double psi_par = 0.0;
				
		        for(int ii=0;ii<2;ii++)
			        for(int jj=0;jj<2;jj++)
				        psi_par		+=	psi[i+ii][j+jj]*weight_square[ii][jj];
                
                if(abs(psi_par) > 1e-8)
                    OUT_OF_BOUDARY = false;
            }
            else
            {
                if(pow(R-R0,2.0)+pow(Z,2.0)<a*a)
                    OUT_OF_BOUDARY = false;
            }


             //guarantee (B0/B/Lambda_0) >0 ,reinject particle in real space if (B0/B/Lambda_0) <0
            if(SPL)
            {
                int i = (int)((marker[p].X[0]-R0+a*a_b)/dR);
		        int j = (int)((marker[p].X[1]+a*a_b)/dZ);
		
		        double R   = R0-a*a_b+i*dR;
		        double Z   =   -a*a_b+j*dZ;

		
		        weight_line_auxiliary[0][0] = (marker[p].X[0] - R)/dR;
		        weight_line_auxiliary[1][0] = (marker[p].X[1] - Z)/dZ;

		        weight_line_auxiliary[0][1] = 1.0 - weight_line_auxiliary[0][0];
		        weight_line_auxiliary[1][1] = 1.0 - weight_line_auxiliary[1][0];
				

		        for(int ii=0;ii<2;ii++)
			        for(int jj=0;jj<2;jj++)
				        weight_square[ii][jj] = weight_line_auxiliary[0][1-ii]*weight_line_auxiliary[1][1-jj];

		        double B_value_par = 0.0;
				
		        for(int ii=0;ii<2;ii++)
			        for(int jj=0;jj<2;jj++)
				        B_value_par		+=	B_value[i+ii][j+jj]*weight_square[ii][jj];

                if(p>=0 && p<N_Lambada_1)
                    Lambda_0_load = Lambda_0_load_1;
                else if(p>=N_Lambada_1 && p<N_Lambada_1+N_Lambada_2)
                    Lambda_0_load = Lambda_0_load_2;

                if((B0/B_value_par/Lambda_0_load - 1) < 0)
                    OUT_OF_BOUDARY = true;
            }


		}while(OUT_OF_BOUDARY);
    }

/*------------------------------------------------------------------------------|
|                                                                               |
|            Random loading particle in velocity space (v_par,v_per^2)          |
|                                                                               |
|------------------------------------------------------------------------------*/
    if(!SPL)//uniform pitch angle loading method
    {
        if(UNIFORM_LOADING)
        {
	        for(int p=0;p<N_p;p++)
	        {
		        //for slowing down case v<v_max
		        do
		        {
			        marker[p].v_per = sqrt(pow(v_per_min,2.0) + (pow(v_per_max,2.0)-pow(v_per_min,2.0))*(rand()*1.0/RAND_MAX));
			        marker[p].v_par = v_par_min + (v_par_max  - v_par_min)*(rand()*1.0/RAND_MAX);
		        }while(sqrt(pow(marker[p].v_par,2.0)+pow(marker[p].v_per,2.0))> v_max);
        //		uniform g
                g[p] = 1.0;
	        }
        }
        else
    /*----------------------------------------------------------------------------------------------------|
    |                                                                                                   |
    |  non-uniform loading particle in velocity space (v_par,v_per^2) using Rejection sampling method   |
    |                                                                                                   |
    |----------------------------------------------------------------------------------------------------*/
        {
            double f_v,g_v,v,M,lambda;
            int pm;
            M = 1.0/pow(v_c,3.0);
            for(int p=0;p<N_p;p++)
	        {
		        do
                {
                    //for slowing down case v<v_max
                
		            do
		            {
			            marker[p].v_per = sqrt(pow(v_per_min,2.0) + (pow(v_per_max,2.0)-pow(v_per_min,2.0))*(rand()*1.0/RAND_MAX));
			            marker[p].v_par = v_par_min + (v_par_max  - v_par_min)*(rand()*1.0/RAND_MAX);
		            }while(sqrt(pow(marker[p].v_par,2.0)+pow(marker[p].v_per,2.0))> v_max);
                
                    v   = sqrt(pow(marker[p].v_par,2.0)+pow(marker[p].v_per,2.0));
                    g_v = rand()*1.0/RAND_MAX;
                    f_v = 1.0/(pow(v_c,3.0) + pow(v,3.0));
                }while(g_v*M > f_v);
        //		non-uniform g
    		    g[p] = f_v;
            }
        }
    }
    else
    /*----------------------------------------------------------------------------------------------------|
    |                                                                                                   |
    |                Specific Pitch angle Loading method : Lambda = mu*B0/E                             |
    |   in order to guarantee (B0/B/Lambda-1)  > 0  ,Lambda_0_load should not be too large!!            |
    |    weight of particles that (B0/B/Lambda-1) < 0 is set to zero !!                                 | 
    |                                                                                                   |
    |----------------------------------------------------------------------------------------------------*/
    {
        for(int p=0;p<N_p;p++)
	    {
            //calc magnetic field
            int i = (int)((marker[p].X[0]-R0+a*a_b)/dR);
		    int j = (int)((marker[p].X[1]+a*a_b)/dZ);
		
		    double R   = R0-a*a_b+i*dR;
		    double Z   =   -a*a_b+j*dZ;

		
		    weight_line_auxiliary[0][0] = (marker[p].X[0] - R)/dR;
		    weight_line_auxiliary[1][0] = (marker[p].X[1] - Z)/dZ;

		    weight_line_auxiliary[0][1] = 1.0 - weight_line_auxiliary[0][0];
		    weight_line_auxiliary[1][1] = 1.0 - weight_line_auxiliary[1][0];
				

		    for(int ii=0;ii<2;ii++)
			    for(int jj=0;jj<2;jj++)
				    weight_square[ii][jj] = weight_line_auxiliary[0][1-ii]*weight_line_auxiliary[1][1-jj];

		    double B_value_par = 0.0;
				
		    for(int ii=0;ii<2;ii++)
			    for(int jj=0;jj<2;jj++)
				    B_value_par		+=	B_value[i+ii][j+jj]*weight_square[ii][jj];

		    //for slowing down case v<v_max
		    do
		    {
			    marker[p].v_par = v_par_min + (v_par_max  - v_par_min)*(rand()*1.0/RAND_MAX);
		    }while(sqrt(pow(marker[p].v_par,2.0))> v_max);
    //		uniform g
            g[p] = 1.0;
//  for multiple pitch angle case
//  N --- [0,1 ... ,N_Lambada_1-1,N_Lambada_1,...N_Lambada_1+N_Lambada_2-1]
            if(p>=0 && p<N_Lambada_1)
                Lambda_0_load = Lambda_0_load_1;
            else if(p>=N_Lambada_1 && p<N_Lambada_1+N_Lambada_2)
                Lambda_0_load = Lambda_0_load_2;
            marker[p].v_per = abs(marker[p].v_par)/sqrt(B0/B_value_par/Lambda_0_load - 1);

            
	    }
    }






/*------------------------------------------------------------------------------|
|                                                                               |
|                       Compute Distrubution Function                           |
|                                                                               |
|------------------------------------------------------------------------------*/

	for(int p=0;p<N_p;p++)
	{
		int i = (int)((marker[p].X[0]-R0+a*a_b)/dR);
		int j = (int)((marker[p].X[1]+a*a_b)/dZ);
		
		double R   = R0-a*a_b+i*dR;
		double Z   =   -a*a_b+j*dZ;

		
		weight_line_auxiliary[0][0] = (marker[p].X[0] - R)/dR;
		weight_line_auxiliary[1][0] = (marker[p].X[1] - Z)/dZ;

		weight_line_auxiliary[0][1] = 1.0 - weight_line_auxiliary[0][0];
		weight_line_auxiliary[1][1] = 1.0 - weight_line_auxiliary[1][0];
				

		for(int ii=0;ii<2;ii++)
			for(int jj=0;jj<2;jj++)
				weight_square[ii][jj] = weight_line_auxiliary[0][1-ii]*weight_line_auxiliary[1][1-jj];

		particle_vector<double> B_par;
		double B_value_par;
		B_par       = 0.0;
		B_value_par = 0.0;
				
		for(int ii=0;ii<2;ii++)
			for(int jj=0;jj<2;jj++)
			{
				B_value_par		+=	B_value[i+ii][j+jj]*weight_square[ii][jj];
				for(int d=0;d<2;d++)
					B_par[d]      +=      B[d][i+ii][j+jj]*weight_square[ii][jj];
			}

        if(NUMERICAL_EQUILIBRIUM)
        {
            double g_eq_par = 0.0;
            for(int ii=0;ii<2;ii++)
			    for(int jj=0;jj<2;jj++)
				    g_eq_par		+=	g_eq[i+ii][j+jj]*weight_square[ii][jj];

            B_par[2] = g_eq_par/marker[p].X[0];
        }
        else
            B_par[2] = B0*R0/marker[p].X[0];

		
		double psi_par =0.0;
		if(SHAFRANOV)
		{
			psi_par = 0;
			for(int ii=0;ii<2;ii++)
				for(int jj=0;jj<2;jj++)
					psi_par += psi[i+ii][j+jj]*weight_square[ii][jj];
		}
		else
		{
			double r2 = pow(marker[p].X[0]-R0,2.0)+pow(marker[p].X[1],2.0);
			psi_par = psi_0*(q_0+1) + 0.5*sqrt(4*pow(psi_0*(q_0+1),2.0)-4*psi_0*B0*(r2-a*a));
			R = marker[p].X[0];
			Z = marker[p].X[1];
			double q = q_0 + (psi_0 - psi_par)/(psi_0-psi_a);
			B_par[0]    =   B0/q/R*Z;
			B_par[1]    = - B0/q/R*(R-R0);
			B_par[2]    =   B0*R0/R;
			B_value_par =   sqrt(B_par[0]*B_par[0]+B_par[1]*B_par[1]+B_par[2]*B_par[2]);
		}

		/***************** use mu *****2011.10.9*************/
		marker[p].mu = 0.5*species[0].mass*pow(marker[p].v_per,2.0)/B_value_par;

        marker[p].P_phi =  (species[0].mass*marker[p].v_par*marker[p].X[0]*(B_par[2]/B_value_par)+species[0].charge*psi_par/alpha);



/*----------------------------------------------------------------------------------------------|
|                                                                                               |
|       f is slowing down distrubution : f ~ 1/(v^3+v_c^3)*exp(-(p_phi/(e*c_1*Delta psi)))      |
|                                                                                               |
|----------------------------------------------------------------------------------------------*/		

        double p_phi = (species[0].mass*marker[p].v_par*marker[p].X[0]*(B_par[2]/B_value_par)+species[0].charge*psi_par/alpha);
	
		double E  = 0.5*species[0].mass*pow(marker[p].v_par,2.0) + marker[p].mu*B_value_par;
		double v = sqrt(2*E/species[0].mass);
		if(AVERAGE_PSI)
		{
			double bracket_psi;
			double mu = marker[p].mu;
			double Lambda = marker[p].mu*B0/E;

			int sgn_v_parallel = marker[p].v_par>0?1:-1;
			if((1-mu*B0/E) > 0)
				bracket_psi = p_phi/species[0].charge - species[0].mass/species[0].charge*sgn_v_parallel*v*R0*sqrt(1-mu*B0/E);
			else
				bracket_psi = p_phi/species[0].charge;
			//add by 2014.10.22
			//Finite width of Gaussian distribution of Pitch angle  : f(Lambda) ~ exp[(Lambda-Lambda_0)^2/Delta_Lambda^2]
			if(FGP)
            {
                if(DFGP)//double pitch angle distribution
                    f[p] = 1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(bracket_psi/(c_1*Delta_psi)))*(c1_FGP*exp(-pow((Lambda-Lambda_1_FGP)/Delta_Lambda_1,2.0)) + c2_FGP*exp(-pow((Lambda-Lambda_2_FGP)/Delta_Lambda_2,2.0)));
                else
                    f[p] = 1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(bracket_psi/(c_1*Delta_psi)))*exp(-pow((Lambda-Lambda_0_FGP)/Delta_Lambda,2.0));

            }
			else
				f[p] = 1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(bracket_psi/(c_1*Delta_psi)));

		}
		else
		{
			double bracket_psi = p_phi/species[0].charge;
			//add by 2013.9.20
			//Finite width of Gaussian distribution of Pitch angle  : f(Lambda) ~ exp[(Lambda-Lambda_0)^2/Delta_Lambda^2]
			double Lambda = marker[p].mu*B0/E;
			if(FGP)
            {
                if(DFGP)//double pitch angle distribution
                    f[p] = 1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(bracket_psi/(c_1*Delta_psi)))*(c1_FGP*exp(-pow((Lambda-Lambda_1_FGP)/Delta_Lambda_1,2.0)) + c2_FGP*exp(-pow((Lambda-Lambda_2_FGP)/Delta_Lambda_2,2.0)));
                else
                    f[p] = 1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(bracket_psi/(c_1*Delta_psi)))*exp(-pow((Lambda-Lambda_0_FGP)/Delta_Lambda,2.0));

            }
			else
				f[p] = 1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(bracket_psi/(c_1*Delta_psi)));
		}
		

		if(myid == 0)
			ff << setiosflags(ios::fixed) << setprecision(8) << f[p] << " " << sqrt(pow(marker[p].X[0]-R0,2.0)+pow(marker[p].X[1],2.0)) << " " << v << " " << p_phi << endl;

		if(myid == 0)
		{
			if(ERF)
				f_vs_v     << 1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A))) << " " << v     << endl;
			else
				f_vs_v     << 1.0/(pow(v,3)+v_c*v_c*v_c) << " " << v     << endl;
			f_vs_p_phi << exp(-(p_phi/(species[0].charge*c_1*(psi_a - psi_0)))) << " " << p_phi << endl;
            f_vs_v.close();
            f_vs_p_phi.close();
		}


/*--------------------------------------------------------------------------------------|
|                                                                                       |
|           calculate volume average fast ion beta in code                              |
|                                                                                       |
|--------------------------------------------------------------------------------------*/

         mybeta_f += (species[0].mass*pow(marker[p].v_par,2.0) + marker[p].mu*B_value_par)/B_value_par/B_value_par*f[p]/g[p];


/*------------------------------------------------------------------------------------------|
|                                                                                           |
|                              compute center fast ion beta                                 |
|                                                                                           |
|------------------------------------------------------------------------------------------*/

		

		dr = 4*dR;
		double r_par;
		if(SHAFRANOV)
		{
/**********************************************************************

		S	:	Shafranov coordinates
		f	:	flux coordinates
		where,
		r_f		= r_S,
		theta_f	= theta_S - 2*eta*sin(theta_S),
		zeta_f	= zeta_S.

		R	= R0 - Dleta + r_S*cos(theta_S),
		Z	= r_S*sin(theta_S),
		phi	= zeta_S.

		R	= R0 - Dleta + r_f*cos(theta_f) + eta*r_f*(cos(2*theta_f)-1),
		Z	= r_f*sin(theta_f) + eta*r_f*sin(2*theta_f),
		phi	= zeta_f.

		x	=	r_f*cos(theta_f)
		y	=	r_f*sin(theta_f)
	
		use Newton's method to find (r_f,theta_f) vs. (R,Z) 

*************************************************************************/
			double x,x_old,y;
		
			R = marker[p].X[0];
			Z = marker[p].X[1];
			double theta_initia;
			if(Z < 0)
				theta_initia = atan2(Z,R-R0)+ 2*PI;
			else
				theta_initia = atan2(Z,R-R0);
			x = sqrt(pow(R-R0,2.0) + pow(Z,2.0))*cos(theta_initia);
			do
			{
				x_old = x;
				double f,f_prime;
				if(SHIFT_ACCURACY)
				{
					f		= -(2*q_0+1)/(48*q_0*q_0*R0*a*a)*pow(x*x+Z*Z,2.0)     - 1.0/(8*R0)*x*x + x + R0 + (a*a)/8/R0 + (2*q_0+1)/(48*q_0*q_0*R0)*a*a - R;
					f_prime = -4*(2*q_0+1)/(48*q_0*q_0*R0*a*a)*(x*x+Z*Z)*x         -1.0/(4*R0)*x   + 1;
				}
				else
				{
					f		= -1.0/(8*R0)*x*x + x + R0 - 1.0/8/R0*(Z*Z-a*a) - R;
					f_prime = -1.0/(4*R0)*x   + 1;
				}
				x = x - f/f_prime;
			}while(abs((x-x_old)/x_old) > 1e-5);
			y = Z;
			r_par		= sqrt(x*x+y*y);
		}
		else
			r_par  = sqrt(pow(marker[p].X[0]-R0,2.0) + pow(marker[p].X[1],2.0));
/*
		if(SHAFRANOV)
		{
			R = marker[p].X[0];
			Z = marker[p].X[1];
			//use Newton's method to find (r_f,theta_f) vs. (R,Z) 
			double x,x_old,y;
			x = R - R0;
			do
			{
				x_old = x;
				double f,f_prime;
				f		= R0 - 1.0/(8*R0)*x*x+x-11.0/(8*R0)*Z*Z/pow(1+5.0/4/R0*x,2.0)-R;
				f_prime =     -1.0/(4*R0)*x  +1+55.0/16/(R0*R0)*Z*Z/pow(1+5.0/4/R0*x,3.0);
				x = x - f/f_prime;
			}while(abs(x-x_old)/x_old > 1e-5);
			y = Z/(1+5.0/4/R0*x);
			r_par		= sqrt(x*x+y*y);
		}
		else
			r_par  = sqrt(pow(marker[p].X[0]-R0,2.0) + pow(marker[p].X[1],2.0));
*/
		if(r_par < dr)
            mybeta_f_center += (species[0].mass*pow(marker[p].v_par,2.0) + marker[p].mu*B_value_par)/B_value_par/B_value_par*f[p]/g[p];

	}

/*------------------------------------------------------------------------------------------|
|                                                                                           |
|                       output volume average fast ion beta radial profile                  |
|                                                                                           |
|------------------------------------------------------------------------------------------*/
	
	double mybeta_out[mR];
	double beta_out[mR];

	for(int i=0;i<mR;i++)
		for(int j=mZ/2;j<=mZ/2;j++)
		{
			mybeta_out[i] = 0.0;
			for(int k=0;k<mymphi;k++)
				if(r[i][j] < a)
					mybeta_out[i] += beta_f[i][j][k];
			mybeta_out[i] /= (mymphi-1);
		}
	MPI::COMM_WORLD.Reduce(&mybeta_out,&beta_out,mR,MPI::DOUBLE,MPI_SUM,0);

	if(myid == 0)
	{
		ofstream file_beta("beta.dat",ios::app);
		for(int i=0;i<mR/2;i++)
			for(int j=mZ/2;j<=mZ/2;j++)
			{
				if(r[i][j] < a)
				{
					beta_out[i] /= numprocs;
					file_beta << beta_out[i] << "  " << r[i][j] << endl;
				}
			}
        file_beta.close();
	}


	MPI::COMM_WORLD.Reduce(&mybeta_f,&beta_f_avg,1,MPI::DOUBLE,MPI_SUM,0);

	MPI::COMM_WORLD.Bcast(&beta_f_avg,1,MPI::DOUBLE,0);


//  beta_code_volume_average = (\sum beta_code*dV)/V

    
    if(NUMERICAL_EQUILIBRIUM)
    {
        double V=0;
        for(int i=0;i<mR;i++)
		    for(int j=0;j<mZ;j++)
				    if(abs(psi[i][j]) > 1e-7)
				    {
					    double R  = R0-a*a_b+i*dR;
					    double dV = R*dR*dZ*2*PI;
                        V += dV;
				    }
        beta_f_avg = beta_f_avg/V;
    }
    else
        beta_f_avg = beta_f_avg/(2*PI*PI*a*a*R0);


//normalize by beta
	c_f = beta/beta_f_avg;


/*------------------------------------------------------------------------------------------|
|                                                                                           |
|                           compute center fast ion beta                                    |
|                                                                                           |
|------------------------------------------------------------------------------------------*/


	MPI::COMM_WORLD.Reduce(&mybeta_f_center,&beta_f_center,1,MPI::DOUBLE,MPI_SUM,0);
	MPI::COMM_WORLD.Reduce(&mybeta_g_center,&beta_g_center,1,MPI::DOUBLE,MPI_SUM,0);

	MPI::COMM_WORLD.Bcast(&beta_f_center,1,MPI::DOUBLE,0);
	MPI::COMM_WORLD.Bcast(&beta_g_center,1,MPI::DOUBLE,0);



//  beta_code_center
	beta_f_center = beta_f_center/(2*PI*PI*dr*dr*R0);
/*****************************************************************************************************
//normalize by beta
	double c_f = beta/beta_f_center;
********************************************************************************************************/
	


	cout << " c_f : " << c_f << " beta_f_center : " << beta_f_center << " beta_f_avg : " << beta_f_avg << " beta_f_center/beta_f_avg : " << beta_f_center/beta_f_avg << " myid : " << myid << endl;
//	exit(0);


	for(int p=0;p<N_p;p++)
	{
		f[p] = f[p]*c_f;
		marker[p].f_over_g = f[p]/g[p];
		marker[p].g = g[p];
	}

/*------------------------------------------------------------------------------------------|
|                                                                                           |
|                   output initial position and distrubution function                       |
|                                                                                           |
|------------------------------------------------------------------------------------------*/
	if(myid==0)
	{
		ofstream f_initial("initial.dat",ios::app);
		for(int p=0;p<N_p;p++)
			f_initial << setiosflags(ios::fixed) << setprecision(8) << marker[p].X[0] << " " << marker[p].X[1] << " " << marker[p].X[2] <<  " " << marker[p].v_par << " " << marker[p].v_per << endl;
		f_initial.close();
		ff.close();
		
		cout << "Particle loading is Ok!!!!!!" << endl;
	}


	if(myid==0)
	{
		cout << "volume average fast ion beta : " << setiosflags(ios::fixed) << setprecision(16) << beta_f_avg << endl;
		cout << "max lamor radius (r/a)       : " << setiosflags(ios::fixed) << setprecision(6)  << species[0].mass*v_0/species[0].charge/B0/a << endl;
	}



/*****************************************************************************************************
                 if choose specific pitch angle case,c_f should be redifined!
*****************************************************************************************************/
    if(SPL && INPUT_C_F)
        c_f = c_f_specific;

}

//clean files
void clean_create_files()
{

	if(myid==0)
	{
		system("rm *.dat");
        system("rm parameter_output");
	}




}

/************************************************************************************************
*                                                                                               *
*                       Initial condition for test particle                                     *
*                                                                                               *
************************************************************************************************/
void initia_single_particle()
{
	if(myid == 0)
    {
		double K0;
		double B_value_par,B_par[3],deltaPhi_par,deltaA_par_par;
		double weight_line_auxiliary[2][2];
		double weight_square[2][2];
		double phi;
		double e,v_par,m;
		int i,j,k;
		double psi_par,q,r2;
		double R,Z;
        double lambda;
        double W_per,W_par;
		//initialization fo #0 particle
		marker[0].X[0]	= SINGLE_R;
		marker[0].X[1]	= SINGLE_Z;
		marker[0].X[2]	= SINGLE_phi;
		marker[0].v_par	= SINGLE_v_par;
		marker[0].v_per	= SINGLE_v_per;
		marker[0].id	= 0;

		q = 1;
		lambda = marker[0].v_par/sqrt(pow(marker[0].v_par,2.0)+pow(marker[0].v_per,2.0));
		cout << "v_wave(m=1) : " << (omega_A[0]/((n[0]-0/q)/marker[0].X[0]))/v_A << "v_A" << endl;
		cout << "v_wave(m=1) : " << (omega_A[0]/((n[0]-2/q)/marker[0].X[0]))/v_A << "v_A" << endl;
		cout << "v_wave(m=2) : " << (omega_A[0]/((n[0]-1/q)/marker[0].X[0]))/v_A << "v_A" << endl;
		cout << "v_wave(m=2) : " << (omega_A[0]/((n[0]-3/q)/marker[0].X[0]))/v_A << "v_A" << endl;
		cout << "v_par[0]    : " << marker[0].v_par/v_A << "v_A" << endl;
		cout << "v_per[0]    : " << marker[0].v_per/v_A << "v_A" << endl;

		i = (int)((marker[0].X[0]-R0+a*a_b)/dR);
		j = (int)((marker[0].X[1]+a*a_b)/dZ);

		R   = R0-a*a_b+i*dR;
		Z   =   -a*a_b+j*dZ;

		weight_line_auxiliary[0][0] = (marker[0].X[0] - R)/dR;
		weight_line_auxiliary[1][0] = (marker[0].X[1] - Z)/dZ;

		weight_line_auxiliary[0][1] = 1.0 - weight_line_auxiliary[0][0];
		weight_line_auxiliary[1][1] = 1.0 - weight_line_auxiliary[1][0];
				

		for(int ii=0;ii<2;ii++)
			for(int jj=0;jj<2;jj++)
				weight_square[ii][jj] = weight_line_auxiliary[0][1-ii]*weight_line_auxiliary[1][1-jj];

		psi_par		= 0.0;
		B_value_par = 0.0;
				
		for(int ii=0;ii<2;ii++)
			for(int jj=0;jj<2;jj++)
			{
				B_value_par	+=	B_value[i+ii][j+jj]*weight_square[ii][jj];
				psi_par		+=	    psi[i+ii][j+jj]*weight_square[ii][jj];
			}

        if(NUMERICAL_EQUILIBRIUM)
        {
            double g_eq_par = 0.0;
            for(int ii=0;ii<2;ii++)
			    for(int jj=0;jj<2;jj++)
				    g_eq_par		+=	g_eq[i+ii][j+jj]*weight_square[ii][jj];

            B_par[2] = g_eq_par/marker[0].X[0];
        }
        else
            B_par[2] = B0*R0/marker[0].X[0];

/*----------------------------------------------------------------------------------------------------|
|                                                                                                   |
|                Specific Pitch angle Choosing : Lambda = mu*B0/E                                   |
|   in order to guarantee (B0/B/Lambda-1)  > 0  ,Lambda_0_load should not be too large!!            |
|    weight of particles that (B0/B/Lambda-1) < 0 is set to zero !!                                 |
|                                                                                                   |
|----------------------------------------------------------------------------------------------------*/
        if(SPL)
        {
            //if B0/B/Lambda_0 - 1 < 0,marker[p].v_per is set to be -1 and weight of that partciels is zero!!
            if((B0/B_value_par/Lambda_0_load_1 - 1) < 0)
            {
                cout << " (B0/B/Lambda-1) < 0 !!!!";
                exit(0);
            }
            else
                marker[0].v_per = abs(marker[0].v_par)/sqrt(B0/B_value_par/Lambda_0_load_1 - 1);
        }


        //if B0/B/Lambda_0 - 1 < 0,marker[p].v_per is set to be -1 and weight of that partciels is zero!!
        //if Lambda !=0
        if(SINGLE_Lambda > 1e-8)
        {
            if((B0/B_value_par/SINGLE_Lambda - 1) < 0)
            {
                cout << " (B0/B/Lambda-1) < 0 !!!!";
                exit(0);
            }
            else
                marker[0].v_per = abs(marker[0].v_par)/sqrt(B0/B_value_par/SINGLE_Lambda - 1);
        }
        else
            marker[0].v_per = 0.0;





		{
			COS_NPHI_PLUS_OMEGA_T = cos(n[0]*marker[0].X[2]-omega_A[0]*t);
			SIN_NPHI_PLUS_OMEGA_T = sin(n[0]*marker[0].X[2]-omega_A[0]*t);


			for(int I=i;I<=i+1;I++)
				for(int J=j;J<=j+1;J++)
				{
					deltaA_par_local[I-i][J-j]					=	X[0]*deltaA_par_cos[0][I][J]*COS_NPHI_PLUS_OMEGA_T - X[0]*deltaA_par_sin[0][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[0]*deltaA_par_sin[0][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[0]*deltaA_par_cos[0][I][J]*SIN_NPHI_PLUS_OMEGA_T;
					deltaPhi_local[I-i][J-j]					=	X[0]*deltaPhi_cos[0][I][J]*COS_NPHI_PLUS_OMEGA_T - X[0]*deltaPhi_sin[0][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[0]*deltaPhi_sin[0][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[0]*deltaPhi_cos[0][I][J]*SIN_NPHI_PLUS_OMEGA_T;
				}
		

			deltaA_par_par = 0.0;
			deltaPhi_par = 0.0;

				
			for(int ii=0;ii<2;ii++)
				for(int jj=0;jj<2;jj++)
				{
					deltaA_par_par	+= deltaA_par_local[ii][jj]*weight_square[ii][jj];
					deltaPhi_par	+= deltaPhi_local[ii][jj]*weight_square[ii][jj];
				}
		}




		//only for passing particles,so W_per = 0.0
		W_par = 0.5*species[0].mass*pow(marker[0].v_par,2.0);
		W_per = 0.5*species[0].mass*pow(marker[0].v_per,2.0);

		marker[0].mu	=	0.5*species[0].mass*pow(marker[0].v_per,2.0)/B_value_par;

		R	  = marker[0].X[0];
		e     = species[0].charge;
		v_par = marker[0].v_par;
		m     = species[0].mass;
		K0 =  e*(R*deltaA_par_par)/alpha*B_par[2]/B_value_par + e*psi_par/alpha + m*v_par*R*B_par[2]/B_value_par - n[0]/omega_A[0]*(W_per+W_par+e*deltaPhi_par);

        cout << "marker[0].v_par  : " << marker[0].v_par/v_A  << "v_A  marker[0].X[0]  : " << marker[0].X[0] << " K0 : " << K0 << endl ;
			
		res_index[0] = marker[0].id;

        if(!COF)
        {
	        srand((unsigned)time(NULL)*(myid+1.0));
		    for(int p=1;p<N_p_diagnosis;p++)
            {
			    //define position
			
			    marker[p].X[0] = SINGLE_R - DELTA_R + 2*DELTA_R*p/(N_p_diagnosis-2);
			    marker[p].X[1] = 0.00;
                marker[p].X[2] = 2*PI*(rand()*1.0/RAND_MAX);
			    marker[p].id   = p;

			    i = (int)((marker[p].X[0]-R0+a*a_b)/dR);
			    j = (int)((marker[p].X[1]+a*a_b)/dZ);


			    R   = R0-a*a_b+i*dR;
			    Z   =   -a*a_b+j*dZ;

			    weight_line_auxiliary[0][0] = (marker[p].X[0] - R)/dR;
			    weight_line_auxiliary[1][0] = (marker[p].X[1] - Z)/dZ;

			    weight_line_auxiliary[0][1] = 1.0 - weight_line_auxiliary[0][0];
			    weight_line_auxiliary[1][1] = 1.0 - weight_line_auxiliary[1][0];
				

			    for(int ii=0;ii<2;ii++)
				    for(int jj=0;jj<2;jj++)
					    weight_square[ii][jj] = weight_line_auxiliary[0][1-ii]*weight_line_auxiliary[1][1-jj];

			    psi_par		= 0.0;
			    B_value_par = 0.0;
				
			    for(int ii=0;ii<2;ii++)
				    for(int jj=0;jj<2;jj++)
				    {
					    B_value_par	+=	B_value[i+ii][j+jj]*weight_square[ii][jj];
					    psi_par		+=	    psi[i+ii][j+jj]*weight_square[ii][jj];
				    }

			    if(NUMERICAL_EQUILIBRIUM)
                {
                    double g_eq_par = 0.0;
                    for(int ii=0;ii<2;ii++)
			            for(int jj=0;jj<2;jj++)
				            g_eq_par		+=	g_eq[i+ii][j+jj]*weight_square[ii][jj];

                    B_par[2] = g_eq_par/marker[p].X[0];
                }
                else
                    B_par[2] = B0*R0/marker[p].X[0];


			    {
				    COS_NPHI_PLUS_OMEGA_T = cos(n[0]*marker[p].X[2]-omega_A[0]*t);
				    SIN_NPHI_PLUS_OMEGA_T = sin(n[0]*marker[p].X[2]-omega_A[0]*t);


				    for(int I=i;I<=i+1;I++)
					    for(int J=j;J<=j+1;J++)
					    {
						    deltaA_par_local[I-i][J-j]  =	X[0]*deltaA_par_cos[0][I][J]*COS_NPHI_PLUS_OMEGA_T - X[0]*deltaA_par_sin[0][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[0]*deltaA_par_sin[0][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[0]*deltaA_par_cos[0][I][J]*SIN_NPHI_PLUS_OMEGA_T;
						    deltaPhi_local[I-i][J-j]    =	X[0]*deltaPhi_cos[0][I][J]*COS_NPHI_PLUS_OMEGA_T - X[0]*deltaPhi_sin[0][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[0]*deltaPhi_sin[0][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[0]*deltaPhi_cos[0][I][J]*SIN_NPHI_PLUS_OMEGA_T;
					    }
		

				    deltaA_par_par = 0.0;
				    deltaPhi_par = 0.0;

				
				    for(int ii=0;ii<2;ii++)
					    for(int jj=0;jj<2;jj++)
					    {
						    deltaA_par_par	+= deltaA_par_local[ii][jj]*weight_square[ii][jj];
						    deltaPhi_par	+=   deltaPhi_local[ii][jj]*weight_square[ii][jj];
					    }
			    }

			    R = marker[p].X[0];
			    Z = marker[p].X[1];

                double coeff_a,coeff_b,coeff_c;

			    if(SINGLE_Lambda < 1e-8)
                    coeff_a = -n[0]/2.0/omega_A[0]*species[0].mass;
                else
                    coeff_a = -n[0]/2.0/omega_A[0]*species[0].mass*(1+1/(B0/B_value_par/SINGLE_Lambda - 1));
			    coeff_b = species[0].mass*marker[p].X[0]*B_par[2]/B_value_par;
                coeff_c = species[0].charge*(R*deltaA_par_par)/alpha*B_par[2]/B_value_par + species[0].charge*psi_par/alpha - n[0]/omega_A[0]*species[0].charge*deltaPhi_par - K0;
			    if(SINGLE_v_par >0)
				    marker[p].v_par = (-coeff_b + sqrt(coeff_b*coeff_b - 4*coeff_a*coeff_c))/2/coeff_a;
			    else
				    marker[p].v_par = (-coeff_b - sqrt(coeff_b*coeff_b - 4*coeff_a*coeff_c))/2/coeff_a;

			    if(SINGLE_Lambda < 1e-8)
                    marker[p].v_per = 0.0;
                else
                    marker[p].v_per = abs(marker[p].v_par)/sqrt(B0/B_value_par/SINGLE_Lambda - 1);

			    marker[p].mu	= 0.5*species[0].mass*pow(marker[p].v_per,2.0)/B_value_par;

			    R	  = marker[p].X[0];
			    e     = species[0].charge;
			    v_par = marker[p].v_par;
			    m     = species[0].mass;
			    W_par = 0.5*species[0].mass*pow(marker[p].v_par,2.0);
			    W_per = 0.5*species[0].mass*pow(marker[p].v_per,2.0);
			
			    double K00 =  e*(R*deltaA_par_par)*B_par[2]/B_value_par/alpha + e*psi_par/alpha + v_par*coeff_b - n[0]/omega_A[0]*(W_per+W_par+e*deltaPhi_par);


			    cout << "marker[p].v_par  : " << marker[p].v_par/v_A  << "v_A  marker[p].v_per " << marker[p].v_per/v_A << "v_A marker[p].X[0]  : " << marker[p].X[0] << " Lambda :" << marker[p].mu*B0/(W_par+W_per) << " K00 : " << K00 << endl ;
			    res_index[p] = marker[p].id;
		    }
        }
    }
	MPI::COMM_WORLD.Bcast(res_index,N_p_diagnosis,MPI::INT,0);

}

/************************************************************************************************
*                                                                                               *
*                       Initial condition for test particle(MPI version)                        *
*                           only use for TEK(Time Evolution of KAM surfaces) case               *
*                                                                                               *
************************************************************************************************/
void initia_single_particle_MPI()
{
	double K0;
	double B_value_par,B_par[3],deltaPhi_par,deltaA_par_par;
	double weight_line_auxiliary[2][2];
	double weight_square[2][2];
	double phi;
	double e,v_par,m;
	int i,j,k;
	double psi_par,q,r2;
	double R,Z;
    double lambda;
    double W_per,W_par;
	//initialization fo #0 particle for main process!!!
    if(myid == 0)
    {
		marker[0].X[0]	= SINGLE_R;
		marker[0].X[1]	= SINGLE_Z;
		marker[0].X[2]	= SINGLE_phi;
		marker[0].v_par	= SINGLE_v_par;
		marker[0].v_per	= SINGLE_v_per;
		marker[0].id	= 0+myid*N_p_diagnosis;

		i = (int)((marker[0].X[0]-R0+a*a_b)/dR);
		j = (int)((marker[0].X[1]+a*a_b)/dZ);

		R   = R0-a*a_b+i*dR;
		Z   =   -a*a_b+j*dZ;

		weight_line_auxiliary[0][0] = (marker[0].X[0] - R)/dR;
		weight_line_auxiliary[1][0] = (marker[0].X[1] - Z)/dZ;

		weight_line_auxiliary[0][1] = 1.0 - weight_line_auxiliary[0][0];
		weight_line_auxiliary[1][1] = 1.0 - weight_line_auxiliary[1][0];
				

		for(int ii=0;ii<2;ii++)
			for(int jj=0;jj<2;jj++)
				weight_square[ii][jj] = weight_line_auxiliary[0][1-ii]*weight_line_auxiliary[1][1-jj];

		psi_par		= 0.0;
		B_value_par = 0.0;
				
		for(int ii=0;ii<2;ii++)
			for(int jj=0;jj<2;jj++)
			{
				B_value_par	+=	B_value[i+ii][j+jj]*weight_square[ii][jj];
				psi_par		+=	    psi[i+ii][j+jj]*weight_square[ii][jj];
			}

        if(NUMERICAL_EQUILIBRIUM)
        {
            double g_eq_par = 0.0;
            for(int ii=0;ii<2;ii++)
			    for(int jj=0;jj<2;jj++)
				    g_eq_par		+=	g_eq[i+ii][j+jj]*weight_square[ii][jj];

            B_par[2] = g_eq_par/marker[0].X[0];
        }
        else
            B_par[2] = B0*R0/marker[0].X[0];

/*----------------------------------------------------------------------------------------------------|
|                                                                                                   |
|                Specific Pitch angle Choosing : Lambda = mu*B0/E                                   |
|                   v_per = abs(v_par)/sqrt(B0/B-Lmabda-1), so (B0/B/Lambda-1)  > 0                 |
|                                                                                                   |
|----------------------------------------------------------------------------------------------------*/
        //if Lambda !=0
        if(SINGLE_Lambda > 1e-8)
        {
            if((B0/B_value_par/SINGLE_Lambda - 1) < 0)
            {
                cout << " (B0/B/Lambda-1) < 0 !!!!";
                exit(0);
            }
            else
                marker[0].v_per = abs(marker[0].v_par)/sqrt(B0/B_value_par/SINGLE_Lambda - 1);
        }
        else
            marker[0].v_per = 0.0;
		{
			COS_NPHI_PLUS_OMEGA_T = cos(n[0]*marker[0].X[2]-omega_A[0]*t);
			SIN_NPHI_PLUS_OMEGA_T = sin(n[0]*marker[0].X[2]-omega_A[0]*t);


			for(int I=i;I<=i+1;I++)
				for(int J=j;J<=j+1;J++)
				{
					deltaA_par_local[I-i][J-j]					=	X[0]*deltaA_par_cos[0][I][J]*COS_NPHI_PLUS_OMEGA_T - X[0]*deltaA_par_sin[0][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[0]*deltaA_par_sin[0][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[0]*deltaA_par_cos[0][I][J]*SIN_NPHI_PLUS_OMEGA_T;
					deltaPhi_local[I-i][J-j]					=	X[0]*deltaPhi_cos[0][I][J]*COS_NPHI_PLUS_OMEGA_T - X[0]*deltaPhi_sin[0][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[0]*deltaPhi_sin[0][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[0]*deltaPhi_cos[0][I][J]*SIN_NPHI_PLUS_OMEGA_T;
				}
		

			deltaA_par_par = 0.0;
			deltaPhi_par = 0.0;

				
			for(int ii=0;ii<2;ii++)
				for(int jj=0;jj<2;jj++)
				{
					deltaA_par_par	+= deltaA_par_local[ii][jj]*weight_square[ii][jj];
					deltaPhi_par	+= deltaPhi_local[ii][jj]*weight_square[ii][jj];
				}
		}


		//only for passing particles,so W_per = 0.0
		W_par = 0.5*species[0].mass*pow(marker[0].v_par,2.0);
		W_per = 0.5*species[0].mass*pow(marker[0].v_per,2.0);

		marker[0].mu	=	0.5*species[0].mass*pow(marker[0].v_per,2.0)/B_value_par;

		R	  = marker[0].X[0];
		e     = species[0].charge;
		v_par = marker[0].v_par;
		m     = species[0].mass;
		K0 =  e*(R*deltaA_par_par)/alpha*B_par[2]/B_value_par + e*psi_par/alpha + m*v_par*R*B_par[2]/B_value_par - n[0]/omega_A[0]*(W_per+W_par+e*deltaPhi_par);

        cout << "marker[0].v_par  : " << marker[0].v_par/v_A  << "v_A  marker[0].X[0]  : " << marker[0].X[0] << " K0 : " << K0 << " myid : " << myid << endl ;
			
		res_index[0] = marker[0].id;

    }

    MPI::COMM_WORLD.Bcast(&K0,1,MPI::DOUBLE,0);


	srand((unsigned)time(NULL)*(myid+1.0));
	for(int p=0;p<N_p_diagnosis;p++)
    {
        //if p=0,myid=0, particle's information do not need to be updated!!
        if(p==0 && myid ==0)
            continue;

		//define position
        //Random in R and phi at midplane(Z=0)
			
		marker[p].X[0] = SINGLE_R - DELTA_R + 2*DELTA_R*(rand()*1.0/RAND_MAX);
		marker[p].X[1] = 0.00;
        marker[p].X[2] = 2*PI*(rand()*1.0/RAND_MAX);
		marker[p].id   = p + myid*N_p_diagnosis;

		i = (int)((marker[p].X[0]-R0+a*a_b)/dR);
		j = (int)((marker[p].X[1]+a*a_b)/dZ);


		R   = R0-a*a_b+i*dR;
		Z   =   -a*a_b+j*dZ;

		weight_line_auxiliary[0][0] = (marker[p].X[0] - R)/dR;
		weight_line_auxiliary[1][0] = (marker[p].X[1] - Z)/dZ;

		weight_line_auxiliary[0][1] = 1.0 - weight_line_auxiliary[0][0];
		weight_line_auxiliary[1][1] = 1.0 - weight_line_auxiliary[1][0];
				

		for(int ii=0;ii<2;ii++)
			for(int jj=0;jj<2;jj++)
				weight_square[ii][jj] = weight_line_auxiliary[0][1-ii]*weight_line_auxiliary[1][1-jj];

		psi_par		= 0.0;
		B_value_par = 0.0;
				
		for(int ii=0;ii<2;ii++)
			for(int jj=0;jj<2;jj++)
			{
				B_value_par	+=	B_value[i+ii][j+jj]*weight_square[ii][jj];
				psi_par		+=	    psi[i+ii][j+jj]*weight_square[ii][jj];
			}

		if(NUMERICAL_EQUILIBRIUM)
        {
            double g_eq_par = 0.0;
            for(int ii=0;ii<2;ii++)
			    for(int jj=0;jj<2;jj++)
				    g_eq_par		+=	g_eq[i+ii][j+jj]*weight_square[ii][jj];

            B_par[2] = g_eq_par/marker[p].X[0];
        }
        else
            B_par[2] = B0*R0/marker[p].X[0];


		{
			COS_NPHI_PLUS_OMEGA_T = cos(n[0]*marker[p].X[2]-omega_A[0]*t);
			SIN_NPHI_PLUS_OMEGA_T = sin(n[0]*marker[p].X[2]-omega_A[0]*t);


			for(int I=i;I<=i+1;I++)
				for(int J=j;J<=j+1;J++)
				{
					deltaA_par_local[I-i][J-j]  =	X[0]*deltaA_par_cos[0][I][J]*COS_NPHI_PLUS_OMEGA_T - X[0]*deltaA_par_sin[0][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[0]*deltaA_par_sin[0][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[0]*deltaA_par_cos[0][I][J]*SIN_NPHI_PLUS_OMEGA_T;
					deltaPhi_local[I-i][J-j]    =	X[0]*deltaPhi_cos[0][I][J]*COS_NPHI_PLUS_OMEGA_T - X[0]*deltaPhi_sin[0][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[0]*deltaPhi_sin[0][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[0]*deltaPhi_cos[0][I][J]*SIN_NPHI_PLUS_OMEGA_T;
				}
		

			deltaA_par_par = 0.0;
			deltaPhi_par = 0.0;

				
			for(int ii=0;ii<2;ii++)
				for(int jj=0;jj<2;jj++)
				{
					deltaA_par_par	+= deltaA_par_local[ii][jj]*weight_square[ii][jj];
					deltaPhi_par	+=   deltaPhi_local[ii][jj]*weight_square[ii][jj];
				}
		}

		R = marker[p].X[0];
		Z = marker[p].X[1];

        double coeff_a,coeff_b,coeff_c;

		if(SINGLE_Lambda < 1e-8)
            coeff_a = -n[0]/2.0/omega_A[0]*species[0].mass;
        else
            coeff_a = -n[0]/2.0/omega_A[0]*species[0].mass*(1+1/(B0/B_value_par/SINGLE_Lambda - 1));
		coeff_b = species[0].mass*marker[p].X[0]*B_par[2]/B_value_par;
        coeff_c = species[0].charge*(R*deltaA_par_par)/alpha*B_par[2]/B_value_par + species[0].charge*psi_par/alpha - n[0]/omega_A[0]*species[0].charge*deltaPhi_par - K0;
		if(SINGLE_v_par >0)
			marker[p].v_par = (-coeff_b + sqrt(coeff_b*coeff_b - 4*coeff_a*coeff_c))/2/coeff_a;
		else
			marker[p].v_par = (-coeff_b - sqrt(coeff_b*coeff_b - 4*coeff_a*coeff_c))/2/coeff_a;

		if(SINGLE_Lambda < 1e-8)
            marker[p].v_per = 0.0;
        else
            marker[p].v_per = abs(marker[p].v_par)/sqrt(B0/B_value_par/SINGLE_Lambda - 1);

		marker[p].mu	= 0.5*species[0].mass*pow(marker[p].v_per,2.0)/B_value_par;

		R	  = marker[p].X[0];
		e     = species[0].charge;
		v_par = marker[p].v_par;
		m     = species[0].mass;
		W_par = 0.5*species[0].mass*pow(marker[p].v_par,2.0);
		W_per = 0.5*species[0].mass*pow(marker[p].v_per,2.0);
			
		double K00 =  e*(R*deltaA_par_par)*B_par[2]/B_value_par/alpha + e*psi_par/alpha + v_par*coeff_b - n[0]/omega_A[0]*(W_per+W_par+e*deltaPhi_par);


//		cout << "marker[p].v_par  : " << marker[p].v_par/v_A  << "v_A  marker[p].v_per " << marker[p].v_per/v_A << "v_A marker[p].X[0]  : " << marker[p].X[0] << " Lambda :" << marker[p].mu*B0/(W_par+W_per) << " K00 : " << K00 << " id = " << marker[p].id <<" myid = " << myid << endl;
			    
        res_index[p] = marker[p].id;

    }
}

void initia_parameter_output()
{

/************************************************************************************************
*                                                                                               *
*               Output initia maximum field quantity (E_R_max,B_R_max)                          *
*                               E_R_max = E_cos_max(t=0)                                        *
*                               B_R_max = B_sin_max(t=0)                                        *
*                                                                                               *
************************************************************************************************/
	
	//output
	ofstream file_par("Tokamak_parameter.dat",ios::app);
    ofstream file_par_mul("Tokamak_parameter_multiple_mode.dat",ios::app);
    ofstream file_par_of("Tokamak_parameter_Calc_Orbit_Frequency.dat",ios::app);
	ofstream file_output("parameter_output",ios::app);
    ofstream file_B_R_max("deltaB_R_max.dat",ios::app);
	double deltaB_R_max = 0.0;
	double deltaE_R_max = 0.0;
    double deltaB_3D_R_max[NUM_MODE] = {0.0};
	if(myid==0)
	{
		for(int i=0;i<mR;i++)
			for(int j=0;j<mZ;j++)
				for(int k=0;k<mymphi;k++)
				{
					double deltaB_R = 0.0;
					double deltaE_R = abs(deltaE_cos_3D[0][0][i][j][k]);
					if(deltaB_R > deltaB_R_max)
						deltaB_R_max = deltaB_R;
					if(deltaE_R > deltaE_R_max)
						deltaE_R_max = deltaE_R;
				}
        //output midplane maxima of B_R at phi = 0 surface!
        for(int s=0;s<NUM_MODE;s++)
            for(int i=0;i<mR;i++)
                deltaB_3D_R_max[s] = abs(deltaB_initial_3D[s][0][i][int(mZ/2)][0])>deltaB_3D_R_max[s] ? (abs(deltaB_initial_3D[s][0][i][int(mZ/2)][0])) : deltaB_3D_R_max[s];

        for(int s=0;s<NUM_MODE;s++)
            file_B_R_max << n[s] << " " << deltaB_3D_R_max[s] << endl;

	//output parameter:R0    a   B0    omega_A    wave amplitude  n deltaB_R_max deltaE_R_max v_A
		file_par << setiosflags(ios::fixed) << setprecision(8) << R0 << " " << a << " " << B0 << " " << omega_A[0] <<  " " << sqrt(X[0]*X[0]+Y[0]*Y[0]) << " " << n[0] << " " << deltaB_R_max << " " << deltaE_R_max << " " << v_A << " " << nu/omega_Alfen << endl;
        file_par_of << NUM_PPHI << " " << NUM_E << " " << NUM_TRACK_PER_PPHI*NUM_TRACK_PER_E << " " << LAMBDA << endl;
        for(int s=0;s<NUM_MODE;s++)
            file_par_mul << setiosflags(ios::fixed) << setprecision(8) << n[s] << " " << omega_A[s] << endl;
	//	output parameter
		file_output << "-------------------------------parameter list-------------------------------" << endl;
		file_output << "Number of particles                                             :   " << numprocs*N_p_initial                           << endl;
		file_output << "omega_c/omega_Alfen                                             :   " << sqrt(alpha1)                                   << endl;
	if(NUM_MODE == 1)
    {
        file_output << "mode frequency(omega/omega_Alfen)                               :   " << omega_A[0]/omega_Alfen                            << endl;
   		file_output << "toridal number                                                  :   " << n[0]                                              << endl;
    }
    else
    {
    for(int s=0;s<NUM_MODE;s++)
    {
        file_output << "mode #"<<s+1<<" frequency(omega/omega_Alfen)                            :   " << omega_A[s]/omega_Alfen                            << endl;
        file_output << "toridal number #"<<s+1<<"                                               :   " << n[s]                                              << endl;
    }
    }
        file_output << "aspect ratio (R/a)                                              :   " << R0/a                                           << endl;
		file_output << "q(0)                                                            :   " << q_0                                            << endl;
		file_output << "q(a)                                                            :   " << q_a                                            << endl;
		file_output << "max lamor radius (rho_h/a)                                      :   " << species[0].mass*v_0/species[0].charge/B0/a     << endl;
		file_output << "volume average fast ion beta                                    :   " << beta                                           << endl;
    if(NUM_MODE == 1)
		file_output << "damping rate (gamma_d/omega_Alfen)                              :   " << gamma_d[0]/omega_Alfen                         << endl;
    else
    {
    for(int s=0;s<NUM_MODE;s++)
        file_output << "damping rate (gamma_d/omega_Alfen) for mode # "<<s+1<<"                 :   " << gamma_d[s]/omega_Alfen                         << endl;
    }
        file_output << "darg rate(nu/omega_Alfen)                                       :   " << nu/omega_Alfen                                 << endl;
    if(SCATTERING)
		file_output << "pitch angle scattering                                          :   " << "included"                                     << endl;
    else
        file_output << "pitch angle scattering                                          :   " << "no"                                           << endl;
    if(SLOWING_DOWN)
		file_output << "slowing down process                                            :   " << "included"                                     << endl;
    else
        file_output << "slowing down process                                            :   " << "no"                                           << endl;
    if(SOURCE_AND_SINK)
		file_output << "source and sink                                                 :   " << "yes"                                                  << endl;
	else
		file_output << "source and sink                                                 :   " << "no"                                                   << endl;
    if(!SPL)
    {
    if(ISOTROPY) 
		if(!FGP)
			file_output << "distrubution in velocity                                        :   " << "isotropy"                                     << endl;
		else
        {
            if(!DFGP)
			file_output << "distrubution in velocity                                        :   " << "isotropy(Finite width of Gaussian distrbution of Pitch angle: " << Lambda_0_FGP << ")" << endl;
            else
            {
                file_output << "distrubution in velocity                                        :   " << "isotropy(Double Finite width of Gaussian distrbution of Pitch angle #1 " << Lambda_1_FGP << "), proportion #1 : " << c1_FGP << endl;
                file_output << "distrubution in velocity                                        :   " << "isotropy(Double Finite width of Gaussian distrbution of Pitch angle #2 " << Lambda_2_FGP << "), proportion #2 : " << c2_FGP << endl;
            }
        }
	else
		file_output << "distrubution in velocity                                        :   " << "anisotropic"                                  << endl;
    }
    else
    {
        file_output << "distrubution in velocity                                        :   " << "anisotropic(specific pitch angle #1)  Lambda = " << Lambda_0_load_1 << " , Number of particles with this pitch angle :" << numprocs*N_Lambada_1 << endl;
        file_output << "                                                                :   " << "anisotropic(specific pitch angle #2)  Lambda = " << Lambda_0_load_2 << " , Number of particles with this pitch angle :" << numprocs*N_Lambada_2 << endl;
    }
    
    if(UNIFORM_RADIAL)
		file_output << "TAE mode structrue uniform in                                   :   " << "r"                                            << endl;
	else
		file_output << "TAE mode structrue uniform in                                   :   " << "sqrt(psi)"                                    << endl;
	if(AVG)
		file_output << "how many wave period to average                                 :   " << WAVE_PERIOD_NUM                                << endl;
	else
		file_output << "how many wave period to average                                 :   " << "no average"                                   << endl;
	if(UNIFORM_LOADING)
        file_output << "particle loading in velocity sapce                              :   " << "uniform"                                   << endl;
    else
        file_output << "particle loading in velocity sapce                              :   " << "non-uniform"                                   << endl;
    if(METHOD == 1)
		file_output << "method used                                                     :   " << "Method 1 : calc current on grid"              << endl;
	else if(METHOD == 2)
		file_output << "method used                                                     :   " << "Method 2 : sum all particles"                 << endl;
	if(DELTAF_METHOD == 1)
		file_output << "delta-f method used                                             :   " << "1. dw/dt = - (f/g - w)dlnf_0/dt"              << endl;
	else if(DELTAF_METHOD == 2)
		file_output << "delta-f method used                                             :   " << "2. dw/dt = - 1/g*df_0/dt"                     << endl;
	if(WEIGHT_CONDITION == 1)
		file_output << "assumption for weight equation                                  :   " << "1. no any assumption"                         << endl;
	else if(WEIGHT_CONDITION == 2)
		file_output << "assumption for weight equation                                  :   " << "2. dp_phidt = n/omega_A*dEdt"                 << endl;
	else if(WEIGHT_CONDITION == 3)
		file_output << "assumption for weight equation                                  :   " << "3. dEdt     = omega_A/n*dp_phidt"             << endl;
	if(A_PAR_INT)
		file_output << "method of geting A_par                                          :   " << "1. A_par = int grad_|| deltaPhi dt"           << endl;
	else
		file_output << "method of geting A_par                                          :   " << "2. A_par = 1/omega_A/B*(-J^{-1}*m+B_phi/R*n * deltaPhi" << endl;
	if(ADIABATIC == 1)
		file_output << "adiabatic term                                                  :   " << "eliminate adiabatic term using candy's method"        << endl;
	else if(ADIABATIC == 2)
		file_output << "adiabatic term                                                  :   " << "eliminate adiabatic term using fu's method"           << endl;
	else if(ADIABATIC == 0)
		file_output << "adiabatic term                                                  :   " << "with adiabatic term"                                  << endl;
	if(RECYCLE_METHOD == 1)
		file_output << "treat bounary particles:                                        :   " << "recycle"                                              << endl;
	else if(RECYCLE_METHOD == 2)
		file_output << "treat bounary particles:                                        :   " << "remove"                                               << endl;
	if(SHAFRANOV)
		file_output << "including Shafranov shift:                                      :   " << "true"                                                 << endl;
	else
		file_output << "including Shafranov shift:                                      :   " << "false"                                                << endl;
	if(AVERAGE_PSI)
		file_output << "different orbit averaging method for two species particles      :   " << "true"                                                 << endl;
	else
		file_output << "different orbit averaging method for two species particles      :   " << "false"                                                << endl;
	if(FLUX_COORDINATE == 1)
		file_output << "flux coordinate type                                            :   " << "theta_f = theta_S - 2*eta*sin(theta_S),zeta_f = phi"  << endl;
	else if(FLUX_COORDINATE == 2)
		file_output << "flux coordinate type                                            :   " << "theta_f = theta_S ,zeta_f = phi-q*delta"              << endl;
	if(JACOBIAN_ACCURACY)
		file_output << "Jocobian approxmation is used                                   :   " << "no"                                                   << endl;
	else
		file_output << "Jocobian approxmation is used                                   :   " << "yes"                                                  << endl;
	if(SHIFT_ACCURACY)
		file_output << "Shafranov shift has q-profile                                   :   " << "real"                                                 << endl;
	else
		file_output << "Shafranov shift has q-profile                                   :   " << "constant"                                             << endl;
    if(NUMERICAL_EQUILIBRIUM)
		file_output << "equilibrium type                                                :   " << "numerical equilibrium"                                << endl;
	else
		file_output << "equilibrium type                                                :   " << "analytical equilibrium"                               << endl;
    if(IMC)
		file_output << "EP current J_h Including Magnetization Current                  :   " << "Yes"                                                  << endl;
	else
		file_output << "EP current J_h Including Magnetization Current                  :   " << "No"                                                   << endl;
    if(COF)
    {
        file_output << endl;
        file_output << "-----------parameter for calculation of particle orbit frequency------------"                                                      << endl;
        if(OPS == 1)
        {
		file_output << "grid on P_phi                                                   :   " << NUM_PPHI                                                 << endl;
        file_output << "grid on Energy                                                  :   " << NUM_E                                                    << endl;
        file_output << "pitch angle variable Lambda=mu*B0/E                             :   " << LAMBDA                                                   << endl;
        file_output << "how many track of particles for diagnostics                     :   " << NUM_TRACK_PER_PPHI*NUM_TRACK_PER_E                       << endl;
        }
        else if (OPS == 2)
        {
        file_output << "grid on Energy                                                  :   " << NUM_E                                                   << endl;
        file_output << "grid on pitch angle variable(Lambda=mu*B0/E)                    :   " << NUM_LAMBDA                                              << endl;
        file_output << "fixed value of P_phi                                            :   " << P_phi_0                                                 << endl;
        file_output << "how many track of particles for diagnostics                     :   " << NUM_TRACK_PER_E*NUM_TRACK_PER_LAMBDA                    << endl;
        }
        if(CO_PASSING)
        file_output << "sign for parallel velocity                                      :   " << "co-passing (+)"                                         << endl;
        else
        file_output << "sign for parallel velocity                                      :   " << "counter-passing (-)"                                    << endl;
        if(CPN)
        file_output << "amplitude(calculate orbit frequency including Nonlinear terms  ):   " << X[0]                                                     << endl;
    }


    if(EDT)
    {
        file_output << endl;
        file_output << "--------------parameter for evolution of distribution function--------------" << endl;
		file_output << "grid on P_phi                                                   :   " << GRID_P_PHI                                               << endl;
        file_output << "grid on Energy                                                  :   " << GRID_E                                                   << endl;
        file_output << "pitch angle variable Lambda0=mu*B0/E for diagnostic             :   " << Lambda_0                                                 << endl;
        file_output << "Energy E0=0.5*m*(v_||)^2+mu*B for diagnostic(/v_A^2)            :   " << E_0/v_A/v_A                                              << endl;
        file_output << "how many time step to record distribution function data         :   " << TIME_STEP_INTERVAL                                       << endl;
        file_output << "size of diagnostic windows for fixed Lambda and E               :   " << windows_width                                            << endl;
        if(DPE && !DPK && !IAL)
        file_output << "output distribution function as function of                     :   " << "f(P_phi,E,Lambda_0)"                                    << endl;
        if(DPE && DPK  && !IAL)
        file_output << "output distribution function as function of                     :   " << "f(P_phi,K,Lambda_0)"                                    << endl;
        if(DPE && !DPK && IAL)
        file_output << "output distribution function as function of                     :   " << "f(P_phi,E) = intergrate f(P_phi,E,Lambda)dLambda "      << endl;
        
    }

    if(DIAGNOSTICS)
    {
        file_output << endl;
        file_output << "--------------parameter for diagnostics of deltaf in 2D phase space--------------" << endl;
        file_output << "time interval region of total simulation for diagnostics        :   " << NUM_TIME_DIAGNOSTICS                                     << endl;
        if(dLP)
        {
        file_output << "output deltaf in 2D phase space                                 :   " << "deltaf(Lambda,P_phi)"                                        << endl;
        file_output << "Energy E0=0.5*m*(v_||)^2+mu*B for diagnostic(/v_A^2)            :   " << E_0_dLP/v_A/v_A                                          << endl;
        file_output << "window width for Energy E0                                      :   " << windows_width_dLP                                          << endl;
        }
        else if(dLE)
        {
        file_output << "output deltaf in 2D phase space                                 :   " << "deltaf(Lambda,E)"                                        << endl;
        file_output << "P_phi0                                                          :   " << P_phi_0_dLE                                                  << endl;
        file_output << "window width for Energy P_phi0                                  :   " << windows_width_dLE                                          << endl;
        }
        else if(dEP)
		if(mIL)
		file_output << "output deltaf in 2D phase space with fixed mu                   :   " << "deltaf(Energy,P_phi,mu=constant)"                                        << endl;
		else
        file_output << "output deltaf in 2D phase space with fixed Lambda               :   " << "deltaf(Energy,P_phi,Lambda=constant)"                                        << endl;
	}


    }

    file_par.close();
    file_par_mul.close();
    file_output.close();
}


/********************************************************************************************************
*                                                                                                       *
*               4th-order  Runge-Kutta Method for push particle                                         *
*   equation of motion  : dX/dt     =  1/B^star_||*[v_||*B^star + b_0 times (1/q*mu*grad(B) - deletaE)] *
*                       : dv_||/dt  =  1/B^star_||/m*B^star*(q*deltaE-mu*grad(B))                       *
*   weight equation     : dwdt      = -1/g*(dP_phi/dt*df_0/dP_phi+dE/dt*df_0/dE)                        *
*                                                                                                       *
********************************************************************************************************/

void stepon_RK4()
{
	particle_vector<double> B_par,b_par,deltaB_par,deltaE_par,curl_b_par,grad_B_par;
	double B_value_par,psi_par;
	double deltaPhi_par,deltaA_par_par;
	particle_vector<double> grad_deltaA_par_par,grad_deltaPhi_par;
	double weight_square[2][2];
	double weight_line[3][2];
	double weight_line_aux[2];
	double weight_line_auxiliary[2][2];

	double R,Z,phi;
	int marker_id=0;

	particle_vector<double> dXdt1,v_d,grad_p_phi;
	double dv_pardt1,dEdt,df0dt,dp_phidt,dLambdadt;
	double df_0dP_phi,df_0dE,df_0dLambda;
	double mu,E,v,p_phi;
    double W_per,W_par;
	df0dt = 0.0;

	particle_vector<double> B_s;
	double B_ss;

//	double E_0=0.5*v_A*v_A;
//	double mu_0 = 0.1*E_0;

	double m = species[0].mass;
	double e = species[0].charge;


    ofstream file_orbit("particle_orbit.dat",ios::app);
    ofstream file_p_phi("p_phi.dat",ios::app);
    ofstream file_E_par("E_par.dat",ios::app);
    ofstream file_B_par("B_par.dat",ios::app);
	ofstream file_KAM("KAM.dat",ios::app);
    ofstream file_KAM_orbit("KAM_orbit.dat",ios::app);

//    //output KAM surfaces vs. time
//    ofstream file_KAM_t;



    TIMESTEP = 0;
    dZ_of_step1 = 0;
    Z_reverse_sign = false;

    //variables related to RK4
	double timestep,coeff;
	double Va,RHS_V,Vold;
	double Wa,RHS_W,Wold;
    //slowing down process
    double Ua,RHS_U,Uold,lambda;
	particle_vector<double> Xa,RHS_X,Xold;
	particle_vector<double> deltaE_X_par[NUM_MODE],deltaE_Y_par[NUM_MODE],grad_deltaPhi_X_par,grad_deltaPhi_Y_par;
	double told;
	int i,j;
	bool OUT_OF_BOUDARY;

    //variables related to J_M
    //add by 2014.05.20
    particle_vector<double> ddeltaB_Xdt_par[NUM_MODE],ddeltaB_Ydt_par[NUM_MODE],dgrad_deltaA_par_Xdt_par,dgrad_deltaA_par_Ydt_par;
    double ddeltaA_par_Xdt_par,ddeltaA_par_Ydt_par;
    double dgrad_deltaA_par_Xdt_local[3][2][2],dgrad_deltaA_par_Ydt_local[3][2][2],ddeltaA_par_Xdt_local[2][2],ddeltaA_par_Ydt_local[2][2];

    double J_M_dot_E_X[NUM_MODE],J_M_dot_E_Y[NUM_MODE];

	while(t < total_time)
	{

        TIMESTEP ++;

/************************************************************************************************
*                                                                                               *
*             output Time Evolution of KAM surfaces(test particle simulation)                   *
*                amplitude and phase should be input form file amplitude.dat                    *
*                                                                                               *
************************************************************************************************/
        if(TEK)
        {
            X[0] = X_input[TIMESTEP-1];
            Y[0] = Y_input[TIMESTEP-1];
        }

/************************************************************************************************
*                                                                                               *
*           Calculate toroidal precession frequency : omega_phi   = Delta_phi/Delta_t           *
*                     poloidal transit    frequency : omega_theta = 2*pi/Delta_t                *
*                                                                                               *
************************************************************************************************/

        if(COF)
        {
            if(calc_oribt_frequency(LOSS,OUTPUT_TRACK))
            {
                cout << "finish calculation of particle orbit frequency!!!";
                break;
            }

/************************************************************************************************
*                                                                                               *
*                          Calculate Phase including Nonlinear terms                            *
*                    amplitude and phase should be input form file amplitude.dat                *
*                                                                                               *
************************************************************************************************/

            if(CPN)
            {
                if(index_t_start_CPN+TIMESTEP-1 >= NUM_LINE_AMP)
                {
                    cout << "particle orbit period in theta direction is longer than input simulation time!!!";
                    exit(0);
                }
                X[0] = X_input[index_t_start_CPN+TIMESTEP-1];
                Y[0] = Y_input[index_t_start_CPN+TIMESTEP-1];
                if(TIMESTEP == 1)
                {
                    alpha_initial = atan2(Y[0],X[0]); 
                    alpha_previous = alpha_initial;
                    Delta_alpha    = 0.0;
                }
                else
                {
                    // region of atan2 is [-pi pi],so if alpha is larger(smaller) than pi(-pi), alpha should be decreased(increased) by pi.
                    // actually we decrease(increase) the alpha_current by pi when Delta alpha is more than 80% of 2*pi.
                    alpha_current  =  atan2(Y[0],X[0]);
                    if(alpha_current > 0 && alpha_previous < 0 && alpha_current - alpha_previous >=  0.8*2*PI)
                        Delta_alpha += (alpha_current - 2*PI - alpha_previous);
                    else if(alpha_current < 0 && alpha_previous > 0 && alpha_current - alpha_previous <=  -0.8*2*PI)
                        Delta_alpha += (alpha_current + 2*PI - alpha_previous);
                    else
                        Delta_alpha += (alpha_current - alpha_previous);
                    alpha_previous =  alpha_current;
                }
            }
        }

/*************************************************************************************************
*                                                                                                *
*                                 diagnosis for resonant particles                               *
*                                                                                                *
*************************************************************************************************/

		if(RES && !TEK)
		{
            /*
            //output time evolution of KAM surface and define file name, add by 2013-10-26
            if(TEK && (int((t-dt)/dt))%(int(1.0/NUM_TIME_DIAGNOSTICS*(total_time/dt))) == 0)
            {
                string fh;
	            stringstream ss,tt;
                int index_of_time_record = ((int((t-dt)/dt))/(int(1.0/NUM_TIME_DIAGNOSTICS*(total_time/dt))));
	            ss << index_of_time_record;
                tt << setiosflags(ios::fixed) << setprecision(2) << SINGLE_Lambda;
	            fh = "KAM_vs_Time_" + ss.str() + "_Lambda=" + tt.str() + ".dat";
	            char *filename = (char*)fh.c_str();
	            file_KAM_t.open(filename,ios::app);
            }
            */

			if(((int)(t/dt))%RES_STEP == 0)
            {
				for(int p=0;p<N_p;p++)
				{
					for(int l=0;l<N_p_diagnosis;l++)
					{
						if(marker[p].id == res_index[l])
                        {
							int i = (int)((marker[p].X[0]-R0+a*a_b)/dR);
							int j = (int)((marker[p].X[1]+a*a_b)/dZ);

							R   = R0-a*a_b+i*dR;
							Z   =   -a*a_b+j*dZ;

							weight_line_auxiliary[0][0] = (marker[p].X[0] - R)/dR;
							weight_line_auxiliary[1][0] = (marker[p].X[1] - Z)/dZ;

							weight_line_auxiliary[0][1] = 1.0 - weight_line_auxiliary[0][0];
							weight_line_auxiliary[1][1] = 1.0 - weight_line_auxiliary[1][0];
				

							for(int ii=0;ii<2;ii++)
								for(int jj=0;jj<2;jj++)
									weight_square[ii][jj] = weight_line_auxiliary[0][1-ii]*weight_line_auxiliary[1][1-jj];

							psi_par		= 0.0;
							B_value_par = 0.0;
				
							for(int ii=0;ii<2;ii++)
								for(int jj=0;jj<2;jj++)
								{
									B_value_par	+=	B_value[i+ii][j+jj]*weight_square[ii][jj];
									psi_par		+=		psi[i+ii][j+jj]*weight_square[ii][jj];
								}

							if(NUMERICAL_EQUILIBRIUM)
                            {
                                double g_eq_par = 0.0;
                                for(int ii=0;ii<2;ii++)
			                        for(int jj=0;jj<2;jj++)
				                        g_eq_par		+=	g_eq[i+ii][j+jj]*weight_square[ii][jj];

                                B_par[2] = g_eq_par/marker[p].X[0];
                            }
                            else
                                B_par[2] = B0*R0/marker[p].X[0];


							{
								COS_NPHI_PLUS_OMEGA_T = cos(n[0]*marker[p].X[2]-omega_A[0]*t);
								SIN_NPHI_PLUS_OMEGA_T = sin(n[0]*marker[p].X[2]-omega_A[0]*t);


								for(int I=i;I<=i+1;I++)
									for(int J=j;J<=j+1;J++)
									{
										deltaA_par_local[I-i][J-j]	=	X[0]*deltaA_par_cos[0][I][J]*COS_NPHI_PLUS_OMEGA_T - X[0]*deltaA_par_sin[0][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[0]*deltaA_par_sin[0][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[0]*deltaA_par_cos[0][I][J]*SIN_NPHI_PLUS_OMEGA_T;
										deltaPhi_local[I-i][J-j]	=	X[0]*deltaPhi_cos[0][I][J]*COS_NPHI_PLUS_OMEGA_T - X[0]*deltaPhi_sin[0][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[0]*deltaPhi_sin[0][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[0]*deltaPhi_cos[0][I][J]*SIN_NPHI_PLUS_OMEGA_T;
									}
		

								deltaA_par_par = 0.0;
								deltaPhi_par = 0.0;

				
								for(int ii=0;ii<2;ii++)
									for(int jj=0;jj<2;jj++)
									{
										deltaA_par_par	+= deltaA_par_local[ii][jj]*weight_square[ii][jj];
										deltaPhi_par	+= deltaPhi_local[ii][jj]*weight_square[ii][jj];
									}
							}


							R = marker[p].X[0];
							Z = marker[p].X[1];

 
							W_per = marker[p].mu*B_value_par;
							W_par = 0.5*m*pow(marker[p].v_par,2.0);

							double p_phi_unperturbed =	e*psi_par/alpha + m*marker[p].v_par*marker[p].X[0]*B_par[2]/B_value_par;
							double K =  e*(marker[p].X[0]*deltaA_par_par)/alpha*B_par[2]/B_value_par + e*psi_par/alpha + m*marker[p].v_par*marker[p].X[0]*B_par[2]/B_value_par - n[0]/omega_A[0]*(W_per+W_par+e*deltaPhi_par);

/************************************************************************************
*                                                                                   *
*       restrictions for particles output : Z = 0(theta = 0)                        *
*                                           K = n/omega*E+P_phi = const             *
*                        plot KAM surface : PHASE = n*phi -omega*t                  *
*                                                                                   *
************************************************************************************/
							if((marker[p].X[1] >0 && marker[p].w <0 && marker[p].X[0] > R0) || (marker[p].X[1] <0 && marker[p].w >0 && marker[p].X[0] > R0))
							{
								double phase = n[0]*marker[p].X[2]-omega_A[0]*t;
								while(phase>PI)
									phase = phase - 2*PI;
								while(phase<-PI)
									phase = phase + 2*PI;
//                                if(!TEK)//if output time evolution of KAM, do not write KAM.dat!!!!
                                    file_KAM << setiosflags(ios::fixed) << setprecision(12) << marker[p].X[0] << "	" << p_phi_unperturbed << "	" << marker[p].v_par << "	"  << phase << "	" << setiosflags(ios::fixed) << setprecision(6) << marker[p].id << " " << t << endl;
/************************************************************************************
*                         output KAM surface v.s. time (2013.10.26)                 *
************************************************************************************/
//                                if(TEK && (int((t-dt)/dt))%(int(1.0/NUM_TIME_DIAGNOSTICS*(total_time/dt))) == 0)
//                                    file_KAM_t << setiosflags(ios::fixed) << setprecision(12) << marker[p].X[0] << "	" << p_phi_unperturbed << "	" << marker[p].v_par << "	"  << phase << "	" << setiosflags(ios::fixed) << setprecision(6) << marker[p].id << endl;
							}
							//update marker_diagnosis info
							//marker[p].w = marker[p].X[1];
/************************************************************************************
*                           output diagnostic particle's P_phi v.s time             *
************************************************************************************/
                            file_KAM_orbit << p_phi_unperturbed << " " << marker[p].id << " " ;
						}
                    }					
                }
            if(myid == 0)
                file_KAM_orbit << t << endl;
            }

//            //close KAM_T file, add by 2013-10-26
//            if(TEK && (int((t-dt)/dt))%(int(1.0/NUM_TIME_DIAGNOSTICS*(total_time/dt))) == 0)
//                file_KAM_t.close();
        }



/************************************************************************************
*                                                                                   *
*                           diagnosis #0 particle                                   *
*                                                                                   *
************************************************************************************/

		if(SPI)
        {
            for(int p=0;p<N_p;p++)
			    if(marker[p].id == 0)
			    {
                    if(TIMESTEP == 1)
                    {
                        ofstream file_SP_parameter("single_particle_parameter.dat",ios::app);
		                file_SP_parameter << marker[p].X[0] << " " << marker[p].X[1] << " " << marker[p].X[2] << " " << marker[p].v_par << " " << marker[p].v_per;
                        if(SPL)
                            file_SP_parameter << " " << Lambda_0_load_1;
                        file_SP_parameter << endl;
		                file_SP_parameter.close();
                    }

				    i = (int)((marker[p].X[0]-R0+a*a_b)/dR);
				    j = (int)((marker[p].X[1]+a*a_b)/dZ);

				    R   = R0-a*a_b+i*dR;
				    Z   =   -a*a_b+j*dZ;

                    weight_line_auxiliary[0][0] = (marker[p].X[0] - R)/dR;
				    weight_line_auxiliary[1][0] = (marker[p].X[1] - Z)/dZ;

				    weight_line_auxiliary[0][1] = 1.0 - weight_line_auxiliary[0][0];
				    weight_line_auxiliary[1][1] = 1.0 - weight_line_auxiliary[1][0];
				

				    for(int ii=0;ii<2;ii++)
					    for(int jj=0;jj<2;jj++)
						    weight_square[ii][jj] = weight_line_auxiliary[0][1-ii]*weight_line_auxiliary[1][1-jj];

				    B_par       = 0.0;
				    curl_b_par  = 0.0;
				    grad_B_par  = 0.0;
				    B_value_par = 0.0;
                    psi_par     = 0.0;
				
				    for(int ii=0;ii<2;ii++)
					    for(int jj=0;jj<2;jj++)
					    {
						    B_value_par		+=	B_value[i+ii][j+jj]*weight_square[ii][jj];
                            psi_par		    +=	    psi[i+ii][j+jj]*weight_square[ii][jj];
						    for(int d=0;d<2;d++)
						    {
							    B_par[d]      +=      B[d][i+ii][j+jj]*weight_square[ii][jj];
							    grad_B_par[d] += grad_B[d][i+ii][j+jj]*weight_square[ii][jj];
						    }

						    for(int d=0;d<3;d++)
							    curl_b_par[d] += curl_b[d][i+ii][j+jj]*weight_square[ii][jj];
					    }

				    if(NUMERICAL_EQUILIBRIUM)
                    {
                        double g_eq_par = 0.0;
                        for(int ii=0;ii<2;ii++)
			                for(int jj=0;jj<2;jj++)
				                g_eq_par		+=	g_eq[i+ii][j+jj]*weight_square[ii][jj];

                        B_par[2] = g_eq_par/marker[p].X[0];
                    }
                    else
                        B_par[2] = B0*R0/marker[p].X[0];

				    grad_B_par[2] = 0.0;


				    for(int d=0;d<3;d++)
					    b_par[d] =   B_par[d]/B_value_par;

				    grad_deltaA_par_par = 0.0;
				    grad_deltaPhi_par   = 0.0;
				    deltaA_par_par      = 0.0;
                    for(int s=0;s<NUM_MODE;s++)
				    {
					    COS_NPHI_PLUS_OMEGA_T = cos(n[s]*marker[p].X[2]-omega_A[s]*t);
					    SIN_NPHI_PLUS_OMEGA_T = sin(n[s]*marker[p].X[2]-omega_A[s]*t);

					    for(int d=0;d<3;d++)
						    for(int I=i;I<=i+1;I++)
							    for(int J=j;J<=j+1;J++)
							    {
								    if(d!=2)
									    grad_deltaA_par_local[d][I-i][J-j]  =       X[s]*grad_deltaA_par_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T - X[s]*grad_deltaA_par_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaA_par_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaA_par_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T;
								    else
									    grad_deltaA_par_local[d][I-i][J-j]  =   (-  X[s]*grad_deltaA_par_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T - X[s]*grad_deltaA_par_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T - Y[s]*grad_deltaA_par_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaA_par_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T)*n[s]/marker[p].X[0];
								    if(d!=2)
									    grad_deltaPhi_local[d][I-i][J-j]    =       X[s]*grad_deltaPhi_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T - X[s]*grad_deltaPhi_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaPhi_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaPhi_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T;
								    else
									    grad_deltaPhi_local[d][I-i][J-j]    =   (-  X[s]*grad_deltaPhi_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T - X[s]*grad_deltaPhi_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T - Y[s]*grad_deltaPhi_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaPhi_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T)*n[s]/marker[p].X[0];

							    }
					    for(int I=i;I<=i+1;I++)
						    for(int J=j;J<=j+1;J++)
							    deltaA_par_local[I-i][J-j]  =   X[s]*deltaA_par_cos[s][I][J]*COS_NPHI_PLUS_OMEGA_T - X[s]*deltaA_par_sin[s][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*deltaA_par_sin[s][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[s]*deltaA_par_cos[s][I][J]*SIN_NPHI_PLUS_OMEGA_T;




				
					    for(int ii=0;ii<2;ii++)
						    for(int jj=0;jj<2;jj++)
						    {
							    for(int d=0;d<3;d++)
							    {
								    grad_deltaA_par_par[d]	+=	grad_deltaA_par_local[d][ii][jj]*weight_square[ii][jj];
								    grad_deltaPhi_par[d]	+=	grad_deltaPhi_local[d][ii][jj]*weight_square[ii][jj];
							    }
							    deltaA_par_par += deltaA_par_local[ii][jj]*weight_square[ii][jj];
						    }


                    }
				    deltaB_par = times(grad_deltaA_par_par,b_par) + deltaA_par_par*curl_b_par;
				    deltaE_par = (b_par*grad_deltaPhi_par)*b_par - grad_deltaPhi_par ;

				    B_s = B_par + deltaB_par + alpha*m*marker[p].v_par/e*curl_b_par;
				    B_ss = B_s*b_par;


				    R = marker[p].X[0];
				    Z = marker[p].X[1];

                    //add by 2013.12.20 to modify W_per = mu*(B+deltaB)
                    double total_B_value_par;
                    total_B_value_par = sqrt((B_par + deltaB_par)*(B_par + deltaB_par));

				    W_par = 0.5*species[0].mass*pow(marker[p].v_par,2.0);
				    //W_per = marker[p].mu*B_value_par;
                    W_per = marker[p].mu*total_B_value_par;

				    double K =  species[0].charge*(marker[p].X[0]*deltaA_par_par)/alpha*B_par[2]/B_value_par + species[0].charge*(psi_par)/alpha + species[0].mass*(marker[p].v_par)*marker[p].X[0]*B_par[2]/B_value_par - n[0]/omega_A[0]*(W_per+W_par+species[0].charge*deltaPhi_par);
				    double p_phi_unperturbed =  species[0].charge*(psi_par)/alpha + species[0].mass*(marker[p].v_par)*marker[p].X[0]*B_par[2]/B_value_par;
	

				    double p_phi1 = species[0].charge*(marker[p].X[0]*deltaA_par_par)/alpha*B_par[2]/B_value_par + species[0].charge*psi_par/alpha;
				    double p_phi2 = species[0].mass*marker[p].v_par*marker[p].X[0]*B_par[2]/B_value_par;
				    double p_phi3 = (W_per+W_par+species[0].charge*deltaPhi_par);
				    double p_phi4 = species[0].charge*(marker[p].X[0]*deltaA_par_par)/alpha*B_par[2]/B_value_par + species[0].charge*(psi_par)/alpha+ species[0].mass*marker[p].v_par*marker[p].X[0]*B_par[2]/B_value_par;
                    double p_phi5 = (W_per+W_par);
				    double p_phi6 = species[0].charge*(psi_par)/alpha+ species[0].mass*marker[p].v_par*marker[p].X[0]*B_par[2]/B_value_par;
				

				    file_orbit << setiosflags(ios::fixed) << setprecision(12) << marker[p].X[0] << " " << marker[p].X[1] << " " << marker[p].X[2] << "  " << marker[p].v_par << " " << p_phi_unperturbed << " " << setprecision(18) << K << "  " << setiosflags(ios::fixed) << setprecision(2) << t << endl;

				    file_p_phi << setiosflags(ios::scientific) << setprecision(15) << p_phi1 << " "  << p_phi2 << " "  << p_phi3 << " "  << p_phi4 << " " << t << " " << p_phi5 << " " << p_phi6 << endl;

				    file_E_par << setiosflags(ios::fixed) << setiosflags(ios::scientific) << setprecision(12) << deltaE_par[0] << "  " << t << endl;

				    file_B_par << setiosflags(ios::fixed) << setiosflags(ios::scientific) << setprecision(12) << deltaB_par[0] << "  " << t << endl;

			    }
         }


/************************************************************************************
*                                                                                   *
*                      Output K at t=0 for every Particle                           *
*                                                                                   *
*  constant motion : K  = ~P_phi + n/omega_A*~E                                     *
*                    here ~P_phi = P_phi + e*R*deltaA_||*B_phi/B                    *
*                         ~E     = E + e*deltaPhi                                   *
*                                                                                   *
************************************************************************************/

        if(TIMESTEP == 1)
        {
		    if(OKP)
            {
                for(int p=0;p<N_p;p++)
                {
                    i = (int)((marker[p].X[0]-R0+a*a_b)/dR);
				    j = (int)((marker[p].X[1]+a*a_b)/dZ);

				    R   = R0-a*a_b+i*dR;
				    Z   =   -a*a_b+j*dZ;

                    weight_line_auxiliary[0][0] = (marker[p].X[0] - R)/dR;
				    weight_line_auxiliary[1][0] = (marker[p].X[1] - Z)/dZ;

				    weight_line_auxiliary[0][1] = 1.0 - weight_line_auxiliary[0][0];
				    weight_line_auxiliary[1][1] = 1.0 - weight_line_auxiliary[1][0];
				

				    for(int ii=0;ii<2;ii++)
					    for(int jj=0;jj<2;jj++)
						    weight_square[ii][jj] = weight_line_auxiliary[0][1-ii]*weight_line_auxiliary[1][1-jj];

				    B_par       = 0.0;
				    curl_b_par  = 0.0;
				    grad_B_par  = 0.0;
				    B_value_par = 0.0;
				
				    for(int ii=0;ii<2;ii++)
					    for(int jj=0;jj<2;jj++)
					    {
						    B_value_par		+=	B_value[i+ii][j+jj]*weight_square[ii][jj];
						    for(int d=0;d<2;d++)
						    {
							    B_par[d]      +=      B[d][i+ii][j+jj]*weight_square[ii][jj];
							    grad_B_par[d] += grad_B[d][i+ii][j+jj]*weight_square[ii][jj];
						    }

						    for(int d=0;d<3;d++)
							    curl_b_par[d] += curl_b[d][i+ii][j+jj]*weight_square[ii][jj];
					    }

				    if(NUMERICAL_EQUILIBRIUM)
                    {
                        double g_eq_par = 0.0;
                        for(int ii=0;ii<2;ii++)
			                for(int jj=0;jj<2;jj++)
				                g_eq_par		+=	g_eq[i+ii][j+jj]*weight_square[ii][jj];

                        B_par[2] = g_eq_par/marker[p].X[0];
                    }
                    else
                        B_par[2] = B0*R0/marker[p].X[0];

				    grad_B_par[2] = 0.0;


				    for(int d=0;d<3;d++)
					    b_par[d] =   B_par[d]/B_value_par;

				    grad_deltaA_par_par = 0.0;
				    grad_deltaPhi_par   = 0.0;
				    deltaA_par_par      = 0.0;
                    for(int s=0;s<NUM_MODE;s++)
				    {
					    COS_NPHI_PLUS_OMEGA_T = cos(n[s]*marker[p].X[2]-omega_A[s]*t);
					    SIN_NPHI_PLUS_OMEGA_T = sin(n[s]*marker[p].X[2]-omega_A[s]*t);

					    for(int d=0;d<3;d++)
						    for(int I=i;I<=i+1;I++)
							    for(int J=j;J<=j+1;J++)
							    {
								    if(d!=2)
									    grad_deltaA_par_local[d][I-i][J-j]  =       X[s]*grad_deltaA_par_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T - X[s]*grad_deltaA_par_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaA_par_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaA_par_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T;
								    else
									    grad_deltaA_par_local[d][I-i][J-j]  =   (-  X[s]*grad_deltaA_par_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T - X[s]*grad_deltaA_par_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T - Y[s]*grad_deltaA_par_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaA_par_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T)*n[s]/marker[p].X[0];
								    if(d!=2)
									    grad_deltaPhi_local[d][I-i][J-j]    =       X[s]*grad_deltaPhi_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T - X[s]*grad_deltaPhi_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaPhi_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaPhi_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T;
								    else
									    grad_deltaPhi_local[d][I-i][J-j]    =   (-  X[s]*grad_deltaPhi_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T - X[s]*grad_deltaPhi_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T - Y[s]*grad_deltaPhi_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaPhi_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T)*n[s]/marker[p].X[0];

							    }
					    for(int I=i;I<=i+1;I++)
						    for(int J=j;J<=j+1;J++)
							    deltaA_par_local[I-i][J-j]  =   X[s]*deltaA_par_cos[s][I][J]*COS_NPHI_PLUS_OMEGA_T - X[s]*deltaA_par_sin[s][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*deltaA_par_sin[s][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[s]*deltaA_par_cos[s][I][J]*SIN_NPHI_PLUS_OMEGA_T;




				
					    for(int ii=0;ii<2;ii++)
						    for(int jj=0;jj<2;jj++)
						    {
							    for(int d=0;d<3;d++)
							    {
								    grad_deltaA_par_par[d]	+=	grad_deltaA_par_local[d][ii][jj]*weight_square[ii][jj];
								    grad_deltaPhi_par[d]	+=	grad_deltaPhi_local[d][ii][jj]*weight_square[ii][jj];
							    }
							    deltaA_par_par += deltaA_par_local[ii][jj]*weight_square[ii][jj];
						    }


                    }
				    deltaB_par = times(grad_deltaA_par_par,b_par) + deltaA_par_par*curl_b_par;
				    deltaE_par = (b_par*grad_deltaPhi_par)*b_par - grad_deltaPhi_par ;

				    B_s = B_par + deltaB_par + alpha*m*marker[p].v_par/e*curl_b_par;
				    B_ss = B_s*b_par;


				    R = marker[p].X[0];
				    Z = marker[p].X[1];

				    W_par = 0.5*species[0].mass*pow(marker[p].v_par,2.0);
				    W_per = marker[p].mu*B_value_par;

				    double K =  species[0].charge*(marker[p].X[0]*deltaA_par_par)/alpha*B_par[2]/B_value_par + species[0].charge*(psi_par)/alpha + species[0].mass*(marker[p].v_par)*marker[p].X[0]*B_par[2]/B_value_par - n[0]/omega_A[0]*(W_per+W_par+species[0].charge*deltaPhi_par);

                    marker[p].K = K;

                }
            }
        }


		//initial J_dot_E
        for(int s=0;s<NUM_MODE;s++)
        {
            myJ_dot_E_X[s] = 0.0;
		    myJ_dot_E_Y[s] = 0.0;
        }


#pragma omp parallel for \
		default(shared) \
		firstprivate(t) \
		private(B_par,b_par,curl_b_par,grad_B_par,B_value_par,psi_par,\
				deltaB_par,deltaE_par,deltaPhi_par,deltaA_par_par,grad_deltaA_par_par,grad_deltaPhi_par,\
				deltaE_X_par,deltaE_Y_par,\
				Va,RHS_V,Vold,\
                Ua,RHS_U,Uold,lambda,\
				Wa,RHS_W,Wold,\
				Xa,RHS_X,Xold,\
				timestep,coeff,told,\
				OUT_OF_BOUDARY,\
				weight_line,weight_square,weight_line_aux,weight_line_auxiliary,\
				i,j,R,Z,phi,B_s,B_ss,W_per,W_par,mu,\
				E,v,v_d,\
				grad_p_phi,dXdt1,dv_pardt1,\
				p_phi,\
				df0dt,dp_phidt,df_0dP_phi,dEdt,df_0dE,dLambdadt,df_0dLambda\
				deltaE_local,deltaE_X_local,deltaE_Y_local,deltaB_local,deltaPhi_local,deltaA_par_local,grad_deltaA_par_local,grad_deltaPhi_local,grad_deltaPhi_X_local,grad_deltaPhi_Y_local,\
                grad_deltaPhi_X_par,grad_deltaPhi_Y_par,\
                dZ_of_step1,Z_reverse_sign,\
                COS_NPHI_PLUS_OMEGA_T,SIN_NPHI_PLUS_OMEGA_T)
		for(int p=0;p<N_p;p++)
        {
			OUT_OF_BOUDARY = false;
			mu = marker[p].mu;
			for(int RK4=1;RK4<=4;RK4++)
			{
				switch(RK4){
				case 1:
					{
						Xold     = marker[p].X;
						Vold     = marker[p].v_par;
						Wold     = marker[p].w;
						Xa       = marker[p].X;
						Va       = marker[p].v_par;
						Wa       = marker[p].w;
                        if(SLOWING_DOWN)
                        {
                            Uold	= marker[p].v_per;
                            Ua      = marker[p].v_per;
                        }
						told     = t;
						timestep = 0.5*dt;
						coeff    = 1.0/6.0;
					}break;
				case 2:
					{
						timestep = 0.5*dt;
						coeff    = 1.0/3.0;
					}break;
				case 3:
					{
						timestep = dt;
						coeff    = 1.0/3.0;
					}break;
				case 4:
					{
						timestep = dt;
						coeff    = 1.0/6.0;
					}break;
				}

				i = (int)((marker[p].X[0]-R0+a*a_b)/dR);
				j = (int)((marker[p].X[1]+a*a_b)/dZ);



				R   = R0-a*a_b+i*dR;
				Z   =   -a*a_b+j*dZ;


 
				weight_line_auxiliary[0][0] = (marker[p].X[0] - R)/dR;
				weight_line_auxiliary[1][0] = (marker[p].X[1] - Z)/dZ;

				weight_line_auxiliary[0][1] = 1.0 - weight_line_auxiliary[0][0];
				weight_line_auxiliary[1][1] = 1.0 - weight_line_auxiliary[1][0];
				

				for(int ii=0;ii<2;ii++)
					for(int jj=0;jj<2;jj++)
						weight_square[ii][jj] = weight_line_auxiliary[0][1-ii]*weight_line_auxiliary[1][1-jj];


				B_par       = 0.0;
				curl_b_par  = 0.0;
				grad_B_par  = 0.0;
				B_value_par = 0.0;
				
				for(int ii=0;ii<2;ii++)
					for(int jj=0;jj<2;jj++)
					{
						B_value_par		+=	B_value[i+ii][j+jj]*weight_square[ii][jj];
						for(int d=0;d<2;d++)
						{
							B_par[d]      +=      B[d][i+ii][j+jj]*weight_square[ii][jj];
							grad_B_par[d] += grad_B[d][i+ii][j+jj]*weight_square[ii][jj];
						}

						for(int d=0;d<3;d++)
							curl_b_par[d] += curl_b[d][i+ii][j+jj]*weight_square[ii][jj];
					}

				if(NUMERICAL_EQUILIBRIUM)
                {
                    double g_eq_par = 0.0;
                    for(int ii=0;ii<2;ii++)
			            for(int jj=0;jj<2;jj++)
				            g_eq_par		+=	g_eq[i+ii][j+jj]*weight_square[ii][jj];

                    B_par[2] = g_eq_par/marker[p].X[0];
                }
                else
                    B_par[2] = B0*R0/marker[p].X[0];

                grad_B_par[2] = 0.0;


				for(int d=0;d<3;d++)
					b_par[d] =   B_par[d]/B_value_par;

				grad_deltaA_par_par = 0.0;
				grad_deltaPhi_par   = 0.0;
				deltaA_par_par      = 0.0;
                for(int s=0;s<NUM_MODE;s++)
				{
					COS_NPHI_PLUS_OMEGA_T = cos(n[s]*marker[p].X[2]-omega_A[s]*t);
					SIN_NPHI_PLUS_OMEGA_T = sin(n[s]*marker[p].X[2]-omega_A[s]*t);

					for(int d=0;d<3;d++)
						for(int I=i;I<=i+1;I++)
							for(int J=j;J<=j+1;J++)
							{
								if(d!=2)
									grad_deltaA_par_local[d][I-i][J-j]  =       X[s]*grad_deltaA_par_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T - X[s]*grad_deltaA_par_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaA_par_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaA_par_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T;
								else
									grad_deltaA_par_local[d][I-i][J-j]  =   (-  X[s]*grad_deltaA_par_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T - X[s]*grad_deltaA_par_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T - Y[s]*grad_deltaA_par_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaA_par_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T)*n[s]/marker[p].X[0];
								if(d!=2)
									grad_deltaPhi_local[d][I-i][J-j]    =       X[s]*grad_deltaPhi_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T - X[s]*grad_deltaPhi_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaPhi_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaPhi_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T;
								else
									grad_deltaPhi_local[d][I-i][J-j]    =   (-  X[s]*grad_deltaPhi_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T - X[s]*grad_deltaPhi_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T - Y[s]*grad_deltaPhi_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaPhi_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T)*n[s]/marker[p].X[0];

							}
					for(int I=i;I<=i+1;I++)
						for(int J=j;J<=j+1;J++)
							deltaA_par_local[I-i][J-j]  =   X[s]*deltaA_par_cos[s][I][J]*COS_NPHI_PLUS_OMEGA_T - X[s]*deltaA_par_sin[s][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*deltaA_par_sin[s][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[s]*deltaA_par_cos[s][I][J]*SIN_NPHI_PLUS_OMEGA_T;



				
					for(int ii=0;ii<2;ii++)
						for(int jj=0;jj<2;jj++)
						{
							for(int d=0;d<3;d++)
							{
								grad_deltaA_par_par[d]	+=	grad_deltaA_par_local[d][ii][jj]*weight_square[ii][jj];
								grad_deltaPhi_par[d]	+=	grad_deltaPhi_local[d][ii][jj]*weight_square[ii][jj];
							}
							deltaA_par_par += deltaA_par_local[ii][jj]*weight_square[ii][jj];
						}


                }
				deltaB_par = times(grad_deltaA_par_par,b_par) + deltaA_par_par*curl_b_par;
				deltaE_par = (b_par*grad_deltaPhi_par)*b_par - grad_deltaPhi_par ;

				B_s = B_par + deltaB_par + alpha*m*marker[p].v_par/e*curl_b_par;
				B_ss = B_s*b_par;



				if(!FULL_F && !RES && !COF)
                {
/*----------------------------------------------------------------------|
|                                                                       |
|               delta-f  weight equation : dw/dt = - 1/g*df_0dt         |
|                                                                       |
-----------------------------------------------------------------------*/

					R = marker[p].X[0];
					Z = marker[p].X[1];

					E	=	0.5*m*pow(marker[p].v_par,2.0)+mu*B_value_par;
					v	=	sqrt(2*E/species[0].mass);

					dXdt1     = 1.0/B_ss*(marker[p].v_par*deltaB_par - times(b_par,deltaE_par));
					dv_pardt1 = e/m/B_ss*(B_s*deltaE_par/alpha   - mu/e*deltaB_par*grad_B_par);

					v_d = 1.0/(e*B_ss)*(alpha*m*pow(marker[p].v_par,2.0)*curl_b_par + mu*times(b_par,grad_B_par));

					if(SHAFRANOV)
					{
						psi_par = 0;
						for(int ii=0;ii<2;ii++)
							for(int jj=0;jj<2;jj++)
								psi_par += psi[i+ii][j+jj]*weight_square[ii][jj];
					}
					else
					{
						double r2 = (R-R0)*(R-R0)+Z*Z;
						psi_par = psi_0*(q_0+1) + 0.5*sqrt(4*pow(psi_0*(q_0+1),2.0)-4*psi_0*B0*(r2-a*a));
					}


                    if(NUMERICAL_EQUILIBRIUM)
                    {
                        particle_vector<double> grad_g_eq_par;
                        grad_g_eq_par = 0.0;
                        for(int ii=0;ii<2;ii++)
						    for(int jj=0;jj<2;jj++)
							    for(int d=0;d<3;d++)
								    grad_g_eq_par[d]	+=	grad_g_eq[d][ii][jj]*weight_square[ii][jj];

                        grad_p_phi[0] =   m*marker[p].v_par/B_value_par*grad_g_eq_par[0] - m*marker[p].v_par*R*B_par[2]/B_value_par/B_value_par*grad_B_par[0] - e*B_par[1]*R/alpha;
					    grad_p_phi[1] =	  m*marker[p].v_par/B_value_par*grad_g_eq_par[1] - m*marker[p].v_par*R*B_par[2]/B_value_par/B_value_par*grad_B_par[1] + e*B_par[0]*R/alpha;
					    grad_p_phi[2] =   0.0;
                    }
                    else
                    {
					    grad_p_phi[0] =   - m*marker[p].v_par*R*B_par[2]/B_value_par/B_value_par*grad_B_par[0] - e*B_par[1]*R/alpha;
					    grad_p_phi[1] =	  - m*marker[p].v_par*R*B_par[2]/B_value_par/B_value_par*grad_B_par[1] + e*B_par[0]*R/alpha;
					    grad_p_phi[2] =   0.0;
                    }



					dp_phidt = dXdt1*grad_p_phi + dv_pardt1*m*R*B_par[2]/B_value_par;
					dEdt     = e*v_d*deltaE_par + e*marker[p].v_par*deltaB_par*deltaE_par/B_ss;

					if(FGP)
						dLambdadt = -marker[p].mu*B0/pow(E,2.0)*dEdt;


/*----------------------------------------------------------------------------------|
|                                                                                   |
|                           dw/dt = -1/g*df_0/dt                                    |
|                                                                                   |			
-----------------------------------------------------------------------------------*/
				
					p_phi = m*marker[p].v_par*marker[p].X[0]*(B_par[2]/B_value_par)+e*psi_par/alpha;

                    E  = 0.5*species[0].mass*pow(marker[p].v_par,2.0) + marker[p].mu*B_value_par;
		            v = sqrt(2*E/species[0].mass);
					double bracket_psi;
					if(AVERAGE_PSI)
					{
                        double mu = marker[p].mu;
                        double Lambda = marker[p].mu*B0/E;
                        double f_0;
						int sgn_v_parallel = marker[p].v_par>0?1:-1;
						if((1-mu*B0/E) > 0)
							bracket_psi = p_phi/species[0].charge - species[0].mass/species[0].charge*sgn_v_parallel*v*R0*sqrt(1-mu*B0/E);
						else
							bracket_psi = p_phi/species[0].charge;

                        //add by 2014.10.22
						//Finite width of Gaussian distribution of Pitch angle  : f(Lambda) ~ exp[(Lambda-Lambda_0)^2/Delta_Lambda^2]


						if(FGP)
                            f_0 = c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(bracket_psi/(c_1*Delta_psi)))*exp(-pow((Lambda-Lambda_0_FGP)/Delta_Lambda,2.0));
                        else
    						f_0 = c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(bracket_psi/(c_1*Delta_psi)));
					
						df_0dP_phi =  -1.0/(e*c_1*Delta_psi)*f_0;
																						 
						if((1-mu*B0/E) > 0)
							df_0dE	=	-3/m*sqrt(2*E/m)/(pow(2*E/m,1.5)+ pow(v_c,3.0))*f_0
										-1/(sqrt(2*PI*E*m)*0.1*v_A*(1+erf((v_0-v)/(0.2*v_A))))*exp(-(pow((v_0-v)/(0.2*v_A),2.0)))*f_0
										+sgn_v_parallel*R0/(species[0].charge*c_1*Delta_psi*sqrt(v*v-2*mu*B0/species[0].mass))*f_0;
						else
							df_0dE	=	-3/m*sqrt(2*E/m)/(pow(2*E/m,1.5)+ pow(v_c,3.0))*f_0
										-1/(sqrt(2*PI*E*m)*0.1*v_A*(1+erf((v_0-v)/(0.2*v_A))))*exp(-(pow((v_0-v)/(0.2*v_A),2.0)))*f_0;

                        if(FGP)
                            df_0dLambda =   -2*(Lambda-Lambda_0_FGP)/pow(Delta_Lambda,2.0)*f_0;
					}
					else
					{
						//add by 2013.9.20
						//Finite width of Gaussian distribution of Pitch angle  : f(Lambda) ~ exp[(Lambda-Lambda_0)^2/Delta_Lambda^2]
						double Lambda = marker[p].mu*B0/E;
						if(FGP)
						{
                            if(DFGP)//double pitch angle distribution
                            {
							    bracket_psi =   p_phi/species[0].charge;
							    double f_0  =   c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(bracket_psi/(c_1*Delta_psi)));

							    df_0dP_phi	=	-1.0/(e*c_1*Delta_psi)*c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(p_phi/(e*c_1*Delta_psi)))*(c1_FGP*exp(-pow((Lambda-Lambda_1_FGP)/Delta_Lambda_1,2.0)) + c2_FGP*exp(-pow((Lambda-Lambda_2_FGP)/Delta_Lambda_2,2.0)));

							    df_0dE		=	(-3/m*sqrt(2*E/m)/(pow(2*E/m,1.5)+ pow(v_c,3.0))*c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(p_phi/(e*c_1*Delta_psi)))
										         -1/(sqrt(2*PI*E*m)*0.1*v_A)*exp(-(pow((v_0-v)/(0.2*v_A),2.0)))*c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*exp(-(p_phi/(e*c_1*Delta_psi))))*(c1_FGP*exp(-pow((Lambda-Lambda_1_FGP)/Delta_Lambda_1,2.0)) + c2_FGP*exp(-pow((Lambda-Lambda_2_FGP)/Delta_Lambda_2,2.0)));

							    df_0dLambda =   (-2*c1_FGP*(Lambda-Lambda_1_FGP)/pow(Delta_Lambda_1,2.0)*exp(-pow((Lambda-Lambda_1_FGP)/Delta_Lambda_1,2.0))
                                                 -2*c2_FGP*(Lambda-Lambda_2_FGP)/pow(Delta_Lambda_2,2.0)*exp(-pow((Lambda-Lambda_2_FGP)/Delta_Lambda_2,2.0)))*f_0;
                            }
                            else
                            {
                                bracket_psi =   p_phi/species[0].charge;
							    double f_0  =   c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(bracket_psi/(c_1*Delta_psi)))*exp(-pow((Lambda-Lambda_0_FGP)/Delta_Lambda,2.0));

							    df_0dP_phi	=	-1.0/(e*c_1*Delta_psi)*c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(p_phi/(e*c_1*Delta_psi)))*exp(-pow((Lambda-Lambda_0_FGP)/Delta_Lambda,2.0));

							    df_0dE		=	(-3/m*sqrt(2*E/m)/(pow(2*E/m,1.5)+ pow(v_c,3.0))*c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(p_phi/(e*c_1*Delta_psi)))
										         -1/(sqrt(2*PI*E*m)*0.1*v_A)*exp(-(pow((v_0-v)/(0.2*v_A),2.0)))*c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*exp(-(p_phi/(e*c_1*Delta_psi))))*exp(-pow((Lambda-Lambda_0_FGP)/Delta_Lambda,2.0));

							    df_0dLambda =   -2*(Lambda-Lambda_0_FGP)/pow(Delta_Lambda,2.0)*f_0;
                            }
						}
						else
						{
							df_0dP_phi	=	-1.0/(e*c_1*Delta_psi)*c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(p_phi/(e*c_1*Delta_psi)));

							df_0dE		=	-3/m*sqrt(2*E/m)/(pow(2*E/m,1.5)+ pow(v_c,3.0))*c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(p_phi/(e*c_1*Delta_psi)))
											-1/(sqrt(2*PI*E*m)*0.1*v_A)*exp(-(pow((v_0-v)/(0.2*v_A),2.0)))*c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*exp(-(p_phi/(e*c_1*Delta_psi)));
						}
					}

					//add by 2013.9.20
					//Finite width of Gaussian distribution of Pitch angle  : f(Lambda) ~ exp[(Lambda-Lambda_0)^2/Delta_Lambda^2]
					if(FGP)
						df0dt	=	dp_phidt*df_0dP_phi + dEdt*df_0dE + dLambdadt*df_0dLambda;
					else
						df0dt	=	dp_phidt*df_0dP_phi + dEdt*df_0dE;



					if(abs(df_0dE) > 1.0/c_f*1e10)
					{
						cout << "df_0dE : " << df_0dE << " v : " << v << " (1e-20 + 1+erf((v_0-v)/(0.2*v_A))) : " << (1e-20 + 1+erf((v_0-v)/(0.2*v_A))) << " exp(-(pow((v_0-v)/(0.2*v_A),2.0))) : " << exp(-(pow((v_0-v)/(0.2*v_A),2.0))) << " myid : " << myid << endl;  
						exit(0);
					}

                }


                if(FULL_F || RES || COF)
				{
					if(SHAFRANOV)
					{
						psi_par = 0;
						for(int ii=0;ii<2;ii++)
							for(int jj=0;jj<2;jj++)
								psi_par += psi[i+ii][j+jj]*weight_square[ii][jj];
					}
					else
					{
						double r2 = (R-R0)*(R-R0)+Z*Z;
						psi_par = psi_0*(q_0+1) + 0.5*sqrt(4*pow(psi_0*(q_0+1),2.0)-4*psi_0*B0*(r2-a*a));
					}
                }


				RHS_X    = (1.0/B_ss*(marker[p].v_par*B_s + times(b_par,alpha*mu/e*grad_B_par - deltaE_par)));
				RHS_V    = (1.0/m/B_ss*B_s*(e*deltaE_par/alpha-mu*grad_B_par));
				RHS_W    = (-df0dt)/marker[p].g;


                if(SLOWING_DOWN)
	            {
		            //record old particle information
		            v		=	sqrt(pow(marker[p].v_par,2.0)+pow(marker[p].v_per,2.0));
		            lambda	=	marker[p].v_par/v;


		            //slowing down particle process
		            RHS_V   +=  lambda*(-nu*(v+pow(v_c,3.0)/pow(v,2.0)));
                    RHS_U   =   sqrt(1-pow(lambda,2.0))*(-nu*(v+pow(v_c,3.0)/pow(v,2.0)));

		            //update particle information
		            marker[p].v_per = Uold + timestep*RHS_U;
                    Ua += coeff*dt*RHS_U;
		
		            //mu should change
		            mu	=	0.5*species[0].mass*pow(marker[p].v_per,2.0)/B_value_par;

	            }

				RHS_X[2] = RHS_X[2]/marker[p].X[0];

				marker[p].X     = Xold + timestep*RHS_X;
				marker[p].v_par = Vold + timestep*RHS_V;
				marker[p].w     = Wold + timestep*RHS_W;
				Xa += coeff*dt*RHS_X;
				Va += coeff*dt*RHS_V;
				Wa += coeff*dt*RHS_W;
				t = told + timestep;


				//recycle particles
				if(RECYCLE_METHOD == 1)
				{
                    if(NUMERICAL_EQUILIBRIUM)
                    {
                        if(abs(psi_par) < 1e-3)
                        {
                            OUT_OF_BOUDARY = true;
						    break;
                        }
                    }
                    else
                    {
                        if(sqrt(pow(marker[p].X[0]-R0,2.0)+pow(marker[p].X[1],2.0)) >= a)
					    {
						    OUT_OF_BOUDARY = true;
						    break;
					    }
                    }
				}
				else if(RECYCLE_METHOD == 2)
				{
                    if(NUMERICAL_EQUILIBRIUM)
                    {
                        if(abs(psi_par) < 1e-8)
                        {
						    marker[p].id = -1;
						    break;
                        }
                    }
                    else
                    {
					    if(sqrt(pow(marker[p].X[0]-R0,2.0)+pow(marker[p].X[1],2.0)) >= a)
					    {
						    marker[p].id = -1;
						    break;
					    }
                    }
				}
            }


			if(OUT_OF_BOUDARY)
			{
				marker[p].v_par = Vold;
				marker[p].X     = Xold;
				marker[p].w     = Wold;
				marker[p].X[1]  =  - marker[p].X[1];
                if(SLOWING_DOWN)
                    marker[p].v_per = Uold;
                //add by 2012-3-14
                if(marker[p].id == 0)
                    LOSS = true ;
			}
			else
			{
				marker[p].X     = Xa;
				marker[p].v_par = Va;
				marker[p].w     = Wa;
                if(SLOWING_DOWN)
                    marker[p].v_per = Ua;
                //add by 2012-3-16
                if(TIMESTEP == 1)
                    dZ_of_step1 = marker[p].X[1] - Xold[1];
                if(marker[p].X[1]*dZ_of_step1 < 0)
                    Z_reverse_sign = true;
			}


			//full f
			if(FULL_F)
                if(!RES)
                    marker[p].w	= marker[p].f_over_g;

            if(RES) // avoid value of marker[p].w is changed, because marker[p].w=Z in RES module for plotting KAM surfaces !!update marker_diagnosis info
                marker[p].w	= Xold[1];

			
			i = (int)((marker[p].X[0]-R0+a*a_b)/dR);
			j = (int)((marker[p].X[1]+a*a_b)/dZ);


			R   = R0-a*a_b+i*dR;
			Z   =   -a*a_b+j*dZ;



			weight_line_auxiliary[0][0] = (marker[p].X[0] - R)/dR;
			weight_line_auxiliary[1][0] = (marker[p].X[1] - Z)/dZ;

			weight_line_auxiliary[0][1] = 1.0 - weight_line_auxiliary[0][0];
			weight_line_auxiliary[1][1] = 1.0 - weight_line_auxiliary[1][0];
				

			for(int ii=0;ii<2;ii++)
				for(int jj=0;jj<2;jj++)
					weight_square[ii][jj] = weight_line_auxiliary[0][1-ii]*weight_line_auxiliary[1][1-jj];

			B_par       = 0.0;
			curl_b_par  = 0.0;
			grad_B_par  = 0.0;
			B_value_par = 0.0;
				
			for(int ii=0;ii<2;ii++)
				for(int jj=0;jj<2;jj++)
				{
					B_value_par		+=	B_value[i+ii][j+jj]*weight_square[ii][jj];
					for(int d=0;d<2;d++)
					{
						B_par[d]      +=      B[d][i+ii][j+jj]*weight_square[ii][jj];
						grad_B_par[d] += grad_B[d][i+ii][j+jj]*weight_square[ii][jj];
					}

					for(int d=0;d<3;d++)
						curl_b_par[d] += curl_b[d][i+ii][j+jj]*weight_square[ii][jj];
				}

			if(NUMERICAL_EQUILIBRIUM)
            {
                double g_eq_par = 0.0;
                for(int ii=0;ii<2;ii++)
			        for(int jj=0;jj<2;jj++)
				        g_eq_par		+=	g_eq[i+ii][j+jj]*weight_square[ii][jj];

                B_par[2] = g_eq_par/marker[p].X[0];
            }
            else
                B_par[2] = B0*R0/marker[p].X[0];

			grad_B_par[2] = 0.0;

			for(int d=0;d<3;d++)
				b_par[d] =   B_par[d]/B_value_par;


            grad_deltaA_par_par = 0.0;
			deltaA_par_par      = 0.0;
			deltaPhi_par        = 0.0;

            for(int s=0;s<NUM_MODE;s++)
			{
                //add by 2014.05.20
                if(IMC)
                {
                    dgrad_deltaA_par_Xdt_par = 0.0;
                    dgrad_deltaA_par_Ydt_par = 0.0;
                    ddeltaA_par_Xdt_par      = 0.0;
                    ddeltaA_par_Ydt_par      = 0.0;
                }

                grad_deltaPhi_X_par = 0.0;
			    grad_deltaPhi_Y_par = 0.0;

				COS_NPHI_PLUS_OMEGA_T = cos(n[s]*marker[p].X[2]-omega_A[s]*t);
				SIN_NPHI_PLUS_OMEGA_T = sin(n[s]*marker[p].X[2]-omega_A[s]*t);

				for(int d=0;d<3;d++)
					for(int I=i;I<=i+1;I++)
						for(int J=j;J<=j+1;J++)
						{

							if(d!=2)
							{
								grad_deltaPhi_X_local[d][I-i][J-j] =    grad_deltaPhi_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T - grad_deltaPhi_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T;
								grad_deltaPhi_Y_local[d][I-i][J-j] =    grad_deltaPhi_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T + grad_deltaPhi_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T;
							}
							else
							{
								grad_deltaPhi_X_local[d][I-i][J-j] =  (-grad_deltaPhi_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T - grad_deltaPhi_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T)*n[s]/marker[p].X[0];
								grad_deltaPhi_Y_local[d][I-i][J-j] =  (-grad_deltaPhi_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T + grad_deltaPhi_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T)*n[s]/marker[p].X[0];
							}
							if(d!=2)
								grad_deltaA_par_local[d][I-i][J-j] =   X[s]*grad_deltaA_par_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T - X[s]*grad_deltaA_par_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaA_par_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaA_par_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T;
							else
								grad_deltaA_par_local[d][I-i][J-j] = (-X[s]*grad_deltaA_par_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T - X[s]*grad_deltaA_par_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T - Y[s]*grad_deltaA_par_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaA_par_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T)*n[s]/marker[p].X[0];
						}

				for(int I=i;I<=i+1;I++)
					for(int J=j;J<=j+1;J++)
					{
						deltaA_par_local[I-i][J-j]  =	X[s]*deltaA_par_cos[s][I][J]*COS_NPHI_PLUS_OMEGA_T - X[s]*deltaA_par_sin[s][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*deltaA_par_sin[s][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[s]*deltaA_par_cos[s][I][J]*SIN_NPHI_PLUS_OMEGA_T;
						deltaPhi_local[I-i][J-j]    =	X[s]*deltaPhi_cos[s][I][J]*COS_NPHI_PLUS_OMEGA_T - X[s]*deltaPhi_sin[s][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*deltaPhi_sin[s][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[s]*deltaPhi_cos[s][I][J]*SIN_NPHI_PLUS_OMEGA_T;
					}
		
				
				for(int ii=0;ii<2;ii++)
					for(int jj=0;jj<2;jj++)
					{
						for(int d=0;d<3;d++)
						{
							grad_deltaA_par_par[d] += grad_deltaA_par_local[d][ii][jj]*weight_square[ii][jj];
							grad_deltaPhi_X_par[d] += grad_deltaPhi_X_local[d][ii][jj]*weight_square[ii][jj];
							grad_deltaPhi_Y_par[d] += grad_deltaPhi_Y_local[d][ii][jj]*weight_square[ii][jj];
							
						}
						deltaA_par_par	+= deltaA_par_local[ii][jj]*weight_square[ii][jj];
						deltaPhi_par	+= deltaPhi_local[ii][jj]*weight_square[ii][jj];
					}

                deltaE_X_par[s] = (b_par*grad_deltaPhi_X_par)*b_par - grad_deltaPhi_X_par;
				deltaE_Y_par[s] = (b_par*grad_deltaPhi_Y_par)*b_par - grad_deltaPhi_Y_par;

                //add by 2014.05.20
                //EP current J_h Including Magnetization Current
                if(IMC)
                {
                    for(int d=0;d<3;d++)
					    for(int I=i;I<=i+1;I++)
						    for(int J=j;J<=j+1;J++)
						    {
	    					    if(d!=2)
                                {
								    dgrad_deltaA_par_Xdt_local[d][I-i][J-j] =   -omega_A[s]*(-grad_deltaA_par_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T - grad_deltaA_par_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T);
                                    dgrad_deltaA_par_Ydt_local[d][I-i][J-j] =   -omega_A[s]*(-grad_deltaA_par_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T + grad_deltaA_par_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T);
                                }
							    else
                                {
								    dgrad_deltaA_par_Xdt_local[d][I-i][J-j] = -omega_A[s]*(-grad_deltaA_par_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T + grad_deltaA_par_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T)*n[s]/marker[p].X[0];
                                    dgrad_deltaA_par_Ydt_local[d][I-i][J-j] = -omega_A[s]*(-grad_deltaA_par_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T - grad_deltaA_par_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T)*n[s]/marker[p].X[0];
                                }
						    }

                    for(int I=i;I<=i+1;I++)
					    for(int J=j;J<=j+1;J++)
					    {
						    ddeltaA_par_Xdt_local[I-i][J-j]  = -omega_A[s]*(-deltaA_par_cos[s][I][J]*SIN_NPHI_PLUS_OMEGA_T - deltaA_par_sin[s][I][J]*COS_NPHI_PLUS_OMEGA_T);
                            ddeltaA_par_Ydt_local[I-i][J-j]  = -omega_A[s]*(-deltaA_par_sin[s][I][J]*SIN_NPHI_PLUS_OMEGA_T + deltaA_par_cos[s][I][J]*COS_NPHI_PLUS_OMEGA_T);
					    }

                    for(int ii=0;ii<2;ii++)
					    for(int jj=0;jj<2;jj++)
					    {
						    for(int d=0;d<3;d++)
						    {
							    dgrad_deltaA_par_Xdt_par[d] += dgrad_deltaA_par_Xdt_local[d][ii][jj]*weight_square[ii][jj];
                                dgrad_deltaA_par_Ydt_par[d] += dgrad_deltaA_par_Ydt_local[d][ii][jj]*weight_square[ii][jj];
							
						    }
						    ddeltaA_par_Xdt_par += ddeltaA_par_Xdt_local[ii][jj]*weight_square[ii][jj];
                            ddeltaA_par_Ydt_par += ddeltaA_par_Ydt_local[ii][jj]*weight_square[ii][jj];
					    }

                        ddeltaB_Xdt_par[s] = times(dgrad_deltaA_par_Xdt_par,b_par) + ddeltaA_par_Xdt_par*curl_b_par;
                        ddeltaB_Ydt_par[s] = times(dgrad_deltaA_par_Ydt_par,b_par) + ddeltaA_par_Ydt_par*curl_b_par;

                }
				
            }
            deltaB_par = times(grad_deltaA_par_par,b_par) + deltaA_par_par*curl_b_par;



            if(SLOWING_DOWN)
                marker[p].mu	=	0.5*species[0].mass*pow(marker[p].v_per,2.0)/B_value_par;

			//update v_per
			marker[p].v_per = sqrt(2*marker[p].mu*B_value_par/m);

            if(SCATTERING)
            {
            	//if add collison effect,mu should be changed!!
	            double v,lambda_old,lambda_new;
	            double nu_d;
	            const double c=0.17*v_0;
	            int pm;
	            
                v				=	sqrt(pow(marker[p].v_par,2.0)+pow(marker[p].v_per,2.0));
		        nu_d			=	nu*pow(v_c,3.0)/2.0/(pow(c,3.0)+pow(v,3.0));

		        lambda_old		=	marker[p].v_par/v;
		        pm				=	(rand()*1.0/RAND_MAX)>0.5?1:-1; 
		        lambda_new		=	lambda_old*(1-2*nu_d*dt) + pm*sqrt((1-lambda_old*lambda_old)*2*nu_d*dt);


		        marker[p].v_par	=	v*lambda_new;
		        marker[p].v_per	=	sqrt(v*v - pow(marker[p].v_par,2.0));

		        //mu should change
		        marker[p].mu	=	0.5*species[0].mass*pow(marker[p].v_per,2.0)/B_value_par;
            }



			B_s = B_par + deltaB_par + alpha*m*marker[p].v_par/e*curl_b_par;
			B_ss = B_s*b_par;

			v_d = 1.0/(e*B_ss)*(alpha*m*pow(marker[p].v_par,2.0)*curl_b_par + marker[p].mu*times(b_par,grad_B_par));



			//update p_phi
			marker[p].P_phi =  (m*marker[p].v_par*marker[p].X[0]*(B_par[2]/B_value_par)+e*psi_par/alpha);
			
			E = 0.5*m*pow(marker[p].v_par,2.0) + marker[p].mu*B_value_par;

			//add by 2011-2.24 
/*************************************************************************************************************************************
*                                                                                                                                    *
*                                   eliminate adiabatic term using                                                                   *
*            0.retain adiabatic term : f_a = 0.0                                                                                     *
*            1.candy's method        : f_a = e*R*deltaA_||*B_phi/B df_0dP_phi + e*deltaPhi*df_0dE - mu*deltaB*B0/E*df_0dLambda       *
*            2.fu's method           : f_a = -xi \cdot \grad f = e*R*deltaA_||*B_phi/B df_0dP_phi                                    *
*                                                                                                                                    *
*************************************************************************************************************************************/
			//adiabatic part of delta-f
			double deltaf_ad = 0.0;
			if(ADIABATIC != 0)
			{
				double psi_par;
				if(SHAFRANOV)
				{
					psi_par = 0;
					for(int ii=0;ii<2;ii++)
						for(int jj=0;jj<2;jj++)
							psi_par += psi[i+ii][j+jj]*weight_square[ii][jj];
				}
				else
				{
					R = marker[p].X[0];
					Z = marker[p].X[1];
					double r2 = (R-R0)*(R-R0)+Z*Z;
					psi_par = psi_0*(q_0+1) + 0.5*sqrt(4*pow(psi_0*(q_0+1),2.0)-4*psi_0*B0*(r2-a*a));
				}

				v = sqrt(2*E/species[0].mass);

				//p_phi =  (species[0].mass*marker[p].v_par*marker[p].X[0]*(B_par[2]/B_value_par)+species[0].charge*psi_par/alpha);
				p_phi = marker[p].P_phi;

                E  = 0.5*species[0].mass*pow(marker[p].v_par,2.0) + marker[p].mu*B_value_par;
		        v = sqrt(2*E/species[0].mass);

				double bracket_psi;	
				if(AVERAGE_PSI)
				{
					double mu = marker[p].mu;
                    double f_0;
                    double Lambda = marker[p].mu*B0/E;
					int sgn_v_parallel = marker[p].v_par>0?1:-1;
					if((1-mu*B0/E) > 0)
						bracket_psi = p_phi/species[0].charge - species[0].mass/species[0].charge*sgn_v_parallel*v*R0*sqrt(1-mu*B0/E);
					else
						bracket_psi = p_phi/species[0].charge;
					
                    //add by 2014.10.22
					//Finite width of Gaussian distribution of Pitch angle  : f(Lambda) ~ exp[(Lambda-Lambda_0)^2/Delta_Lambda^2]
					if(FGP)
                        f_0 = c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(bracket_psi/(c_1*Delta_psi)))*exp(-pow((Lambda-Lambda_0_FGP)/Delta_Lambda,2.0));
                    else
    					f_0 = c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(bracket_psi/(c_1*Delta_psi)));
					
					df_0dP_phi =  -1.0/(e*c_1*Delta_psi)*f_0;
																						 
					if((1-mu*B0/E) > 0)
						df_0dE	=	-3/m*sqrt(2*E/m)/(pow(2*E/m,1.5)+ pow(v_c,3.0))*f_0
									-1/(sqrt(2*PI*E*m)*0.1*v_A*(1+erf((v_0-v)/(0.2*v_A))))*exp(-(pow((v_0-v)/(0.2*v_A),2.0)))*f_0
									+sgn_v_parallel*R0/(species[0].charge*c_1*Delta_psi*sqrt(v*v-2*mu*B0/species[0].mass))*f_0;
					else
						df_0dE	=	-3/m*sqrt(2*E/m)/(pow(2*E/m,1.5)+ pow(v_c,3.0))*f_0
									-1/(sqrt(2*PI*E*m)*0.1*v_A*(1+erf((v_0-v)/(0.2*v_A))))*exp(-(pow((v_0-v)/(0.2*v_A),2.0)))*f_0;

                    if(FGP)
                        df_0dLambda =   -2*(Lambda-Lambda_0_FGP)/pow(Delta_Lambda,2.0)*f_0;
						
					df_0dP_phi =  -1.0/(e*c_1*Delta_psi)*f_0;
				}
				else
				{
                    //add by 2013.9.20
			        //Finite width of Gaussian distribution of Pitch angle  : f(Lambda) ~ exp[(Lambda-Lambda_0)^2/Delta_Lambda^2]
			        double Lambda = marker[p].mu*B0/E;

                    if(FGP)
                    {
                        if(DFGP)//double pitch angle distribution
                        {
							    bracket_psi =   p_phi/species[0].charge;
							    double f_0  =   c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(bracket_psi/(c_1*(psi_a - psi_0))));

							    df_0dP_phi	=	-1.0/(e*c_1*Delta_psi)*c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(p_phi/(e*c_1*Delta_psi)))*(c1_FGP*exp(-pow((Lambda-Lambda_1_FGP)/Delta_Lambda_1,2.0)) + c2_FGP*exp(-pow((Lambda-Lambda_2_FGP)/Delta_Lambda_2,2.0)));

							    df_0dE		=	(-3/m*sqrt(2*E/m)/(pow(2*E/m,1.5)+ pow(v_c,3.0))*c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(p_phi/(e*c_1*Delta_psi)))
										         -1/(sqrt(2*PI*E*m)*0.1*v_A)*exp(-(pow((v_0-v)/(0.2*v_A),2.0)))*c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*exp(-(p_phi/(e*c_1*Delta_psi))))*(c1_FGP*exp(-pow((Lambda-Lambda_1_FGP)/Delta_Lambda_1,2.0)) + c2_FGP*exp(-pow((Lambda-Lambda_2_FGP)/Delta_Lambda_2,2.0)));

							    df_0dLambda =   (-2*c1_FGP*(Lambda-Lambda_1_FGP)/pow(Delta_Lambda_1,2.0)*exp(-pow((Lambda-Lambda_1_FGP)/Delta_Lambda_1,2.0))
                                                 -2*c2_FGP*(Lambda-Lambda_2_FGP)/pow(Delta_Lambda_2,2.0)*exp(-pow((Lambda-Lambda_2_FGP)/Delta_Lambda_2,2.0)))*f_0;
                        }
                        else
                        {
						    double f_0  =   c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(bracket_psi/(c_1*(psi_a - psi_0))))*exp(-pow((Lambda-Lambda_0_FGP)/Delta_Lambda,2.0));

						    df_0dP_phi	=	-1.0/(e*c_1*Delta_psi)*c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(p_phi/(e*c_1*Delta_psi)))*exp(-pow((Lambda-Lambda_0_FGP)/Delta_Lambda,2.0));

						    df_0dE		=	(-3/m*sqrt(2*E/m)/(pow(2*E/m,1.5)+ pow(v_c,3.0))*c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(p_phi/(e*c_1*Delta_psi)))
										        -1/(sqrt(2*PI*E*m)*0.1*v_A)*exp(-(pow((v_0-v)/(0.2*v_A),2.0)))*c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*exp(-(p_phi/(e*c_1*Delta_psi))))*exp(-pow((Lambda-Lambda_0_FGP)/Delta_Lambda,2.0));

						    df_0dLambda =   -2*(Lambda-Lambda_0_FGP)/pow(Delta_Lambda,2.0)*f_0;
                        }

                    }
                    else
                    {
                        df_0dP_phi	=	-1.0/(e*c_1*Delta_psi)*c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(p_phi/(e*c_1*Delta_psi)));
						df_0dE		=	-3/m*sqrt(2*E/m)/(pow(2*E/m,1.5)+ pow(v_c,3.0))*c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*(1+erf((v_0-v)/(0.2*v_A)))*exp(-(p_phi/(e*c_1*Delta_psi)))
									-1/(sqrt(2*PI*E*m)*0.1*v_A)*exp(-(pow((v_0-v)/(0.2*v_A),2.0)))*c_f*1.0/(pow(v,3)+v_c*v_c*v_c)*exp(-(p_phi/(e*c_1*Delta_psi)));

                    }

				}
                if(ADIABATIC == 1)
                    //add by 2013.9.20
                    if(FGP)
                    {
                        double deltaB_value_par;
                        deltaB_value_par = sqrt(deltaB_par*deltaB_par);
                        deltaf_ad = (e*marker[p].X[0]*deltaA_par_par*B_par[2]/B_value_par/alpha*df_0dP_phi + e*deltaPhi_par*df_0dE -  marker[p].mu*deltaB_value_par*B0/B_value_par/E*df_0dLambda)/marker[p].g;
                    }
                    else
                        deltaf_ad = (e*marker[p].X[0]*deltaA_par_par*B_par[2]/B_value_par/alpha*df_0dP_phi + e*deltaPhi_par*df_0dE)/marker[p].g;
                 else if(ADIABATIC == 2)
                    deltaf_ad = (e*marker[p].X[0]*deltaA_par_par*B_par[2]/B_value_par/alpha*df_0dP_phi)/marker[p].g;
					
            }
			
			t = told;


			if(marker[p].id != -1)
			{
				#pragma omp critical
				{
                    for(int s=0;s<NUM_MODE;s++)
                    {
					    myJ_dot_E_X[s] += e*(v_d+marker[p].v_par*deltaB_par/B_ss)*deltaE_X_par[s]*(marker[p].w - deltaf_ad);
					    myJ_dot_E_Y[s] += e*(v_d+marker[p].v_par*deltaB_par/B_ss)*deltaE_Y_par[s]*(marker[p].w - deltaf_ad);


                        // add by 2014.05.20
                        //EP current J_h Including Magnetization Current
                        if(IMC)
                        {
                            myJ_dot_E_X[s] += marker[p].mu*b_par*ddeltaB_Xdt_par[s]*(marker[p].w - deltaf_ad);
					        myJ_dot_E_Y[s] += marker[p].mu*b_par*ddeltaB_Ydt_par[s]*(marker[p].w - deltaf_ad);
                        }
                    }
				}
			}

//boundary condition
//for toroidal frequency phi should not be used period boundary
            if(!COF)
            {
                while(marker[p].X[2] > 2*PI)
				    marker[p].X[2] -= 2*PI;
			    while(marker[p].X[2] < 0)
				    marker[p].X[2] += 2*PI;
            }


    }

//add by 2011-6-20
//collision processing
//#1 : pitch angle scatering
//	if(SCATTERING)
//		scattering();
//#2 : slowing down particle
//	if(SLOWING_DOWN)
//		slowing_down();
//add by 2011-7-30
//add source and sink
	    if(SOURCE_AND_SINK)
		    source_and_sink();

/************************************************************************************************
*                                                                                               *
*                           diagnostics for P_phi-E phase space structrue                       *
*                                                                                               *
************************************************************************************************/
        if(DIAGNOSTICS && (int(t/dt))%(int(1.0/NUM_TIME_DIAGNOSTICS*(total_time/dt))) == 0)
        {
            if(ISOTROPY)
            {
                //add by 2013-11-21
                //output deltaf in P_phi and Lambda phase space
                if(dLP)
                    diagnostics_output_deltaf_vs_P_phi_and_Lambda(((int(t/dt))/(int(1.0/NUM_TIME_DIAGNOSTICS*(total_time/dt)))),E_0_dLP,windows_width_dLP);
                
                //add by 2013-11-21
                //output deltaf in P_phi and Lambda phase space
                if(dLE)
                    diagnostics_output_deltaf_vs_E_and_Lambda(((int(t/dt))/(int(1.0/NUM_TIME_DIAGNOSTICS*(total_time/dt)))),P_phi_0_dLE,windows_width_dLE);

                //output deltaf in Energy and P_phi phase space
                if(dEP)
                {
                    if(FGP && !DFGP)
                    {
						//add by 2014-07-29
						// fixed mu instead of Lambda to plot deltaf(P_phi,E)
						if (mIL)
						{
							diagnostics_movie_use_MPI_Gather_fixed_mu(((int(t / dt)) / (int(1.0 / NUM_TIME_DIAGNOSTICS*(total_time / dt)))), mu_0_1);
							diagnostics_movie_use_MPI_Gather_fixed_mu(((int(t / dt)) / (int(1.0 / NUM_TIME_DIAGNOSTICS*(total_time / dt)))), mu_0_2);
							diagnostics_movie_use_MPI_Gather_fixed_mu(((int(t / dt)) / (int(1.0 / NUM_TIME_DIAGNOSTICS*(total_time / dt)))), mu_0_3);
							diagnostics_movie_use_MPI_Gather_fixed_mu(((int(t / dt)) / (int(1.0 / NUM_TIME_DIAGNOSTICS*(total_time / dt)))), mu_0_4);
						}
						else
							diagnostics_movie_use_MPI_Gather(((int(t/dt))/(int(1.0/NUM_TIME_DIAGNOSTICS*(total_time/dt)))),Lambda_0_FGP);
                    }
                    else
                    {
                        for(int lambda_index=0;lambda_index<5;lambda_index++)
                        {
                            double Lambda_0_ISOTROPY;
                            switch(lambda_index){
                                case 0:Lambda_0_ISOTROPY = 0.01;break;
                                case 1:Lambda_0_ISOTROPY = 0.25;break;
                                case 2:Lambda_0_ISOTROPY = 0.50;break;
                                case 3:Lambda_0_ISOTROPY = 0.75;break;
                                case 4:Lambda_0_ISOTROPY = 0.99;break;
                            }
		                    diagnostics_movie_use_MPI_Gather(((int(t/dt))/(int(1.0/NUM_TIME_DIAGNOSTICS*(total_time/dt)))),Lambda_0_ISOTROPY);
                        }
                    }
                }
            }
            else
            {
                if(dEP)
                    diagnostics_movie_use_MPI_Gather(((int(t/dt))/(int(1.0/NUM_TIME_DIAGNOSTICS*(total_time/dt)))),0.0);
            }
        }


/************************************************************************************************
*                                                                                               *
*                           output Time Evolution of KAM surfaces                               *
*                                                                                               *
************************************************************************************************/
        if(TEK && (int(t/dt))%(int(1.0/NUM_TIME_DIAGNOSTICS*(total_time/dt))) == 0)
            output_Time_Evolution_of_KAM(((int(t/dt))/(int(1.0/NUM_TIME_DIAGNOSTICS*(total_time/dt)))));


/************************************************************************************************
*                                                                                               *
*                  Evolution of Distribution function as function of Time                       *
*                                                                                               *
************************************************************************************************/
	    if(EDT)
        {
		    if(ISOTROPY)
            {
                diagnostics_output_f_vs_P_phi(0.01,E_0,windows_width);
                diagnostics_output_f_vs_P_phi1(0.25,E_0,windows_width);
                //add by 2013.9.23
			    //Finite width of Gaussian distribution of Pitch angle  : f(Lambda) ~ exp[(Lambda-Lambda_0)^2/Delta_Lambda^2]
                if(DPE)
                {
                    if(FGP && !DFGP)
                        diagnostics_output_f_vs_P_phi_and_E(Lambda_0,windows_width);
                    else
                        diagnostics_output_f_vs_P_phi_and_E(Lambda_0,windows_width);
                }
            }
            else
            {
                diagnostics_output_f_vs_P_phi(Lambda_0,E_0,windows_width);
                if(DPE)
                    diagnostics_output_f_vs_P_phi_and_E(Lambda_0,windows_width);
            }

//            diagnostics_output_f_vs_P_phi(Lambda_0,E_0,5e-2);
//            diagnostics_output_f_vs_P_phi(Lambda_0,E_0,2e-2);
//            diagnostics_output_f_vs_P_phi(Lambda_0,E_0,1e-2);
//            diagnostics_output_f_vs_P_phi(Lambda_0,E_0,5e-3);
        }



	    if(RECYCLE_METHOD == 2)
		    kick_particle();
	    if(AMP)
		    amplitude();

//add by 2013-10-27
        if(TEK && !AMP)
        {
            if(myid == 0)
            {
                ofstream amp_verified("amplitude_verified.dat",ios::app);
                amp_verified << setiosflags(ios::fixed) << setiosflags(ios::scientific) << setprecision(12) << setw(18) << X[0] << " " << setw(18) << Y[0] << resetiosflags(ios::scientific) << resetiosflags(ios::showpoint) << resetiosflags(ios::fixed) << " " << setw(2) << T_INITIA << " " << setw(2) << T_RECORD << " " << setw(9) << t;
                amp_verified << endl;
            }
        }

//add by 2013-10-29
/*************************************************************************************
*             output perturbed field(deltaB_R,deltaE_R) v.s. time                    *
*              at position (R_field,Z_field,phi_field)  (2013.10.29)                 *
*************************************************************************************/
        if(myid == 0)
            if(!COF)
                output_perturbed_field_vs_time(R_field,Z_field,phi_field);

        t = t + dt;

    }
	
}






void amplitude() 
{
	ofstream amp_file("amplitude.dat",ios::app);


    for(int s=0;s<NUM_MODE;s++)
    {
	    MPI::COMM_WORLD.Reduce(&myJ_dot_E_X[s], &J_dot_E_X[s], 1, MPI::DOUBLE, MPI_SUM, 0);
	    MPI::COMM_WORLD.Reduce(&myJ_dot_E_Y[s], &J_dot_E_Y[s], 1, MPI::DOUBLE, MPI_SUM, 0);
    

//	amplitude evlution equation

	    double RHS_X,RHS_Y;
	
	    if(myid == 0)
        {
		    //volume average
		    J_dot_E_X[s] /= (2*PI*PI*a*a*R0);
		    J_dot_E_Y[s] /= (2*PI*PI*a*a*R0);

		    RHS_X = (-J_dot_E_X[s]/(2*W[s]) - gamma_d[s]*X[s])*dt;
		    RHS_Y = (-J_dot_E_Y[s]/(2*W[s]) - gamma_d[s]*Y[s])*dt;

    //add by 2011-4-26
		    if(AMP_MOD_TERM)
		    {
			    RHS_X = RHS_X + (-0.5*omega_A[s]*species[0].charge*Y[s]*mod_X_RHS[s]/(2*W[s])*dt);
			    RHS_Y = RHS_Y + (-0.5*omega_A[s]*species[0].charge*X[s]*mod_Y_RHS[s]/(2*W[s])*dt);
		    }

		

		    if(AVG)
		    {
			    if(T_INITIA < T_AVG_NUM)
			    {
				    X_AVG += RHS_X;
				    Y_AVG += RHS_Y;
				    X_RECORD[T_INITIA] = RHS_X;
				    Y_RECORD[T_INITIA] = RHS_Y;
				    T_INITIA++;
			    }
			    else
			    {	
				    X_AVG = X_AVG - X_RECORD[T_RECORD] + RHS_X;
				    Y_AVG = Y_AVG - Y_RECORD[T_RECORD] + RHS_Y;
				    X_RECORD[T_RECORD] = RHS_X;
				    Y_RECORD[T_RECORD] = RHS_Y;
				    T_RECORD++;
				    if(T_RECORD == T_AVG_NUM)
					    T_RECORD = 0;
				    X[s] = X[s] + X_AVG/(T_AVG_NUM);
				    Y[s] = Y[s] + Y_AVG/(T_AVG_NUM);
			    }
		    }
		    else
		    {
			    X[s] = X[s] + RHS_X;
			    Y[s] = Y[s] + RHS_Y;
		    }
		    if(PHASE_FIXED)
			    Y[s] = 0.0;

            //output        
            if((int)(t/dt)%100 == 0)
			    cout << setiosflags(ios::fixed) << setiosflags(ios::scientific) << setprecision(6) << setiosflags(ios::showpoint) << "X: " << X[0] << "  J_dot_E_X: " << J_dot_E_X[0] << "  w[0]:" << marker[0].w << resetiosflags(ios::scientific) << resetiosflags(ios::showpoint) << resetiosflags(ios::fixed) << " N_p : " << N_p  << " time : " << t << endl;
		    if(s == 0)
                amp_file << setiosflags(ios::fixed) << setiosflags(ios::scientific) << setprecision(12) << setw(18) << X[s] << " " << setw(18) << Y[s] << resetiosflags(ios::scientific) << resetiosflags(ios::showpoint) << resetiosflags(ios::fixed) << " " << setw(2) << T_INITIA << " " << setw(2) << T_RECORD << " " << setw(9) << t;
            else
                amp_file << " " << setiosflags(ios::fixed) << setiosflags(ios::scientific) << setprecision(12) << setw(18) << X[s] << " " << setw(18) << Y[s];
        }
    }

    if(myid == 0)
        amp_file << endl;

	MPI::COMM_WORLD.Bcast(X,NUM_MODE,MPI::DOUBLE,0);
	MPI::COMM_WORLD.Bcast(Y,NUM_MODE,MPI::DOUBLE,0);




 
}

/********************************************************************************************
*                                                                                           *
*           remove particles that out of boundary and their ID=-1                           *
*                                                                                           *
********************************************************************************************/
void kick_particle()
{
	int myN_lost;
	int index_lost[deltaN];
	myN_lost = 0;
	N_lost = 0;
	for(int p=0;p<N_p;p++)
		if(marker[p].id == -1)
		{
			index_lost[myN_lost] = p;
			myN_lost++;
		}

	MPI::COMM_WORLD.Reduce(&myN_lost, &N_lost, 1, MPI::INT, MPI_SUM, 0);
	MPI::COMM_WORLD.Reduce(&N_p, &N_total, 1, MPI::INT, MPI_SUM, 0);

	if(myid ==0 && N_lost>0)
		cout << " N_lost : " << N_lost << " N_total : " << N_total << endl;

	int p=N_p-1;
	int k=0;
	for(int i=0;i<myN_lost;i++)
	{
		while(marker[p].id == -1)
		{
			p--;
			k++;
		}
		marker[index_lost[i]] = marker[p];
		if((i+1+k) >= myN_lost)
			break;
		p--;
	}
	N_p -= myN_lost;
}


double get_B_value(const particle_vector<double> &X)
{
    double B_value_par;
    double weight_line_auxiliary[2][2];
    double weight_square[2][2];
    double R,Z;
    int i,j;
                
    i = (int)((X[0]-R0+a*a_b)/dR);
    j = (int)((X[1]+a*a_b)/dZ);

    R   = R0-a*a_b+i*dR;
    Z   =   -a*a_b+j*dZ;
                
    weight_line_auxiliary[0][0] = (X[0] - R)/dR;
    weight_line_auxiliary[1][0] = (X[1] - Z)/dZ;

    weight_line_auxiliary[0][1] = 1.0 - weight_line_auxiliary[0][0];
    weight_line_auxiliary[1][1] = 1.0 - weight_line_auxiliary[1][0];

    for(int ii=0;ii<2;ii++)
	    for(int jj=0;jj<2;jj++)
		    weight_square[ii][jj] = weight_line_auxiliary[0][1-ii]*weight_line_auxiliary[1][1-jj];

    B_value_par = 0.0;
				
    for(int ii=0;ii<2;ii++)
	    for(int jj=0;jj<2;jj++)
		    B_value_par		+=	B_value[i+ii][j+jj]*weight_square[ii][jj];

	return B_value_par;
}

//add by 2011-6-20
/********************************************************************************************
*                                                                                           *
*           pitch angle scatering process                                                   *
*           lambda_new = lambda_old*(1-2*nu_d*dt)+-sqrt[(1-lambda_old^2)*2*nu_d*dt]         *
*           (using Monte Carlo method to perform "+-" randomly)                             *
*                                                                                           *
********************************************************************************************/
void scattering()
{
	//if add collison effect,mu should be changed!!
	double v,lambda_old,lambda_new;
	double B_value_par;
	double nu_d;
	const double c=0.17*v_0;
	int pm;
	#pragma omp parallel for \
		default(shared) \
		private(B_value_par,v,nu_d,lambda_old,pm,lambda_new)
	for(int p=0;p<N_p;p++)
	{
		B_value_par		=	0.5*species[0].mass*pow(marker[p].v_per,2.0)/marker[p].mu;

		v				=	sqrt(pow(marker[p].v_par,2.0)+pow(marker[p].v_per,2.0));
		nu_d			=	nu*pow(v_c,3.0)/2.0/(pow(c,3.0)+pow(v,3.0));

		lambda_old		=	marker[p].v_par/v;
		pm				=	(rand()*1.0/RAND_MAX)>0.5?1:-1; 
		lambda_new		=	lambda_old*(1-2*nu_d*dt) + pm*sqrt((1-lambda_old*lambda_old)*2*nu_d*dt);


		marker[p].v_par	=	v*lambda_new;
		marker[p].v_per	=	sqrt(v*v - pow(marker[p].v_par,2.0));

		//mu should change
		marker[p].mu	=	0.5*species[0].mass*pow(marker[p].v_per,2.0)/B_value_par;
	}
}

/********************************************************************************************
*                                                                                           *
*               particle slowing down process                                               *
*               dv/dt = -nu*[v+v_c*v_c*v_c/(v*v)]     (pitch angle is not changed)          *
*                                                                                           *
********************************************************************************************/
void slowing_down()
{
	double v,lambda;
	double B_value_par;
	#pragma omp parallel for \
		default(shared) \
		private(B_value_par,v,lambda)
	for(int p=0;p<N_p;p++)
	{
		//record old particle information
		B_value_par	=	0.5*species[0].mass*pow(marker[p].v_per,2.0)/marker[p].mu; 
		v			=	sqrt(pow(marker[p].v_par,2.0)+pow(marker[p].v_per,2.0));
		lambda		=	marker[p].v_par/v;


		//slowing down particle process
		v			=   v + (-nu*(v+pow(v_c,3.0)/pow(v,2.0)))*dt;

		//update particle information
		marker[p].v_par = lambda*v;
		marker[p].v_per = sqrt(v*v - pow(marker[p].v_par,2.0));
		
		//mu should change
		marker[p].mu	=	0.5*species[0].mass*pow(marker[p].v_per,2.0)/B_value_par;

	}
}

//add by 2011-7-30
/********************************************************************************************
*                                                                                           *
*        source and sink                                                                    *
*       (when v<v_cut,these particles are removed and reinjected back into plasma           *
*        uniformly in real space with a fixed speed v=v_max and a random pitch angle)       *
*                                                                                           *
********************************************************************************************/
void source_and_sink()
{
	double v,lambda,B_value_par;
    double weight_line_auxiliary[2][2],weight_square[2][2];
	#pragma omp parallel for \
		default(shared) \
		private(B_value_par,v,lambda,weight_line_auxiliary,weight_square)
	for(int p=0;p<N_p;p++)
	{
		v	=	sqrt(pow(marker[p].v_par,2.0)+pow(marker[p].v_per,2.0));
		if(v < v_cut)
		{
			v = v_max;
			//random pitch angle
			lambda = (rand()*1.0/RAND_MAX);
			marker[p].v_par = lambda*v;
			marker[p].v_per = sqrt(v*v - marker[p].v_par*marker[p].v_par);

			//uniform in real space
			double R2,R,Z,phi;
			do
			{
				R2  =  pow(R0-a,2.0) + (pow(R0+a,2.0)-pow(R0-a,2.0))*(rand()*1.0/RAND_MAX);
				Z   = -a + 2*a*(rand()*1.0/RAND_MAX);
				phi =  2*PI*(rand()*1.0/(RAND_MAX+1.0));
				R   =  sqrt(R2);
				marker[p].X[0] = R;
				marker[p].X[1] = Z;
				marker[p].X[2] = phi;
			}while(pow(R-R0,2.0)+pow(Z,2.0)>=a*a);
            //get B_value_par
            int i = (int)((marker[p].X[0]-R0+a*a_b)/dR);
			int j = (int)((marker[p].X[1]+a*a_b)/dZ);

			R   = R0-a*a_b+i*dR;
			Z   =   -a*a_b+j*dZ;

			weight_line_auxiliary[0][0] = (marker[p].X[0] - R)/dR;
			weight_line_auxiliary[1][0] = (marker[p].X[1] - Z)/dZ;

			weight_line_auxiliary[0][1] = 1.0 - weight_line_auxiliary[0][0];
			weight_line_auxiliary[1][1] = 1.0 - weight_line_auxiliary[1][0];

			for(int ii=0;ii<2;ii++)
				for(int jj=0;jj<2;jj++)
					weight_square[ii][jj] = weight_line_auxiliary[0][1-ii]*weight_line_auxiliary[1][1-jj];

			B_value_par = 0.0;
				
			for(int ii=0;ii<2;ii++)
				for(int jj=0;jj<2;jj++)
					B_value_par	+=	B_value[i+ii][j+jj]*weight_square[ii][jj];

			//mu should change
			marker[p].mu	=	0.5*species[0].mass*pow(marker[p].v_per,2.0)/B_value_par;

            //zero weight
            marker[p].w     =   0.0;

            //redefine g(if non-uniform loading)
//            if(!UNIFORM_LOADING)
//                marker[p].g     =   1.0/(pow(v_c,3.0) + pow(v,3.0));
        }
	}
}

/********************************************************************************************************************
*                                                                                                                   *
*     diagnostics for deltaf vs. E&P_phi to output movie    (fixed Lambda)                                          *
*                                                           (where Lambda = mu*B_0/E = 0.01,0.25,0.50,0.75,0.99)    *
*                                                           (but for anisotropic case , Lambda = mu*B_0/E = 0)      *
*                                                                                                                   *
*             NUM_TIME_DIAGNOSTICS define number of time interval region for diagnostics                            *
*                   use MPI_Gatherv function to replace Write_shared function                                       *
*                                                                                                                   *
********************************************************************************************************************/
void diagnostics_movie_use_MPI_Gather(int index_of_time_record,double Lambda_0)
{
	double Lambda,v,E,E_0,bracket_psi;//Lambda = mu*B_0/E
    double windows_width = 2e-2;
    string fh;
	stringstream ss,tt;
	ss << index_of_time_record;
    tt << setiosflags(ios::fixed) << setprecision(2) << Lambda_0;
	fh = "delta-f_vs_E&P_phi_vs_Time_" + ss.str() + "_Lambda=" + tt.str() + ".dat";
	char *filename = (char*)fh.c_str();
	ofstream deltaf_out_vs_time(filename);
    int myN_diagnostics = 0;
    int N_diagnostics;
    const int num_var = 6;


    //first step: confirm total number of particles for diagnostics every process
    for(int p=0;p<N_p;p++)
	{
		E		=	0.5*species[0].mass*pow(marker[p].v_par,2.0) + 0.5*species[0].mass*pow(marker[p].v_per,2.0);
		v		=	sqrt(2*E/species[0].mass);
        Lambda  =   marker[p].mu*B0/E;
//for anisotropic case,use 5% particles for diagnostics
        if(Lambda_0 < 1e-7)
        {
            if(p < int(N_p*0.05))
                myN_diagnostics++;
        }
        else
            if(abs(Lambda-Lambda_0) < windows_width)
                myN_diagnostics++;
	}

        
    MPI::COMM_WORLD.Reduce(&myN_diagnostics,&N_diagnostics,1,MPI::INT,MPI_SUM,0);
    MPI::COMM_WORLD.Bcast(&N_diagnostics,1,MPI::INT,0);

    //guarantee myN_diagnostics not be zero!!
    if(myN_diagnostics == 0)
        myN_diagnostics ++;
    const int mybuffer_N = myN_diagnostics;
    const int buffer_N   =   N_diagnostics;
    double myData[mybuffer_N][num_var];
    double Data[buffer_N][num_var]; 
    myN_diagnostics = 0;
    int *recvcounts  = new int[numprocs];
    int *displs      = new int[numprocs];

    for(int p=0;p<N_p;p++)
	{
		E		=	0.5*species[0].mass*pow(marker[p].v_par,2.0) + 0.5*species[0].mass*pow(marker[p].v_per,2.0);
		v		=	sqrt(2*E/species[0].mass);
        Lambda  =   marker[p].mu*B0/E;
//for anisotropic case,use 5% particles for diagnostics
        if(Lambda_0 < 1e-7 || !ISOTROPY)
        {
            if(p < int(N_p*0.05))
            {

                int sgn_v_parallel = marker[p].v_par>0?1:-1;
		        if((1-Lambda) > 0)
			        bracket_psi = marker[p].P_phi/species[0].charge - species[0].mass/species[0].charge*sgn_v_parallel*v*R0*sqrt(1-Lambda);
		        else
			        bracket_psi = marker[p].P_phi/species[0].charge;

                myData[myN_diagnostics][0] = marker[p].w*marker[p].g;
                myData[myN_diagnostics][1] = E;
                myData[myN_diagnostics][2] = marker[p].P_phi;
                myData[myN_diagnostics][3] = bracket_psi;
                myData[myN_diagnostics][4] = marker[p].X[0];
                myData[myN_diagnostics][5] = marker[p].f_over_g*marker[p].g;
                myN_diagnostics++;
            }
        }
        else
            if(abs(Lambda-Lambda_0) < windows_width)
            {
                int sgn_v_parallel = marker[p].v_par>0?1:-1;
		        if((1-Lambda) > 0)
			        bracket_psi = marker[p].P_phi/species[0].charge - species[0].mass/species[0].charge*sgn_v_parallel*v*R0*sqrt(1-Lambda);
		        else
			        bracket_psi = marker[p].P_phi/species[0].charge;

                myData[myN_diagnostics][0] = marker[p].w*marker[p].g;
                myData[myN_diagnostics][1] = E;
                myData[myN_diagnostics][2] = marker[p].P_phi;
                myData[myN_diagnostics][3] = bracket_psi;
                myData[myN_diagnostics][4] = marker[p].X[0];
                myData[myN_diagnostics][5] = marker[p].f_over_g*marker[p].g;
                myN_diagnostics++;
            }
	}

    MPI::COMM_WORLD.Gather(&myN_diagnostics,1,MPI::INT,recvcounts,1,MPI::INT,0);
    for(int i=0;i<numprocs;i++)
        recvcounts[i] = recvcounts[i]*num_var;
    displs[0] = 0;
    for(int i=1;i<numprocs;i++)
        displs[i] = displs[i-1] + recvcounts[i-1];
    
    MPI::COMM_WORLD.Gatherv(myData,myN_diagnostics*num_var,MPI::DOUBLE,Data,recvcounts,displs,MPI::DOUBLE,0);

    if(myid == 0)
    {
	    for(int i=0;i<N_diagnostics;i++)
        {
            for(int j=0;j<num_var;j++)
                deltaf_out_vs_time << Data[i][j] << " ";
            deltaf_out_vs_time << endl;
        }
    }
    delete [] recvcounts;
    delete [] displs;
    deltaf_out_vs_time.close();
}


/************************************************************************************************
*                                                                                               *
*           Calculate toroidal precession frequency : omega_phi   = Delta_phi/Delta_t           *
*                     poloidal transit    frequency : omega_theta = 2*pi/Delta_t                *
*                                                                                               *
************************************************************************************************/
bool calc_oribt_frequency(bool LOSS,bool OUTPUT_TRACK)
{
    bool cof = false;
    double Lambda;

    double eps_COF;
    if(CPN)//for perturbed case, increase eps!
        eps_COF = 3e-2;
    else
        eps_COF = 1e-2;
    for(int p=0;p<N_p;p++)
        if(marker[p].id == 0)
        {
            ofstream oribt_track("oribt_track.dat",ios::app);
            if(OUTPUT_TRACK)
                oribt_track << setiosflags(ios::fixed) << setprecision(6) << setw(9) << setiosflags(ios::internal) << marker[p].X[0] << " " << setw(9) << marker[p].X[1] << " " << setw(9) << marker[p].X[2] << " " << setw(9) << marker[p].v_par << " " << setw(9) << marker[p].v_per << " "  << setprecision(8) << setw(11) << marker[p].P_phi << " " <<setprecision(2) << setw(7) << t << endl;
            if(abs(marker[p].X[1]-SINGLE_Z) < eps_COF*a  && abs(marker[p].X[0]-SINGLE_R) < eps_COF*a && Z_reverse_sign )
            {
                ofstream oribt_frequency("oribt_frequency.dat",ios::app);
                ofstream oribt_frequency_CPN("oribt_frequency_CPN.dat",ios::app);
                
                double omega_theta,omega_phi,E;
                double omega_alpha;//including nonlinear term, frequency of phase alpha should be included! 
                omega_theta = 2*PI/t;
                omega_phi   = (marker[p].X[2]-SINGLE_phi)/t;
                E           = 0.5*species[0].mass*pow(marker[p].v_par,2.0) + 0.5*species[0].mass*pow(marker[p].v_per,2.0);
                Lambda      = marker[p].mu*B0/E;
                if(CPN)
                    omega_alpha = Delta_alpha/t;
                
                //for particle out of boundary,set orbit frequency to be zero!
                if(LOSS)
                {
                    omega_theta = 0.0;
                    omega_phi   = 0.0;
                    omega_alpha = 0.0;
                }
                if(OPS == 1)
                    oribt_frequency     << setiosflags(ios::fixed) << setprecision(8) << setw(11) << setiosflags(ios::internal) << marker[p].P_phi << " " << E << " " << omega_phi << " " << omega_theta << " " << TIMESTEP << endl;
                else if(OPS == 2)
                    oribt_frequency     << setiosflags(ios::fixed) << setprecision(8) << setw(11) << setiosflags(ios::internal) << E << " " << Lambda << " " << omega_phi << " " << omega_theta << " " << TIMESTEP << endl;

                if(CPN)
                    if(OPS == 1)
                        oribt_frequency_CPN << setiosflags(ios::fixed) << setprecision(8) << setw(11) << setiosflags(ios::internal) << marker[p].P_phi << " " << E << " " << omega_alpha << " " << Delta_alpha << " " << alpha_initial << " " << TIMESTEP << endl;
                    else
                        oribt_frequency_CPN << setiosflags(ios::fixed) << setprecision(8) << setw(11) << setiosflags(ios::internal) << E << " " << Lambda << " " << omega_alpha << " " << Delta_alpha << " " << alpha_initial << " " << TIMESTEP << endl;
               
                cof = true;
                oribt_frequency.close();
                if(CPN)
                    oribt_frequency_CPN.close();
                double EOT = -1.0;
                if(OUTPUT_TRACK)
                    oribt_track << setiosflags(ios::fixed) << setprecision(6) << setw(9) << setiosflags(ios::internal) << EOT << " " << setw(9) << EOT << " " << setw(9) << EOT << " " << setw(9) << EOT << " " << setw(9) << EOT << " "  << setprecision(8) << setw(11) << EOT << " " <<setprecision(2) << setw(7) << t << endl;
            }
            oribt_track.close();
        }
    return cof;
}

/************************************************************************************************
*                                                                                               *
*                       Numerical equilibrium form NOVA code(j-solver)                          *
*                     import data : psi(R,Z) in psi.dat                                         *
*                                   psi , g , g' , q , q' , p , p' in equilibrium_data.dat      *                  
*                                                                                               *
************************************************************************************************/
void numerical_equilibrium()
{

    int j_max,i_max;
    double R0_NOVA;
    ifstream fin("psi");

	if(!fin)
	{
		cout << "!!!!!!!!!!!psi file does not exsit!!!!!!!!!!" << endl;
		exit(1);
	}
    string str;
	getline(fin,str);//1 j = xxx i = xxx
    fin >> j_max >> i_max;
    getline(fin,str);//empty line
	getline(fin,str);//     psi          x           z

	double **psi_NOVA = new double* [(i_max-2)*2];
    double **R_NOVA   = new double* [(i_max-2)*2];
    double **Z_NOVA   = new double* [(i_max-2)*2];
	for(int i=0;i<(i_max-2)*2;i++)
    {
        psi_NOVA[i] = new double[j_max-2];
        R_NOVA[i]   = new double[j_max-2];
        Z_NOVA[i]   = new double[j_max-2];
    }

    for(int j=0;j<j_max;j++)
	    for(int i=0;i<i_max;i++)
	    {
            double R_line,Z_line,psi_line;
		    fin >> psi_line >> R_line >> Z_line;
            //psi(i,j) -- > psi(theta,r),where theta:0->pi , r:0->a
            //define R0 form NOVA code!
            if(i == 0 && j == 0)
                R0_NOVA = R_line;
            if(i!=0 && i!=i_max-1 && j!=0 && j!=j_max-1)
            {
                psi_NOVA[i-1][j-1]              =  psi_line/R0_NOVA/R0_NOVA;
                  R_NOVA[i-1][j-1]              =  R_line/R0_NOVA;
                  Z_NOVA[i-1][j-1]              =  Z_line/R0_NOVA;
                psi_NOVA[(i_max-2)*2-i][j-1]    =   psi_line/R0_NOVA/R0_NOVA;
                    R_NOVA[(i_max-2)*2-i][j-1]  =   R_line/R0_NOVA;
                    Z_NOVA[(i_max-2)*2-i][j-1]  =  -Z_line/R0_NOVA;
            }
        }
    //identify data is correct
	//output psi(R,Z)
	if(myid == 0)
	{
		ofstream psi_RZ_output("psi_RZ_NOVA.dat",ios::app);
		ofstream R_output("R_NOVA.dat",ios::app);
		ofstream Z_output("Z_NOVA.dat",ios::app);
		for(int i=0;i<(i_max-2)*2;i++)
		{
			for(int j=0;j<j_max-2;j++)
			{
				psi_RZ_output << psi_NOVA[i][j] << " ";
				R_output << R_NOVA[i][j] << " ";
				Z_output << Z_NOVA[i][j] << " ";
			}
			psi_RZ_output << endl;
			R_output << endl;
			Z_output << endl;
		}
        psi_RZ_output.close();
        R_output.close();
        Z_output.close();
	}


    if(myid ==0)
        cout << "Inputing form NOVA code is OK!!!!!!!!!!!!!" << endl;

///////////////////////////////////////////////////////////////////////////////
    psi_NOVA[0][0] = 0.0;
	for(int i=0;i<mR;i++)
		for(int j=0;j<mZ;j++)
        {
            double R   = R0-a*a_b+i*dR;
	        double Z   =   -a*a_b+j*dZ;
            double r;
            double theta;
            double theta1,theta2,theta3,theta4;
            double r1,r2,r3,r4;

            int index_theta=0,index_r=0;
            
            double err;

            err = 100;
	        for(int ii=0;ii<(i_max-2)*2;ii++)
		        for(int jj=0;jj<j_max-2;jj++)
                {
                    if(sqrt(pow(R_NOVA[ii][jj]-R,2.0)+pow(Z_NOVA[ii][jj]-Z,2.0)) < err)
                    {
                        err = sqrt(pow(R_NOVA[ii][jj]-R,2.0)+pow(Z_NOVA[ii][jj]-Z,2.0));
                        index_theta = ii;
                        index_r = jj;
                    }
                }

            r1 = sqrt(pow(R_NOVA[index_theta][index_r]-R0,2.0) + pow(Z_NOVA[index_theta][index_r],2.0));

            if(Z_NOVA[index_theta][index_r] < 0)
                theta1 = atan2(Z_NOVA[index_theta][index_r],R_NOVA[index_theta][index_r]-R0)+ 2*PI;
            else
                theta1 = atan2(Z_NOVA[index_theta][index_r],R_NOVA[index_theta][index_r]-R0);

            r = sqrt(pow(R-R0,2.0) + pow(Z,2.0));

            if(Z < 0)
                theta = atan2(Z,R-R0)+ 2*PI;
            else
                theta = atan2(Z,R-R0);

            int II,JJ;
            if(r>r1)
                JJ = index_r + 1;
            else
                JJ = index_r - 1;

            if(index_r == 0 || index_r == j_max-3)
                JJ = index_r;

            if(theta>theta1)
                II = index_theta + 1;
            else
                II = index_theta - 1;

            if(II == -1)
                II = (i_max-2)*2-1;
            else if(II == (i_max-2)*2)
                II = 0;

            double total_area;
            double area1,area2,area3,area4;

            r2 = sqrt(pow(R_NOVA[II][index_r]-R0,2.0) + pow(Z_NOVA[II][index_r],2.0));

            if(Z_NOVA[II][index_r] < 0)
                theta2 = atan2(Z_NOVA[II][index_r],R_NOVA[II][index_r]-R0)+ 2*PI;
            else
                theta2 = atan2(Z_NOVA[II][index_r],R_NOVA[II][index_r]-R0);

            r3 = sqrt(pow(R_NOVA[index_theta][JJ]-R0,2.0) + pow(Z_NOVA[index_theta][JJ],2.0));

            if(Z_NOVA[index_theta][JJ] < 0)
                theta3 = atan2(Z_NOVA[index_theta][JJ],R_NOVA[index_theta][JJ]-R0)+ 2*PI;
            else
                theta3 = atan2(Z_NOVA[index_theta][JJ],R_NOVA[index_theta][JJ]-R0);

            r4 = sqrt(pow(R_NOVA[II][JJ]-R0,2.0) + pow(Z_NOVA[II][JJ],2.0));

            if(Z_NOVA[II][JJ] < 0)
                theta4 = atan2(Z_NOVA[II][JJ],R_NOVA[II][JJ]-R0)+ 2*PI;
            else
                theta4 = atan2(Z_NOVA[II][JJ],R_NOVA[II][JJ]-R0);

            double dtheta1,dtheta2,dtheta3,dtheta4;
            dtheta1 = abs(theta1-theta);
            dtheta2 = abs(theta2-theta);
            dtheta3 = abs(theta3-theta);
            dtheta4 = abs(theta4-theta);

            if(dtheta1 > PI)
                dtheta1 = abs(dtheta1-2*PI);
            if(dtheta2 > PI)
                dtheta2 = abs(dtheta2-2*PI);
            if(dtheta3 > PI)
                dtheta3 = abs(dtheta3-2*PI);
            if(dtheta4 > PI)
                dtheta4 = abs(dtheta4-2*PI);

            area1 = abs(PI*r1*r1*(dtheta1/2*PI) - PI*r*r*(dtheta1/2*PI));
            area2 = abs(PI*r2*r2*(dtheta2/2*PI) - PI*r*r*(dtheta2/2*PI));
            area3 = abs(PI*r3*r3*(dtheta3/2*PI) - PI*r*r*(dtheta3/2*PI));
            area4 = abs(PI*r4*r4*(dtheta4/2*PI) - PI*r*r*(dtheta4/2*PI));

            total_area = area1+area2+area3+area4;

            psi[i][j] = (psi_NOVA[index_theta][index_r]*area4 + psi_NOVA[II][index_r]*area3 + psi_NOVA[index_theta][JJ]*area2 + psi_NOVA[II][JJ]*area1)/total_area;
            //psi[i][j] = psi_NOVA[index_theta][index_r];
        }

    double psi_min = 0.0;
    for(int i=0;i<mR;i++)
		for(int j=0;j<mZ;j++)
            psi_min = psi_min>psi[i][j]?psi[i][j]:psi_min;

    //psi(a) should change!!!!!
    psi_0 = psi_min;
    psi_a = 0.0;
    Delta_psi = psi_a - psi_0;
    if(myid == 0)
        cout << "!!!psi(0) form NOVA is :" << psi_min << endl;


    for(int i=0;i<(i_max-2)*2;i++)
    {
		delete [] psi_NOVA[i];
        delete [] R_NOVA[i];
        delete [] Z_NOVA[i];
    }
	delete [] psi_NOVA;
    delete [] R_NOVA;
    delete [] Z_NOVA;

    fin.close();
        
    fin.open("equilibrium_data");

	if(!fin)
	{
		cout << "!!!!!!!!!!!equilibrium data file does not exsit!!!!!!!!!!" << endl;
		exit(1);
	}

    int radial_grid = 0;
    while(getline(fin,str))
    {
        radial_grid++;
    }
    radial_grid -= 2;

    double *psi_1d = new double [radial_grid];
    double *g_1d   = new double [radial_grid];
    double *gp_1d  = new double [radial_grid];
    double *q_1d   = new double [radial_grid];
    double *qp_1d  = new double [radial_grid];
    double *p_1d   = new double [radial_grid];
    double *pp_1d  = new double [radial_grid];
    fin.clear(); 
    fin.seekg(0);
    getline(fin,str);//1  j     psi           g           gp          q           qp          p          pp
    for(int i=0;i<radial_grid;i++)
    {
        int j;
        fin >> j >> psi_1d[i] >> g_1d[i] >> gp_1d[i] >> q_1d[i] >> qp_1d[i] >> p_1d[i] >> pp_1d[i];
        psi_1d[i]   = psi_1d[i]/R0_NOVA/R0_NOVA;
        g_1d[i]     = g_1d[i];
        gp_1d[i]    = gp_1d[i]*R0_NOVA;
        qp_1d[i]    = qp_1d[i]*R0_NOVA;
        pp_1d[i]    = pp_1d[i]*R0_NOVA;
    }


    

    for(int i=0;i<mR;i++)
		for(int j=0;j<mZ;j++)
        {
            double weight;
            int index_psi;
            double err;
            err = 100;
            for(int ii=0;ii<radial_grid;ii++)
            {
                if(abs(psi[i][j]-psi_1d[ii]) < err)
                {
                    err = abs(psi[i][j]-psi_1d[ii]);
                    index_psi = ii;
                }
           }
           if(psi[i][j]-psi_1d[index_psi] > 0)
           {
               weight = (psi[i][j]-psi_1d[index_psi])/(psi_1d[index_psi+1]-psi_1d[index_psi]);
               g_eq[i][j] = g_1d[index_psi]*(1.0-weight) + g_1d[index_psi+1]*weight;
           }
           else
           {
               weight = (psi[i][j]-psi_1d[index_psi-1])/(psi_1d[index_psi]-psi_1d[index_psi-1]);
               g_eq[i][j] = g_1d[index_psi-1]*(1.0-weight) + g_1d[index_psi]*weight;
           }
           if(abs(psi[i][j]) < 1e-10)
               g_eq[i][j] = B0*R0;
        }
    //normalized by g at magnetic axis!
    for(int i=0;i<mR;i++)
		for(int j=0;j<mZ;j++)
            g_eq[i][j] = g_eq[i][j]/g_1d[0];

    delete [] psi_1d;
    delete [] g_1d;
    delete [] gp_1d;
    delete [] q_1d;
    delete [] qp_1d;
    delete [] p_1d;
    delete [] pp_1d;


}


/********************************************************************************************************************
*                                                                                                                   *
*     diagnostics for distribution function vs. P_phi with fiexd pitch angle Lambda and Energy E                    *
*                                                           (where Lambda = mu*B_0/E            = Lambda_0)         *
*                                                           (      E      = 0.5*m*(v_||)^2+mu*B = E_0	  )         *
*                                                                                                                   *
*                                                                                                                   *
********************************************************************************************************************/
void diagnostics_output_f_vs_P_phi(double Lambda_0,double E_0,double windows_width)
{
	double weight_line[2];
	int    index_P_phi;
	double P_phi;
	double Lambda,E;//Lambda = mu*B_0/E


    //get maxima and minima of P_phi at t=0
	if(TIMESTEP == 1)
	{
		double my_P_phi_min,my_P_phi_max;
		my_P_phi_min =  1e9;
        my_P_phi_max = -1e9;

		for(int p=0;p<N_p;p++)
		{
            E		=	0.5*species[0].mass*pow(marker[p].v_par,2.0) + 0.5*species[0].mass*pow(marker[p].v_per,2.0);
            Lambda  =   marker[p].mu*B0/E;

            if(ISOTROPY)
            {
                //to avoid denominator tend to be zero when Lambda_0 is too small!! 
                if(Lambda_0 < 1.5e-2)
                {
                    if(abs(Lambda-Lambda_0) < 1e-2)
			            if(abs((E-E_0)/E_0) < windows_width)
                        {
			                my_P_phi_max = marker[p].P_phi>my_P_phi_max?marker[p].P_phi:my_P_phi_max;
			                my_P_phi_min = marker[p].P_phi<my_P_phi_min?marker[p].P_phi:my_P_phi_min;
                        }
                }
                else
                {
                    if(abs((Lambda-Lambda_0)/Lambda_0) < windows_width)
			            if(abs((E-E_0)/E_0) < windows_width)
                        {
			                my_P_phi_max = marker[p].P_phi>my_P_phi_max?marker[p].P_phi:my_P_phi_max;
			                my_P_phi_min = marker[p].P_phi<my_P_phi_min?marker[p].P_phi:my_P_phi_min;
                        }
                }
            }
            else
            {
                if(abs((E-E_0)/E_0) < windows_width)
                {
			        my_P_phi_max = marker[p].P_phi>my_P_phi_max?marker[p].P_phi:my_P_phi_max;
			        my_P_phi_min = marker[p].P_phi<my_P_phi_min?marker[p].P_phi:my_P_phi_min;
                }
            }
		}

		MPI::COMM_WORLD.Reduce(&my_P_phi_max,&P_phi_max,1,MPI::DOUBLE,MPI_MAX,0);
		MPI::COMM_WORLD.Reduce(&my_P_phi_min,&P_phi_min,1,MPI::DOUBLE,MPI_MIN,0);

        MPI::COMM_WORLD.Bcast(&P_phi_max,1,MPI::DOUBLE,0);
        MPI::COMM_WORLD.Bcast(&P_phi_min,1,MPI::DOUBLE,0);

        dP_phi = (P_phi_max - P_phi_min)/(GRID_P_PHI-1);

        for(int i=0;i<GRID_P_PHI;i++)
		{
			myf_vs_P_phi[i]         = 0.0;
            mydeltaf_vs_P_phi[i]    = 0.0;
			mynum[i]                = 0.0;
		}

	}

    for(int p=0;p<N_p;p++)
	{
		E		=	0.5*species[0].mass*pow(marker[p].v_par,2.0) + 0.5*species[0].mass*pow(marker[p].v_per,2.0);
        Lambda  =   marker[p].mu*B0/E;

        
        bool Lambda_flag = false;
        //to avoid denominator tend to be zero when Lambda_0 is too small!! 
        if(Lambda_0 < 1.5e-2)
        {
            if(abs(Lambda-Lambda_0) < 1e-2)
                Lambda_flag = true;
        }
        else
        {
            if(abs((Lambda-Lambda_0)/Lambda_0) < windows_width)
                Lambda_flag = true;
        }


        if(ISOTROPY)
        {
		    if(Lambda_flag)
			    if(abs((E-E_0)/E_0) < windows_width)
				    if(marker[p].P_phi > P_phi_min && marker[p].P_phi < P_phi_max)
				    {
					    index_P_phi    = int((marker[p].P_phi-P_phi_min)/dP_phi);
					    P_phi          = P_phi_min + index_P_phi*dP_phi;
					    weight_line[0] = (marker[p].P_phi - P_phi)/dP_phi;
					    weight_line[1] = 1 - weight_line[0];

					    for(int ii=0;ii<2;ii++)
					    {
                            if(TIMESTEP <= TIME_STEP_INTERVAL+1)
                                mydeltaf_vs_P_phi[index_P_phi+ii]    += marker[p].w*marker[p].g*weight_line[1-ii];
                        
						    myf_vs_P_phi[index_P_phi+ii]    += marker[p].f_over_g*marker[p].g*weight_line[1-ii];
						    mynum[index_P_phi+ii]           += 1.0*weight_line[1-ii];
					    }
                    }
        }
        else
        {
			if(abs((E-E_0)/E_0) < windows_width)
				if(marker[p].P_phi > P_phi_min && marker[p].P_phi < P_phi_max)
				{
					index_P_phi    = int((marker[p].P_phi-P_phi_min)/dP_phi);
					P_phi          = P_phi_min + index_P_phi*dP_phi;
					weight_line[0] = (marker[p].P_phi - P_phi)/dP_phi;
					weight_line[1] = 1 - weight_line[0];

					for(int ii=0;ii<2;ii++)
					{
                        if(TIMESTEP <= TIME_STEP_INTERVAL+1)
                            mydeltaf_vs_P_phi[index_P_phi+ii]    += marker[p].w*marker[p].g*weight_line[1-ii];
                        
						myf_vs_P_phi[index_P_phi+ii]    += marker[p].f_over_g*marker[p].g*weight_line[1-ii];
						mynum[index_P_phi+ii]           += 1.0*weight_line[1-ii];
					}
                }
        }
    }


	if(TIMESTEP%TIME_STEP_INTERVAL == 0 && TIMESTEP != 1)
	{
		string fh;
		stringstream ss,tt,gg;
		ss << setiosflags(ios::fixed) << setprecision(2) << E_0/v_A/v_A;
		tt << setiosflags(ios::fixed) << setprecision(2) << Lambda_0;
        gg << setiosflags(ios::fixed) << setprecision(3) << windows_width;
		fh = "f_vs_P_phi_vs_Time_E=" + ss.str() + "_Lambda=" + tt.str() + "_windows_width=" + gg.str() + ".dat";
		char *filename = (char*)fh.c_str();
		ofstream f_out_vs_time(filename,ios::app);

		MPI::COMM_WORLD.Reduce(myf_vs_P_phi,f_vs_P_phi,GRID_P_PHI,MPI::DOUBLE,MPI_SUM,0);
		MPI::COMM_WORLD.Reduce(mynum,num,GRID_P_PHI,MPI::DOUBLE,MPI_SUM,0);
        
        if(TIMESTEP <= TIME_STEP_INTERVAL+1)
            MPI::COMM_WORLD.Reduce(mydeltaf_vs_P_phi,deltaf_vs_P_phi,GRID_P_PHI,MPI::DOUBLE,MPI_SUM,0);

		if(myid == 0)
		{
            if(TIMESTEP <= TIME_STEP_INTERVAL+1)
            {
                for(int i=0;i<GRID_P_PHI;i++)
				    f_out_vs_time << P_phi_min + i*dP_phi << " " ;
			    f_out_vs_time << 0.0 << endl;


                for(int i=0;i<GRID_P_PHI;i++)
                    if(num[i] == 0)
                        f_out_vs_time << 0.0 << " " ;
                    else
                        f_out_vs_time << (f_vs_P_phi[i]-deltaf_vs_P_phi[i])/num[i] << " " ;
			    f_out_vs_time << 0.0 << endl;
            }

			for(int i=0;i<GRID_P_PHI;i++)
                if(num[i] == 0)
                    f_out_vs_time << 0.0 << " " ;
                else
                    f_out_vs_time << f_vs_P_phi[i]/num[i] << " " ;
			f_out_vs_time << t << endl;
            
            
            double num_total = 0;
            for(int i=0;i<GRID_P_PHI;i++)
                num_total += num[i];

            cout << "Total number in 1-d f(P_phi) : " << num_total << " Lambda_0 = " << Lambda_0 << endl;
		}

		f_out_vs_time.close();



		for(int i=0;i<GRID_P_PHI;i++)
		{
			myf_vs_P_phi[i]   = 0.0;
			mynum[i] = 0.0;
		}

	}

}

void diagnostics_output_f_vs_P_phi1(double Lambda_0,double E_0,double windows_width)
{
	double weight_line[2];
	int    index_P_phi;
	double P_phi;
	double Lambda,E;//Lambda = mu*B_0/E


    //get maxima and minima of P_phi at t=0
	if(TIMESTEP == 1)
	{
		double my_P_phi_min1,my_P_phi_max1;
		my_P_phi_min1 =  1e9;
        my_P_phi_max1 = -1e9;

		for(int p=0;p<N_p;p++)
		{
            E		=	0.5*species[0].mass*pow(marker[p].v_par,2.0) + 0.5*species[0].mass*pow(marker[p].v_per,2.0);
            Lambda  =   marker[p].mu*B0/E;

            if(ISOTROPY)
            {
                //to avoid denominator tend to be zero when Lambda_0 is too small!! 
                if(Lambda_0 < 1.5e-2)
                {
                    if(abs(Lambda-Lambda_0) < 1e-2)
			            if(abs((E-E_0)/E_0) < windows_width)
                        {
			                my_P_phi_max1 = marker[p].P_phi>my_P_phi_max1?marker[p].P_phi:my_P_phi_max1;
			                my_P_phi_min1 = marker[p].P_phi<my_P_phi_min1?marker[p].P_phi:my_P_phi_min1;
                        }
                }
                else
                {
                    if(abs((Lambda-Lambda_0)/Lambda_0) < windows_width)
			            if(abs((E-E_0)/E_0) < windows_width)
                        {
			                my_P_phi_max1 = marker[p].P_phi>my_P_phi_max1?marker[p].P_phi:my_P_phi_max1;
			                my_P_phi_min1 = marker[p].P_phi<my_P_phi_min1?marker[p].P_phi:my_P_phi_min1;
                        }
                }
            }
            else
            {
                if(abs((E-E_0)/E_0) < windows_width)
                {
			        my_P_phi_max1 = marker[p].P_phi>my_P_phi_max1?marker[p].P_phi:my_P_phi_max1;
			        my_P_phi_min1 = marker[p].P_phi<my_P_phi_min1?marker[p].P_phi:my_P_phi_min1;
                }
            }
		}

		MPI::COMM_WORLD.Reduce(&my_P_phi_max1,&P_phi_max1,1,MPI::DOUBLE,MPI_MAX,0);
		MPI::COMM_WORLD.Reduce(&my_P_phi_min1,&P_phi_min1,1,MPI::DOUBLE,MPI_MIN,0);

        MPI::COMM_WORLD.Bcast(&P_phi_max1,1,MPI::DOUBLE,0);
        MPI::COMM_WORLD.Bcast(&P_phi_min1,1,MPI::DOUBLE,0);

        dP_phi1 = (P_phi_max1 - P_phi_min1)/(GRID_P_PHI-1);

        for(int i=0;i<GRID_P_PHI;i++)
		{
			myf_vs_P_phi1[i]         = 0.0;
            mydeltaf_vs_P_phi1[i]    = 0.0;
			mynum1[i]                = 0.0;
		}

    }

    for(int p=0;p<N_p;p++)
	{
		E		=	0.5*species[0].mass*pow(marker[p].v_par,2.0) + 0.5*species[0].mass*pow(marker[p].v_per,2.0);
        Lambda  =   marker[p].mu*B0/E;

        
        bool Lambda_flag = false;
        //to avoid denominator tend to be zero when Lambda_0 is too small!! 
        if(Lambda_0 < 1.5e-2)
        {
            if(abs(Lambda-Lambda_0) < 1e-2)
                Lambda_flag = true;
        }
        else
        {
            if(abs((Lambda-Lambda_0)/Lambda_0) < windows_width)
                Lambda_flag = true;
        }


        if(ISOTROPY)
        {
		    if(Lambda_flag)
			    if(abs((E-E_0)/E_0) < windows_width)
				    if(marker[p].P_phi > P_phi_min1 && marker[p].P_phi < P_phi_max1)
				    {
					    index_P_phi    = int((marker[p].P_phi-P_phi_min1)/dP_phi1);
					    P_phi          = P_phi_min1 + index_P_phi*dP_phi1;
					    weight_line[0] = (marker[p].P_phi - P_phi)/dP_phi1;
					    weight_line[1] = 1 - weight_line[0];

					    for(int ii=0;ii<2;ii++)
					    {
                            if(TIMESTEP <= TIME_STEP_INTERVAL+1)
                                mydeltaf_vs_P_phi1[index_P_phi+ii]    += marker[p].w*marker[p].g*weight_line[1-ii];
                        
						    myf_vs_P_phi1[index_P_phi+ii]    += marker[p].f_over_g*marker[p].g*weight_line[1-ii];
						    mynum1[index_P_phi+ii]           += 1.0*weight_line[1-ii];
					    }
                    }
        }
        else
        {
			if(abs((E-E_0)/E_0) < windows_width)
				if(marker[p].P_phi > P_phi_min1 && marker[p].P_phi < P_phi_max1)
				{
					index_P_phi    = int((marker[p].P_phi-P_phi_min1)/dP_phi1);
					P_phi          = P_phi_min1 + index_P_phi*dP_phi1;
					weight_line[0] = (marker[p].P_phi - P_phi)/dP_phi1;
					weight_line[1] = 1 - weight_line[0];

					for(int ii=0;ii<2;ii++)
					{
                        if(TIMESTEP <= TIME_STEP_INTERVAL+1)
                            mydeltaf_vs_P_phi1[index_P_phi+ii]    += marker[p].w*marker[p].g*weight_line[1-ii];
                        
						myf_vs_P_phi1[index_P_phi+ii]    += marker[p].f_over_g*marker[p].g*weight_line[1-ii];
						mynum1[index_P_phi+ii]           += 1.0*weight_line[1-ii];
					}
                }
        }
    }


	if(TIMESTEP%TIME_STEP_INTERVAL == 0 && TIMESTEP != 1)
	{
		string fh;
		stringstream ss,tt,gg;
		ss << setiosflags(ios::fixed) << setprecision(2) << E_0/v_A/v_A;
		tt << setiosflags(ios::fixed) << setprecision(2) << Lambda_0;
        gg << setiosflags(ios::fixed) << setprecision(3) << windows_width;
		fh = "f_vs_P_phi_vs_Time_E=" + ss.str() + "_Lambda=" + tt.str() + "_windows_width=" + gg.str() + ".dat";
		char *filename = (char*)fh.c_str();
		ofstream f_out_vs_time(filename,ios::app);

		MPI::COMM_WORLD.Reduce(myf_vs_P_phi1,f_vs_P_phi1,GRID_P_PHI,MPI::DOUBLE,MPI_SUM,0);
		MPI::COMM_WORLD.Reduce(mynum1,num1,GRID_P_PHI,MPI::DOUBLE,MPI_SUM,0);
        
        if(TIMESTEP <= TIME_STEP_INTERVAL+1)
            MPI::COMM_WORLD.Reduce(mydeltaf_vs_P_phi1,deltaf_vs_P_phi1,GRID_P_PHI,MPI::DOUBLE,MPI_SUM,0);

		if(myid == 0)
		{
            if(TIMESTEP <= TIME_STEP_INTERVAL+1)
            {
                for(int i=0;i<GRID_P_PHI;i++)
				    f_out_vs_time << P_phi_min1 + i*dP_phi1 << " " ;
			    f_out_vs_time << 0.0 << endl;


                for(int i=0;i<GRID_P_PHI;i++)
                    if(num1[i] == 0)
                        f_out_vs_time << 0.0 << " " ;
                    else
                        f_out_vs_time << (f_vs_P_phi1[i]-deltaf_vs_P_phi1[i])/num1[i] << " " ;
			    f_out_vs_time << 0.0 << endl;
            }

			for(int i=0;i<GRID_P_PHI;i++)
                if(num1[i] == 0)
                    f_out_vs_time << 0.0 << " " ;
                else
                    f_out_vs_time << f_vs_P_phi1[i]/num1[i] << " " ;
			f_out_vs_time << t << endl;
            
            
            double num_total = 0;
            for(int i=0;i<GRID_P_PHI;i++)
                num_total += num1[i];

            cout << "Total number in 1-d f(P_phi) : " << num_total << " Lambda_0 = " << Lambda_0 << endl;
		}

		f_out_vs_time.close();



		for(int i=0;i<GRID_P_PHI;i++)
		{
			myf_vs_P_phi1[i]   = 0.0;
			mynum1[i] = 0.0;
		}

	}

}


/********************************************************************************************************************
*                                                                                                                   *
*     diagnostics for distribution function vs. P_phi and Energy E with fiexd pitch angle Lambda                    *
*                                                           (where Lambda = mu*B_0/E            = Lambda_0)         *
*                                                           (      E      = 0.5*m*(v_||)^2+mu*B = E_0	  )         *
*                                                                                                                   *
*    #file structrue : f_vs_P_phi_and_E_vs_Time_Lambda=XXX_windows_width=XXX.dat                                    *
*                         E E E                                                                                     *
*                       P X X X                                                                                     *
*                       P X X X                                                                                     *
*                       P X X X                                                                                     *
*                       p O O O                                                                                     *
*                       p O O O                                                                                     *
*                       P O O O                                                                                     *
*                                                                                                                   *
*            here    :  E : Energy  P : P_phi    X : t=t0    O : t=t0+dt(next time step)                            *
*                                                                                                                   *
*  add by 2012.9.29                                                                                                 *
*  if DPK is true, it will output f(P_phi,K), where K = ~P_phi-n/omega*~E (actually,I replace E with K directly)    *
*                                                                                                                   *
*  add by 2013.12.27                                                                                                *
*  if IAL is true, it will output f(P_phi,E) and deltaf(P_phi,E) by Intergrating All Lambda                         *                                                                                                                *
*                                                                                                                   *
********************************************************************************************************************/
void diagnostics_output_f_vs_P_phi_and_E(double Lambda_0,double windows_width)
{
	double weight_line[2][2];
    double weight_square[2][2];
	int    index_P_phi,index_E;
	double P_phi_floor,E_floor;
	double Lambda,E;//Lambda = mu*B_0/E


    //get maxima and minima of P_phi and Energy at t=0
	if(TIMESTEP == 1)
	{
		double my_P_phi_min,my_P_phi_max;
        double my_E_min,my_E_max;
        my_P_phi_min =  1e9;
        my_P_phi_max = -1e9;
        my_E_min     =  1e9;
        my_E_max     = -1e9;

		for(int p=0;p<N_p;p++)
		{
            E		=	0.5*species[0].mass*pow(marker[p].v_par,2.0) + 0.5*species[0].mass*pow(marker[p].v_per,2.0);
            Lambda  =   marker[p].mu*B0/E;
            if(DPK)
                E		=	marker[p].K;

            if(ISOTROPY)
            {
                //Intergrating All Lambda!!
                if(IAL)
                {
                    my_P_phi_max = marker[p].P_phi>my_P_phi_max?marker[p].P_phi:my_P_phi_max;
			        my_P_phi_min = marker[p].P_phi<my_P_phi_min?marker[p].P_phi:my_P_phi_min;
                    my_E_max     = E>my_E_max?E:my_E_max;
			        my_E_min     = E<my_E_min?E:my_E_min;
                }
                else
                {
                    if(Lambda_0 < 1.5e-2)
                    {
                        //to avoid denominator tend to be zero when Lambda_0 is too small!! 
                        if(abs(Lambda-Lambda_0) < 1e-2)
                        {
			                my_P_phi_max = marker[p].P_phi>my_P_phi_max?marker[p].P_phi:my_P_phi_max;
			                my_P_phi_min = marker[p].P_phi<my_P_phi_min?marker[p].P_phi:my_P_phi_min;
                            my_E_max     = E>my_E_max?E:my_E_max;
			                my_E_min     = E<my_E_min?E:my_E_min;
                        }
                    }
                    else
                    {
                        if(abs((Lambda-Lambda_0)/Lambda_0) < windows_width)
                        {
			                my_P_phi_max = marker[p].P_phi>my_P_phi_max?marker[p].P_phi:my_P_phi_max;
			                my_P_phi_min = marker[p].P_phi<my_P_phi_min?marker[p].P_phi:my_P_phi_min;
                            my_E_max     = E>my_E_max?E:my_E_max;
			                my_E_min     = E<my_E_min?E:my_E_min;
                        }
                    }
                }
            }
            else
            {
			    my_P_phi_max = marker[p].P_phi>my_P_phi_max?marker[p].P_phi:my_P_phi_max;
			    my_P_phi_min = marker[p].P_phi<my_P_phi_min?marker[p].P_phi:my_P_phi_min;
                my_E_max     = E>my_E_max?E:my_E_max;
			    my_E_min     = E<my_E_min?E:my_E_min;
            }
		}

		MPI::COMM_WORLD.Reduce(&my_P_phi_max,&P_phi_max,1,MPI::DOUBLE,MPI_MAX,0);
		MPI::COMM_WORLD.Reduce(&my_P_phi_min,&P_phi_min,1,MPI::DOUBLE,MPI_MIN,0);

        MPI::COMM_WORLD.Bcast(&P_phi_max,1,MPI::DOUBLE,0);
        MPI::COMM_WORLD.Bcast(&P_phi_min,1,MPI::DOUBLE,0);

        dP_phi = (P_phi_max - P_phi_min)/(GRID_P_PHI-1);

        MPI::COMM_WORLD.Reduce(&my_E_max,&E_max,1,MPI::DOUBLE,MPI_MAX,0);
		MPI::COMM_WORLD.Reduce(&my_E_min,&E_min,1,MPI::DOUBLE,MPI_MIN,0);

        MPI::COMM_WORLD.Bcast(&E_max,1,MPI::DOUBLE,0);
        MPI::COMM_WORLD.Bcast(&E_min,1,MPI::DOUBLE,0);

        dE = (E_max - E_min)/(GRID_E-1);

    }


    for(int i=0;i<GRID_P_PHI;i++)
        for(int j=0;j<GRID_E;j++)
		{
			myf_vs_P_phi_and_E[i][j]         = 0.0;
            mydeltaf_vs_P_phi_and_E[i][j]    = 0.0;
			mynum_2d[i][j]                   = 0.0;
		}

    for(int p=0;p<N_p;p++)
	{
	    E		=	0.5*species[0].mass*pow(marker[p].v_par,2.0) + 0.5*species[0].mass*pow(marker[p].v_per,2.0);
        Lambda  =   marker[p].mu*B0/E;
        if(DPK)
            E		=	marker[p].K;

        
        bool Lambda_flag = false;
        if(IAL)
            Lambda_flag = true;
        else
        {
            //to avoid denominator tend to be zero when Lambda_0 is too small!! 
            if(Lambda_0 < 1.5e-2)
            {
                if(abs(Lambda-Lambda_0) < 1e-2)
                    Lambda_flag = true;
            }
            else
            {
                if(abs((Lambda-Lambda_0)/Lambda_0) < windows_width)
                    Lambda_flag = true;
            }
        }


        if(ISOTROPY)
        {
		    if(Lambda_flag)
			    if(E > E_min && E < E_max)
				    if(marker[p].P_phi > P_phi_min && marker[p].P_phi < P_phi_max)
				    {
					    index_P_phi        = int((marker[p].P_phi-P_phi_min)/dP_phi);
					    P_phi_floor        = P_phi_min + index_P_phi*dP_phi;
					    weight_line[0][0]  = (marker[p].P_phi - P_phi_floor)/dP_phi;
					    weight_line[0][1]  = 1 - weight_line[0][0];
                        index_E            = int((E-E_min)/dE);
					    E_floor            = E_min + index_E*dE;
                        weight_line[1][0]  = (E - E_floor)/dE;
					    weight_line[1][1]  = 1 - weight_line[1][0];

					    for(int ii=0;ii<2;ii++)
                            for(int jj=0;jj<2;jj++)
					        {
                                weight_square[ii][jj]                                      = weight_line[0][1-ii]*weight_line[1][1-jj];

                                mydeltaf_vs_P_phi_and_E[index_P_phi+ii][index_E+jj]       +=        marker[p].w*marker[p].g*weight_square[ii][jj];
                                     myf_vs_P_phi_and_E[index_P_phi+ii][index_E+jj]       += marker[p].f_over_g*marker[p].g*weight_square[ii][jj];
						                       mynum_2d[index_P_phi+ii][index_E+jj]       +=                            1.0*weight_square[ii][jj];
					        }
                    }
        }
        else
        {
			if(E > E_min && E < E_max)
				if(marker[p].P_phi > P_phi_min && marker[p].P_phi < P_phi_max)
				{
					index_P_phi        = int((marker[p].P_phi-P_phi_min)/dP_phi);
					P_phi_floor        = P_phi_min + index_P_phi*dP_phi;
					weight_line[0][0]  = (marker[p].P_phi - P_phi_floor)/dP_phi;
					weight_line[0][1]  = 1 - weight_line[0][0];
                    index_E            = int((E-E_min)/dE);
					E_floor            = E_min + index_E*dE;
                    weight_line[1][0]  = (E - E_floor)/dE;
					weight_line[1][1]  = 1 - weight_line[1][0];

					for(int ii=0;ii<2;ii++)
                        for(int jj=0;jj<2;jj++)
					    {
                            weight_square[ii][jj]                                      = weight_line[0][1-ii]*weight_line[1][1-jj];

                            mydeltaf_vs_P_phi_and_E[index_P_phi+ii][index_E+jj]       +=        marker[p].w*marker[p].g*weight_square[ii][jj];
                                    myf_vs_P_phi_and_E[index_P_phi+ii][index_E+jj]    += marker[p].f_over_g*marker[p].g*weight_square[ii][jj];
						                    mynum_2d[index_P_phi+ii][index_E+jj]      +=                            1.0*weight_square[ii][jj];
					    }
                }
        }
    }


	if(TIMESTEP%TIME_STEP_INTERVAL == 0 && TIMESTEP != 1)
	{
		string fh,fh_P,fh_E,fh_t,fh_d;
		stringstream tt,gg;
		tt << setiosflags(ios::fixed) << setprecision(2) << Lambda_0;
        gg << setiosflags(ios::fixed) << setprecision(3) << windows_width;
        if(DPK)
        {
            fh   = "f_vs_P_phi_and_K_vs_Time_Lambda=" + tt.str() + "_windows_width=" + gg.str() + ".dat";
            fh_P = "f_vs_P_phi_and_K_vs_Time_list_of_P_phi.dat";
            fh_E = "f_vs_P_phi_and_K_vs_Time_list_of_K.dat";
            fh_t = "f_vs_P_phi_and_K_vs_Time_list_of_time.dat";
        }
        else if(IAL)
        {
		    fh   = "f_vs_P_phi_and_E_vs_Time_AllLambda.dat";
            fh_d = "deltaf_vs_P_phi_and_E_vs_Time_AllLambda.dat";
            fh_P = "f_vs_P_phi_and_E_vs_Time_list_of_P_phi.dat";
            fh_E = "f_vs_P_phi_and_E_vs_Time_list_of_E.dat";
            fh_t = "f_vs_P_phi_and_E_vs_Time_list_of_time.dat";
        }
        else
        {
		    fh   =      "f_vs_P_phi_and_E_vs_Time_Lambda=" + tt.str() + "_windows_width=" + gg.str() + ".dat";
            fh_d = "deltaf_vs_P_phi_and_E_vs_Time_Lambda=" + tt.str() + "_windows_width=" + gg.str() + ".dat";
            fh_P = "f_vs_P_phi_and_E_vs_Time_list_of_P_phi.dat";
            fh_E = "f_vs_P_phi_and_E_vs_Time_list_of_E.dat";
            fh_t = "f_vs_P_phi_and_E_vs_Time_list_of_time.dat";
        }
		char *filename  = (char*)fh.c_str();
        char *filename1 = (char*)fh_P.c_str();
        char *filename2 = (char*)fh_E.c_str();
        char *filename3 = (char*)fh_t.c_str();
        char *filename4 = (char*)fh_d.c_str();
		ofstream f_out_vs_time(filename,ios::app);
        ofstream f_out_P_phi(filename1,ios::app);
        ofstream f_out_E(filename2,ios::app);
        ofstream f_out_t(filename3,ios::app);
        ofstream deltaf_out_vs_time(filename4,ios::app);

		MPI::COMM_WORLD.Reduce(myf_vs_P_phi_and_E,f_vs_P_phi_and_E,GRID_P_PHI*GRID_E,MPI::DOUBLE,MPI_SUM,0);
		MPI::COMM_WORLD.Reduce(mynum_2d,num_2d,GRID_P_PHI*GRID_E,MPI::DOUBLE,MPI_SUM,0);
        
        MPI::COMM_WORLD.Reduce(mydeltaf_vs_P_phi_and_E,deltaf_vs_P_phi_and_E,GRID_P_PHI*GRID_E,MPI::DOUBLE,MPI_SUM,0);

		if(myid == 0)
		{
            if(TIMESTEP <= TIME_STEP_INTERVAL+1)
            {
                for(int i=0;i<GRID_P_PHI;i++)
				    f_out_P_phi << P_phi_min + i*dP_phi << " " ;
                f_out_P_phi.close();

                for(int i=0;i<GRID_E;i++)
				    f_out_E << E_min + i*dE << " " ;
                f_out_E.close();


                for(int i=0;i<GRID_P_PHI;i++)
                {
                    for(int j=0;j<GRID_E;j++)
                    {
                        if(num_2d[i][j] == 0)
                            f_out_vs_time << 0.0 << " " ;
                        else
                            f_out_vs_time << (f_vs_P_phi_and_E[i][j]-deltaf_vs_P_phi_and_E[i][j])/num_2d[i][j] << " " ;
                    }
                    f_out_vs_time << endl;
                }
                //output initia time
                f_out_t << 0.0 << " " ;
                f_out_t.close();
            }

			for(int i=0;i<GRID_P_PHI;i++)
            {
                for(int j=0;j<GRID_E;j++)
                {
                    if(num_2d[i][j] == 0)
                        f_out_vs_time << 0.0 << " " ;
                    else
                        f_out_vs_time << f_vs_P_phi_and_E[i][j]/num_2d[i][j] << " " ;
                }
			    f_out_vs_time << endl;
            }

            //output deltaf vs time add by 2103.12.27
            for(int i=0;i<GRID_P_PHI;i++)
            {
                for(int j=0;j<GRID_E;j++)
                {
                    if(num_2d[i][j] == 0)
                        deltaf_out_vs_time << 0.0 << " " ;
                    else
                        deltaf_out_vs_time << deltaf_vs_P_phi_and_E[i][j]/num_2d[i][j] << " " ;
                }
			    deltaf_out_vs_time << endl;
            }


            f_out_t << t << " " ;
            f_out_t.close();
            
            
            double num_total = 0;
            for(int i=0;i<GRID_P_PHI;i++)
                num_total += num_2d[i][int(GRID_E/2)];

            cout << "Total number in 1-d f(P_phi) : " << num_total << " Lambda_0 = " << Lambda_0 << endl;
		}

		f_out_vs_time.close();
        deltaf_out_vs_time.close();

	}

}

/********************************************************************************************************************
*                                                                                                                   *
*               diagnostics for number of particles with specific pitch angle Lambda                                *
*                                                                                                                   *
********************************************************************************************************************/
void diagnostics_on_number_of_particle()
{
    double Lambda,E;
    int my_number_Lambda_1 = 0;
    int my_number_Lambda_2 = 0;
    int number_Lambda_1,number_Lambda_2;

    for(int p=0;p<N_p;p++)
    {
        E		=	0.5*species[0].mass*pow(marker[p].v_par,2.0) + 0.5*species[0].mass*pow(marker[p].v_per,2.0);
        Lambda  =   marker[p].mu*B0/E;
        if(Lambda >  Lambda_0_load_1-0.01 && Lambda <  Lambda_0_load_1+0.01)
            my_number_Lambda_1 ++ ;
        if(Lambda >  Lambda_0_load_2-0.01 && Lambda <  Lambda_0_load_2+0.01)
            my_number_Lambda_2 ++ ;

        MPI::COMM_WORLD.Reduce(&my_number_Lambda_1,&number_Lambda_1,1,MPI::INT,MPI_SUM,0);
        MPI::COMM_WORLD.Reduce(&my_number_Lambda_2,&number_Lambda_2,1,MPI::INT,MPI_SUM,0);
    }

    if(myid == 0)
    {
        cout << "==========Diagnostics for number of particles with specific pitch angle Lambda==========" << endl;
        cout << "Number of particles for with specific pitch angle (Lambda=" << setiosflags(ios::fixed)  << setprecision(2) << Lambda_0_load_1 << ") : " << number_Lambda_1 << endl;
        cout << "Number of particles for with specific pitch angle (Lambda=" << setiosflags(ios::fixed)  << setprecision(2) << Lambda_0_load_2 << ") : " << number_Lambda_2 << endl;
        cout << "Normalization factors for beta of energetic particles (c_f)     : " << setiosflags(ios::fixed)  << setprecision(6) << setiosflags(ios::scientific) << c_f << endl;
        cout << "========================================================================================" << endl;
    }
}

/********************************************************************************************************************
*                                                                                                                   *
*                  load amplitude.dat file into code and do test particle simulation                                *
*                                                                                                                   *
********************************************************************************************************************/
int load_amp_file()
{
    ifstream fin("amp_file/amplitude.dat");
    ifstream fin1("amp_file/amplitude.dat");
    string str;
    int num_line = 0;
	if(!fin)
	{
		cout << "!!!!!!!!!!!input amplitude file does not exsit!!!!!!!!!!" << endl;
		exit(1);
	}
    else
    {
        while(getline(fin,str))
        {
            num_line ++;
        }
    }
    fin.clear();
    fin.seekg(0,fin.beg);


    X_input = new double[num_line];
    Y_input = new double[num_line];
    t_input = new double[num_line];

    for(int i=0;i<num_line;i++)
    {
        double temp1,temp2;
        fin >> X_input[i] >> Y_input[i] >> temp1 >> temp2 >> t_input[i];
    }
    fin.close();
    return num_line;

}


/************************************************************************************
*                                                                                   *
*                  plot Time Evolution of KAM surfaces                              *
*                                                                                   *
*       restrictions for particles output : Z      = 0(theta = 0)                   *
*                                           K      = n/omega*E+P_phi = const        *
*                                           Lambda = SINGLE_Lambda                  *
*                        plot KAM surface : PHASE  = n*phi -omega*t                 *
*                                                                                   *
*           p.s. use MPI_Gatherv function to replace Write_shared function          *
*                                                                                   *
*   file structure : KAM_vs_Time_XX_Lambda\=0.XX.dat                                *
*                    R   P_phi   v_||   Phase   phi                                 *
*                                                                                   *
************************************************************************************/
void output_Time_Evolution_of_KAM(int index_of_time_record)
{
    string fh;
	stringstream ss,tt;
	ss << index_of_time_record;
    tt << setiosflags(ios::fixed) << setprecision(2) << SINGLE_Lambda;
	fh = "KAM_vs_Time_" + ss.str() + "_Lambda=" + tt.str() + ".dat";
	char *filename = (char*)fh.c_str();
	ofstream file_KAM_t(filename);
    int myN_diagnostics = 0;
    int N_diagnostics;
    const int num_var = 5;

//first step: confirm total number of particles for diagnostics every process
    for(int p=0;p<N_p;p++)
		for(int l=0;l<N_p_diagnosis;l++)
			if(marker[p].id == res_index[l])
//				if((marker[p].X[1] >0 && marker[p].w <0 && marker[p].X[0] > R0) || (marker[p].X[1] <0 && marker[p].w >0 && marker[p].X[0] > R0)) 
// for trapped particles, marker[p].X[0] is always greater than R0, so condition should be modified!!! //add by 2014-09-19
                if(marker[p].X[1] >0 && marker[p].w <0)
                    myN_diagnostics++;
    
    MPI::COMM_WORLD.Reduce(&myN_diagnostics,&N_diagnostics,1,MPI::INT,MPI_SUM,0);
    MPI::COMM_WORLD.Bcast(&N_diagnostics,1,MPI::INT,0);

    //guarantee myN_diagnostics not be zero!!
    if(myN_diagnostics == 0)
        myN_diagnostics ++;
    const int mybuffer_N = myN_diagnostics;
    const int buffer_N   =   N_diagnostics;
    double myData[mybuffer_N][num_var];
    double Data[buffer_N][num_var]; 
    myN_diagnostics = 0;
    int *recvcounts  = new int[numprocs];
    int *displs      = new int[numprocs];


    for(int p=0;p<N_p;p++)
		for(int l=0;l<N_p_diagnosis;l++)
			if(marker[p].id == res_index[l])
//				if((marker[p].X[1] >0 && marker[p].w <0 && marker[p].X[0] > R0) || (marker[p].X[1] <0 && marker[p].w >0 && marker[p].X[0] > R0)) 
// for trapped particles, marker[p].X[0] is always greater than R0, so condition should be modified!!! //add by 2014-09-19
                if(marker[p].X[1] >0 && marker[p].w <0)
				{
					double phase = n[0]*marker[p].X[2]-omega_A[0]*t;
					while(phase>PI)
						phase = phase - 2*PI;
					while(phase<-PI)
						phase = phase + 2*PI;
/************************************************************************************
*                         output KAM surface v.s. time (2013.10.26)                 *
************************************************************************************/
                    myData[myN_diagnostics][0] = marker[p].X[0];
                    myData[myN_diagnostics][1] = marker[p].P_phi;
                    myData[myN_diagnostics][2] = marker[p].v_par;
                    myData[myN_diagnostics][3] = phase;
                    myData[myN_diagnostics][4] = marker[p].X[2];
                    myN_diagnostics++;
                }
    
                            
    MPI::COMM_WORLD.Gather(&myN_diagnostics,1,MPI::INT,recvcounts,1,MPI::INT,0);
    for(int i=0;i<numprocs;i++)
        recvcounts[i] = recvcounts[i]*num_var;
    displs[0] = 0;
    for(int i=1;i<numprocs;i++)
        displs[i] = displs[i-1] + recvcounts[i-1];
    
    MPI::COMM_WORLD.Gatherv(myData,myN_diagnostics*num_var,MPI::DOUBLE,Data,recvcounts,displs,MPI::DOUBLE,0);

    if(myid == 0)
    {
	    for(int i=0;i<N_diagnostics;i++)
        {
            for(int j=0;j<num_var;j++)
                file_KAM_t << setiosflags(ios::fixed) << setprecision(12) << Data[i][j] << " ";
            file_KAM_t << endl;
        }
    }
    delete [] recvcounts;
    delete [] displs;

    file_KAM_t.close();

}

/*************************************************************************************
*             output perturbed field(deltaB_R,deltaE_R) v.s. time                    *
*              at position (R_field,Z_field,phi_field)  (2013.10.29)                 *
*************************************************************************************/
void output_perturbed_field_vs_time(double R_field,double Z_field,double phi_field)
{
    particle_vector<double> B_par,b_par,deltaB_par,deltaE_par,curl_b_par,grad_B_par;
	double B_value_par,psi_par;
	double deltaPhi_par,deltaA_par_par;
	particle_vector<double> grad_deltaA_par_par,grad_deltaPhi_par;
	double weight_square[2][2];
	double weight_line[3][2];
	double weight_line_aux[2];
	double weight_line_auxiliary[2][2];
    int i,j;

	double R,Z,phi;


    i = (int)((R_field-R0+a*a_b)/dR);
	j = (int)((Z_field+a*a_b)/dZ);



	R   = R0-a*a_b+i*dR;
	Z   =   -a*a_b+j*dZ;


 
	weight_line_auxiliary[0][0] = (R_field - R)/dR;
	weight_line_auxiliary[1][0] = (Z_field - Z)/dZ;

	weight_line_auxiliary[0][1] = 1.0 - weight_line_auxiliary[0][0];
	weight_line_auxiliary[1][1] = 1.0 - weight_line_auxiliary[1][0];
				

	for(int ii=0;ii<2;ii++)
		for(int jj=0;jj<2;jj++)
			weight_square[ii][jj] = weight_line_auxiliary[0][1-ii]*weight_line_auxiliary[1][1-jj];


	B_par       = 0.0;
	curl_b_par  = 0.0;
	grad_B_par  = 0.0;
	B_value_par = 0.0;
				
	for(int ii=0;ii<2;ii++)
		for(int jj=0;jj<2;jj++)
		{
			B_value_par		+=	B_value[i+ii][j+jj]*weight_square[ii][jj];
			for(int d=0;d<2;d++)
			{
				B_par[d]      +=      B[d][i+ii][j+jj]*weight_square[ii][jj];
				grad_B_par[d] += grad_B[d][i+ii][j+jj]*weight_square[ii][jj];
			}

			for(int d=0;d<3;d++)
				curl_b_par[d] += curl_b[d][i+ii][j+jj]*weight_square[ii][jj];
		}

	if(NUMERICAL_EQUILIBRIUM)
    {
        double g_eq_par = 0.0;
        for(int ii=0;ii<2;ii++)
			for(int jj=0;jj<2;jj++)
				g_eq_par		+=	g_eq[i+ii][j+jj]*weight_square[ii][jj];

        B_par[2] = g_eq_par/R_field;
    }
    else
        B_par[2] = B0*R0/R_field;

    grad_B_par[2] = 0.0;


	for(int d=0;d<3;d++)
		b_par[d] =   B_par[d]/B_value_par;

	grad_deltaA_par_par = 0.0;
	grad_deltaPhi_par   = 0.0;
	deltaA_par_par      = 0.0;
    for(int s=0;s<NUM_MODE;s++)
	{
		COS_NPHI_PLUS_OMEGA_T = cos(n[s]*phi_field-omega_A[s]*t);
		SIN_NPHI_PLUS_OMEGA_T = sin(n[s]*phi_field-omega_A[s]*t);

		for(int d=0;d<3;d++)
			for(int I=i;I<=i+1;I++)
				for(int J=j;J<=j+1;J++)
				{
					if(d!=2)
						grad_deltaA_par_local[d][I-i][J-j]  =       X[s]*grad_deltaA_par_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T - X[s]*grad_deltaA_par_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaA_par_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaA_par_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T;
					else
						grad_deltaA_par_local[d][I-i][J-j]  =   (-  X[s]*grad_deltaA_par_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T - X[s]*grad_deltaA_par_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T - Y[s]*grad_deltaA_par_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaA_par_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T)*n[s]/R_field;
					if(d!=2)
						grad_deltaPhi_local[d][I-i][J-j]    =       X[s]*grad_deltaPhi_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T - X[s]*grad_deltaPhi_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaPhi_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaPhi_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T;
					else
						grad_deltaPhi_local[d][I-i][J-j]    =   (-  X[s]*grad_deltaPhi_cos[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T - X[s]*grad_deltaPhi_sin[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T - Y[s]*grad_deltaPhi_sin[s][d][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*grad_deltaPhi_cos[s][d][I][J]*COS_NPHI_PLUS_OMEGA_T)*n[s]/R_field;

				}
		for(int I=i;I<=i+1;I++)
			for(int J=j;J<=j+1;J++)
				deltaA_par_local[I-i][J-j]  =   X[s]*deltaA_par_cos[s][I][J]*COS_NPHI_PLUS_OMEGA_T - X[s]*deltaA_par_sin[s][I][J]*SIN_NPHI_PLUS_OMEGA_T + Y[s]*deltaA_par_sin[s][I][J]*COS_NPHI_PLUS_OMEGA_T + Y[s]*deltaA_par_cos[s][I][J]*SIN_NPHI_PLUS_OMEGA_T;



				
		for(int ii=0;ii<2;ii++)
			for(int jj=0;jj<2;jj++)
			{
				for(int d=0;d<3;d++)
				{
					grad_deltaA_par_par[d]	+=	grad_deltaA_par_local[d][ii][jj]*weight_square[ii][jj];
					grad_deltaPhi_par[d]	+=	grad_deltaPhi_local[d][ii][jj]*weight_square[ii][jj];
				}
				deltaA_par_par += deltaA_par_local[ii][jj]*weight_square[ii][jj];
			}


    }
	deltaB_par = times(grad_deltaA_par_par,b_par) + deltaA_par_par*curl_b_par;
	deltaE_par = (b_par*grad_deltaPhi_par)*b_par - grad_deltaPhi_par ;

    //out_put field deltaB_R,deltaE_R
    ofstream file_perturbed_field_t("perturbed_field_vs_time.dat",ios::app);
    file_perturbed_field_t << setiosflags(ios::fixed) << setiosflags(ios::scientific) << setprecision(12) << setw(18) <<deltaB_par[0] << " " << setw(18) << deltaE_par[0] << resetiosflags(ios::scientific) << resetiosflags(ios::showpoint) << resetiosflags(ios::fixed) << " " << setw(9) << t << endl;

}


/********************************************************************************************************************
*                                                                                                                   *
*     diagnostics for deltaf vs. P_phi&Lambda to output movie    (fixed E(Energy))                                  *
*                                                                                                                   *
*             NUM_TIME_DIAGNOSTICS define number of time interval region for diagnostics                            *
*                   use MPI_Gatherv function to replace Write_shared function      (2013.11.21)                     *
*                                                                                                                   *
********************************************************************************************************************/
void diagnostics_output_deltaf_vs_P_phi_and_Lambda(int index_of_time_record,double E_0,double windows_width)
{
	double Lambda,v,E,bracket_psi;//Lambda = mu*B_0/E
    string fh;
	stringstream ss,tt;
	ss << index_of_time_record;
    tt << setiosflags(ios::fixed) << setprecision(2) << E_0/v_A/v_A;
	fh = "delta-f_vs_P_phi&Lambda_vs_Time_" + ss.str() + "_E0=" + tt.str() + ".dat";
	char *filename = (char*)fh.c_str();
	ofstream deltaf_out_vs_time(filename);
    int myN_diagnostics = 0;
    int N_diagnostics;
    const int num_var = 6;


    //first step: confirm total number of particles for diagnostics every process
    for(int p=0;p<N_p;p++)
	{
		E		=	0.5*species[0].mass*pow(marker[p].v_par,2.0) + 0.5*species[0].mass*pow(marker[p].v_per,2.0);
		v		=	sqrt(2*E/species[0].mass);
        Lambda  =   marker[p].mu*B0/E;

        if(abs((E-E_0)/E) < windows_width)
            myN_diagnostics++;
	}

        
    MPI::COMM_WORLD.Reduce(&myN_diagnostics,&N_diagnostics,1,MPI::INT,MPI_SUM,0);
    MPI::COMM_WORLD.Bcast(&N_diagnostics,1,MPI::INT,0);

    //guarantee myN_diagnostics not be zero!!
    if(myN_diagnostics == 0)
        myN_diagnostics ++;
    const int mybuffer_N = myN_diagnostics;
    const int buffer_N   =   N_diagnostics;
    double myData[mybuffer_N][num_var];
    double Data[buffer_N][num_var]; 
    myN_diagnostics = 0;
    int *recvcounts  = new int[numprocs];
    int *displs      = new int[numprocs];

    for(int p=0;p<N_p;p++)
	{
		E		=	0.5*species[0].mass*pow(marker[p].v_par,2.0) + 0.5*species[0].mass*pow(marker[p].v_per,2.0);
		v		=	sqrt(2*E/species[0].mass);
        Lambda  =   marker[p].mu*B0/E;
        if(abs((E-E_0)/E) < windows_width)
        {
            int sgn_v_parallel = marker[p].v_par>0?1:-1;
		    if((1-Lambda) > 0)
			    bracket_psi = marker[p].P_phi/species[0].charge - species[0].mass/species[0].charge*sgn_v_parallel*v*R0*sqrt(1-Lambda);
		    else
			    bracket_psi = marker[p].P_phi/species[0].charge;

            myData[myN_diagnostics][0] = marker[p].w*marker[p].g;
            myData[myN_diagnostics][1] = Lambda;
            myData[myN_diagnostics][2] = marker[p].P_phi;
            myData[myN_diagnostics][3] = bracket_psi;
            myData[myN_diagnostics][4] = marker[p].X[0];
            myData[myN_diagnostics][5] = marker[p].f_over_g*marker[p].g;
            myN_diagnostics++;
        }
	}

    MPI::COMM_WORLD.Gather(&myN_diagnostics,1,MPI::INT,recvcounts,1,MPI::INT,0);
    for(int i=0;i<numprocs;i++)
        recvcounts[i] = recvcounts[i]*num_var;
    displs[0] = 0;
    for(int i=1;i<numprocs;i++)
        displs[i] = displs[i-1] + recvcounts[i-1];
    
    MPI::COMM_WORLD.Gatherv(myData,myN_diagnostics*num_var,MPI::DOUBLE,Data,recvcounts,displs,MPI::DOUBLE,0);

    if(myid == 0)
    {
	    for(int i=0;i<N_diagnostics;i++)
        {
            for(int j=0;j<num_var;j++)
                deltaf_out_vs_time << Data[i][j] << " ";
            deltaf_out_vs_time << endl;
        }
    }
    delete [] recvcounts;
    delete [] displs;
    deltaf_out_vs_time.close();
}


/********************************************************************************************************************
*                                                                                                                   *
*     diagnostics for deltaf vs. E&Lambda to output movie    (fixed P_phi)                                          *
*                                                                                                                   *
*             NUM_TIME_DIAGNOSTICS define number of time interval region for diagnostics                            *
*                   use MPI_Gatherv function to replace Write_shared function      (2013.11.21)                     *
*                                                                                                                   *
********************************************************************************************************************/
void diagnostics_output_deltaf_vs_E_and_Lambda(int index_of_time_record,double P_phi_0,double windows_width)
{
	double Lambda,v,E,bracket_psi;//Lambda = mu*B_0/E
    string fh;
	stringstream ss,tt;
	ss << index_of_time_record;
    tt << setiosflags(ios::fixed) << setprecision(3) << P_phi_0;
	fh = "delta-f_vs_E&Lambda_vs_Time_" + ss.str() + "_P_phi0=" + tt.str() + ".dat";
	char *filename = (char*)fh.c_str();
	ofstream deltaf_out_vs_time(filename);
    int myN_diagnostics = 0;
    int N_diagnostics;
    const int num_var = 6;


    //first step: confirm total number of particles for diagnostics every process
    for(int p=0;p<N_p;p++)
	{
		E		=	0.5*species[0].mass*pow(marker[p].v_par,2.0) + 0.5*species[0].mass*pow(marker[p].v_per,2.0);
		v		=	sqrt(2*E/species[0].mass);
        Lambda  =   marker[p].mu*B0/E;

        if(abs((marker[p].P_phi-P_phi_0)/marker[p].P_phi) < windows_width)
            myN_diagnostics++;
	}

        
    MPI::COMM_WORLD.Reduce(&myN_diagnostics,&N_diagnostics,1,MPI::INT,MPI_SUM,0);
    MPI::COMM_WORLD.Bcast(&N_diagnostics,1,MPI::INT,0);

    //guarantee myN_diagnostics not be zero!!
    if(myN_diagnostics == 0)
        myN_diagnostics ++;
    const int mybuffer_N = myN_diagnostics;
    const int buffer_N   =   N_diagnostics;
    double myData[mybuffer_N][num_var];
    double Data[buffer_N][num_var]; 
    myN_diagnostics = 0;
    int *recvcounts  = new int[numprocs];
    int *displs      = new int[numprocs];

    for(int p=0;p<N_p;p++)
	{
		E		=	0.5*species[0].mass*pow(marker[p].v_par,2.0) + 0.5*species[0].mass*pow(marker[p].v_per,2.0);
		v		=	sqrt(2*E/species[0].mass);
        Lambda  =   marker[p].mu*B0/E;

        if(abs((marker[p].P_phi-P_phi_0)/marker[p].P_phi) < windows_width)
        {
            int sgn_v_parallel = marker[p].v_par>0?1:-1;
		    if((1-Lambda) > 0)
			    bracket_psi = marker[p].P_phi/species[0].charge - species[0].mass/species[0].charge*sgn_v_parallel*v*R0*sqrt(1-Lambda);
		    else
			    bracket_psi = marker[p].P_phi/species[0].charge;

            myData[myN_diagnostics][0] = marker[p].w*marker[p].g;
            myData[myN_diagnostics][1] = Lambda;
            myData[myN_diagnostics][2] = E;
            myData[myN_diagnostics][3] = bracket_psi;
            myData[myN_diagnostics][4] = marker[p].X[0];
            myData[myN_diagnostics][5] = marker[p].f_over_g*marker[p].g;
            myN_diagnostics++;
        }
	}

    MPI::COMM_WORLD.Gather(&myN_diagnostics,1,MPI::INT,recvcounts,1,MPI::INT,0);
    for(int i=0;i<numprocs;i++)
        recvcounts[i] = recvcounts[i]*num_var;
    displs[0] = 0;
    for(int i=1;i<numprocs;i++)
        displs[i] = displs[i-1] + recvcounts[i-1];
    
    MPI::COMM_WORLD.Gatherv(myData,myN_diagnostics*num_var,MPI::DOUBLE,Data,recvcounts,displs,MPI::DOUBLE,0);

    if(myid == 0)
    {
	    for(int i=0;i<N_diagnostics;i++)
        {
            for(int j=0;j<num_var;j++)
                deltaf_out_vs_time << Data[i][j] << " ";
            deltaf_out_vs_time << endl;
        }
    }
    delete [] recvcounts;
    delete [] displs;
    deltaf_out_vs_time.close();
}


/********************************************************************************************************************
*                                                                                                                   *
*                          output particle region in Lambda-<psi> 2d space                                          *
*                                                                                                                   *
*                   use MPI_Gatherv function to replace Write_shared function      (2013.12.06)                     *
*                                                                                                                   *
********************************************************************************************************************/
void output_particle_region()
{
	double Lambda,v,E,bracket_psi;//Lambda = mu*B_0/E
    string fh;
	fh = "particle_region.dat";
	char *filename = (char*)fh.c_str();
	ofstream particle_region(filename);
    int myN_diagnostics = N_p_initial;
    int N_diagnostics;
    const int num_var = 5;


    //first step: confirm total number of particles for diagnostics every process
    MPI::COMM_WORLD.Reduce(&myN_diagnostics,&N_diagnostics,1,MPI::INT,MPI_SUM,0);
    MPI::COMM_WORLD.Bcast(&N_diagnostics,1,MPI::INT,0);

    const int mybuffer_N = myN_diagnostics;
    const int buffer_N   =   N_diagnostics;
    double myData[mybuffer_N][num_var];
    double Data[buffer_N][num_var]; 
    myN_diagnostics = 0;
    int *recvcounts  = new int[numprocs];
    int *displs      = new int[numprocs];

    for(int p=0;p<N_p;p++)
	{
		E		      =   0.5*species[0].mass*pow(marker[p].v_par,2.0) + 0.5*species[0].mass*pow(marker[p].v_per,2.0);
		v		      =   sqrt(2*E/species[0].mass);
        Lambda        =   marker[p].mu*B0/E;

        int sgn_v_parallel = marker[p].v_par>0?1:-1;
		if((1-Lambda) > 0)
			bracket_psi = marker[p].P_phi/species[0].charge - species[0].mass/species[0].charge*sgn_v_parallel*v*R0*sqrt(1-Lambda);
		else
			bracket_psi = marker[p].P_phi/species[0].charge;

        myData[myN_diagnostics][0] = bracket_psi;
        myData[myN_diagnostics][1] = Lambda;
        myData[myN_diagnostics][2] = marker[p].X[0];
        myData[myN_diagnostics][3] = marker[p].P_phi;
        myData[myN_diagnostics][4] = E;
        myN_diagnostics++;

	}

    MPI::COMM_WORLD.Gather(&myN_diagnostics,1,MPI::INT,recvcounts,1,MPI::INT,0);
    for(int i=0;i<numprocs;i++)
        recvcounts[i] = recvcounts[i]*num_var;
    displs[0] = 0;
    for(int i=1;i<numprocs;i++)
        displs[i] = displs[i-1] + recvcounts[i-1];
    
    MPI::COMM_WORLD.Gatherv(myData,myN_diagnostics*num_var,MPI::DOUBLE,Data,recvcounts,displs,MPI::DOUBLE,0);

    if(myid == 0)
    {
	    for(int i=0;i<N_diagnostics;i++)
        {
            for(int j=0;j<num_var;j++)
                particle_region << Data[i][j] << " ";
            particle_region << endl;
        }
    }
    delete [] recvcounts;
    delete [] displs;
    particle_region.close();
}


//add by 2014-07-29
/********************************************************************************************************************
*                                                                                                                   *
*             diagnostics for deltaf vs. E&P_phi to output movie    (fixed mu)                                      *
*                                                                                                                   *
*             NUM_TIME_DIAGNOSTICS define number of time interval region for diagnostics                            *
*                   use MPI_Gatherv function to replace Write_shared function                                       *
*                                                                                                                   *
********************************************************************************************************************/
void diagnostics_movie_use_MPI_Gather_fixed_mu(int index_of_time_record, double mu_0)
{
	double Lambda, v, E, E_0, bracket_psi;
	double windows_width = 5e-2;
	string fh;
	stringstream ss, tt;
	ss << index_of_time_record;
	tt << setiosflags(ios::fixed) << setprecision(2) << mu_0/v_A/v_A*B0;
	fh = "delta-f_vs_E&P_phi_vs_Time_" + ss.str() + "_mu=" + tt.str() + "v_A_2.dat";
	char *filename = (char*)fh.c_str();
	ofstream deltaf_out_vs_time(filename);
	int myN_diagnostics = 0;
	int N_diagnostics;
	const int num_var = 6;


	//first step: confirm total number of particles for diagnostics every process
	for (int p = 0; p<N_p; p++)
	{
		if (abs((marker[p].mu - mu_0)/mu_0) < windows_width)
			myN_diagnostics++;
	}


	MPI::COMM_WORLD.Reduce(&myN_diagnostics, &N_diagnostics, 1, MPI::INT, MPI_SUM, 0);
	MPI::COMM_WORLD.Bcast(&N_diagnostics, 1, MPI::INT, 0);

	//guarantee myN_diagnostics not be zero!!
	if (myN_diagnostics == 0)
		myN_diagnostics++;
	const int mybuffer_N = myN_diagnostics;
	const int buffer_N = N_diagnostics;
	double myData[mybuffer_N][num_var];
	double Data[buffer_N][num_var];
	myN_diagnostics = 0;
	int *recvcounts = new int[numprocs];
	int *displs = new int[numprocs];

	for (int p = 0; p<N_p; p++)
	{
		E = 0.5*species[0].mass*pow(marker[p].v_par, 2.0) + 0.5*species[0].mass*pow(marker[p].v_per, 2.0);
		v = sqrt(2 * E / species[0].mass);
		Lambda = marker[p].mu*B0 / E;

		if (abs((marker[p].mu - mu_0) / mu_0) < windows_width)
		{
			int sgn_v_parallel = marker[p].v_par>0 ? 1 : -1;
			if ((1 - Lambda) > 0)
				bracket_psi = marker[p].P_phi / species[0].charge - species[0].mass / species[0].charge*sgn_v_parallel*v*R0*sqrt(1 - Lambda);
			else
				bracket_psi = marker[p].P_phi / species[0].charge;

			myData[myN_diagnostics][0] = marker[p].w*marker[p].g;
			myData[myN_diagnostics][1] = E;
			myData[myN_diagnostics][2] = marker[p].P_phi;
			myData[myN_diagnostics][3] = bracket_psi;
			myData[myN_diagnostics][4] = marker[p].X[0];
			myData[myN_diagnostics][5] = marker[p].f_over_g*marker[p].g;
			myN_diagnostics++;
		}
	}

	MPI::COMM_WORLD.Gather(&myN_diagnostics, 1, MPI::INT, recvcounts, 1, MPI::INT, 0);
	for (int i = 0; i<numprocs; i++)
		recvcounts[i] = recvcounts[i] * num_var;
	displs[0] = 0;
	for (int i = 1; i<numprocs; i++)
		displs[i] = displs[i - 1] + recvcounts[i - 1];

	MPI::COMM_WORLD.Gatherv(myData, myN_diagnostics*num_var, MPI::DOUBLE, Data, recvcounts, displs, MPI::DOUBLE, 0);

	if (myid == 0)
	{
		for (int i = 0; i<N_diagnostics; i++)
		{
			for (int j = 0; j<num_var; j++)
				deltaf_out_vs_time << Data[i][j] << " ";
			deltaf_out_vs_time << endl;
		}
	}
	delete[] recvcounts;
	delete[] displs;
	deltaf_out_vs_time.close();
}

/************************************************************************************************
*                                                                                               *
*              Initial condition for test particle with fixed mu(MPI version)                   *
*                    only use for TEK(Time Evolution of KAM surfaces) case                      *
*                                                                                               *
************************************************************************************************/
void initia_single_particle_MPI_fix_mu(double E_0, double P_phi_0, double mu_0, int sigma)
{
	double K0;
	double B_value_par, B_par[3];
	double weight_line_auxiliary[2][2];
	double weight_square[2][2];
	double phi;
	double e, v_par, m;
	int i, j, k;
	double psi_par, q, r2;
	double R, Z;
	double W_per, W_par;

	if(myid == 0)
		K0 = P_phi_0 - n[0]/omega_A[0]*E_0;

	MPI::COMM_WORLD.Bcast(&K0, 1, MPI::DOUBLE, 0);


	srand((unsigned)time(NULL)*(myid + 1.0));
	for (int p = 0; p<N_p_diagnosis; p++)
	{
		//define position
		//Random in R and phi at midplane(Z=0)

		marker[p].X[0] = SINGLE_R - DELTA_R + 2 * DELTA_R*(rand()*1.0/RAND_MAX);
		marker[p].X[1] = 0.00;
		marker[p].X[2] = 2*PI*(rand()*1.0/RAND_MAX);
		marker[p].id   = p + myid*N_p_diagnosis;

		i = (int)((marker[p].X[0] - R0 + a*a_b)/dR);
		j = (int)((marker[p].X[1] + a*a_b)/dZ);


		R = R0 - a*a_b + i*dR;
		Z = -a*a_b + j*dZ;

		weight_line_auxiliary[0][0] = (marker[p].X[0] - R)/dR;
		weight_line_auxiliary[1][0] = (marker[p].X[1] - Z)/dZ;

		weight_line_auxiliary[0][1] = 1.0 - weight_line_auxiliary[0][0];
		weight_line_auxiliary[1][1] = 1.0 - weight_line_auxiliary[1][0];


		for (int ii = 0; ii<2; ii++)
			for (int jj = 0; jj<2; jj++)
				weight_square[ii][jj] = weight_line_auxiliary[0][1-ii] * weight_line_auxiliary[1][1-jj];

		psi_par = 0.0;
		B_value_par = 0.0;

		for (int ii=0;ii<2;ii++)
			for (int jj=0;jj<2;jj++)
			{
				B_value_par += B_value[i + ii][j + jj]*weight_square[ii][jj];
				psi_par += psi[i + ii][j + jj]*weight_square[ii][jj];
			}

		if (NUMERICAL_EQUILIBRIUM)
		{
			double g_eq_par = 0.0;
			for (int ii = 0; ii<2; ii++)
				for (int jj = 0; jj<2; jj++)
					g_eq_par += g_eq[i + ii][j + jj]*weight_square[ii][jj];

			B_par[2] = g_eq_par/marker[p].X[0];
		}
		else
			B_par[2] = B0*R0/marker[p].X[0];


		R = marker[p].X[0];
		Z = marker[p].X[1];

		double coeff_a, coeff_b, coeff_c;

		coeff_a = -n[0]/2.0/omega_A[0]*species[0].mass;
		coeff_b = species[0].mass*marker[p].X[0]*B_par[2]/B_value_par;
		coeff_c = species[0].charge*psi_par/alpha - n[0]/omega_A[0]*mu_0*B_value_par - K0;
		if(sigma >0)
			marker[p].v_par = (-coeff_b + sqrt(coeff_b*coeff_b - 4*coeff_a*coeff_c))/2/coeff_a;
		else
			marker[p].v_par = (-coeff_b - sqrt(coeff_b*coeff_b - 4*coeff_a*coeff_c))/2/coeff_a;

		marker[p].mu    = mu_0;
		marker[p].v_per = sqrt(2*mu_0*B_value_par/species[0].mass);

		R = marker[p].X[0];
		e = species[0].charge;
		v_par = marker[p].v_par;
		m = species[0].mass;
		W_par = 0.5*species[0].mass*pow(marker[p].v_par, 2.0);
		W_per = 0.5*species[0].mass*pow(marker[p].v_per, 2.0);

		double K00 = e*psi_par/alpha + v_par*coeff_b - n[0]/omega_A[0]*(W_per + W_par);


//		cout << "marker[p].v_par  : " << marker[p].v_par/v_A  << "v_A  marker[p].v_per " << marker[p].v_per/v_A << "v_A marker[p].X[0]  : " << marker[p].X[0] << " mu :" << marker[p].mu  << " P_phi :" << species[0].mass*marker[p].v_par*marker[p].X[0]*B_par[2]/B_value_par + species[0].charge*psi_par/alpha << " K00 : " << K00 << " id = " << marker[p].id <<" myid = " << myid << " K0 : " << K0 << endl;

		res_index[p] = marker[p].id;

	}
}
#endif
