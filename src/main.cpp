#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include "parameter.h"
#include "tokamak.h"
#include "test1.h"

/*
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
*/

#include <mpi.h>

using namespace std;



int main(int argc,char *argv[])
{
    if(!ISOTROPY)
        v_per_max = 0.0;

	for(int s=0;s<NUM_MODE;s++)
    {
        Xi_initia_cos[s].initia(mR,mZ);
        Xi_initia_sin[s].initia(mR,mZ);
        deltaA_par_sin[s].initia(mR,mZ);
        deltaA_par_cos[s].initia(mR,mZ);
        deltaPhi_sin[s].initia(mR,mZ);
        deltaPhi_cos[s].initia(mR,mZ);
        grad_deltaA_par_sin[s].initia(mR,mZ);
        grad_deltaA_par_cos[s].initia(mR,mZ);
        grad_deltaPhi_sin[s].initia(mR,mZ);
        grad_deltaPhi_cos[s].initia(mR,mZ);
        deltaPhi_sin_3D[s].initia(mR,mZ,mymphi);
        deltaPhi_cos_3D[s].initia(mR,mZ,mymphi);
        deltaE_sin_3D[s].initia(mR,mZ,mymphi);
        deltaE_cos_3D[s].initia(mR,mZ,mymphi);

        //out_put deltaB_R_max
        deltaA_par_cos_3D[s].initia(mR,mZ,mymphi);
        deltaB_initial_3D[s].initia(mR,mZ,mymphi);
        switch(s){
		case 0:
			{
//                X[s]            = B0*0.0005112;
//                X[s]            = B0*0.002223;
//                X[s]            = B0*0.001;
                X[s]            = B0*0.00001;
                Y[s]            = B0*0.00000;
               input_list[s]   = "outzpsi/n=1_Jacobian=R^2/outzpsi";
//              input_list[s]   = "outzpsi/n=2_Jacobian=R^2/outzpsi";
//               input_list[s]   = "outzpsi/n=2_Jacobian=R^2_artificial_frequency/outzpsi";
			}break;
		case 1:
			{
//                X[s]            = B0*0.0003800;
                X[s]            = B0*0.00001;
                Y[s]            = B0*0.00000;
                input_list[s]   = "outzpsi/n=3_Jacobian=R^2/outzpsi";
//                input_list[s]   = "outzpsi/n=3_Jacobian=R^2_artificial_frequency/outzpsi";
			}break;
        } 
    }


	MPI::Init(argc,argv);
	MPI::Status status;
	myid = MPI::COMM_WORLD.Get_rank();
	numprocs =  MPI::COMM_WORLD.Get_size();

    dphi= 2*PI/(mymphi-1);
	
    phi_left  = 0.0;
	phi_right = 2*PI;
						
	clean_create_files();
	r2R_accuracy();
    equilibrium_field_2d();
	initia_output();
	initia_TAE_NOVA();
	initia_parameter_output();
	calc_wave_energy();
	loading_particle();

    //output Time Evolution of KAM surfaces
    if(TEK)
    {
        NUM_LINE_AMP = load_amp_file();
        RES = true;
    }

    if(COF || RES)
    {
        AMP         = false;
        DIAGNOSTICS = false;
        ADIABATIC   = 0;
        FULL_F      = true;
        if(COF)
        {
            //for CPN = true, nonlinear term should be retained to calc orbit frequency and Amplitude should not be zero!
            if(!CPN)
                for(int s=0;s<NUM_MODE;s++)
                {
                    X[s]  = 0.0;
                    Y[s]  = 0.0;
                }
            else
            {
                NUM_LINE_AMP = load_amp_file();
                index_t_start_CPN = (int)(t_start_CPN*sqrt(alpha1)/(t_input[1]-t_input[0]));
            }
        }
    }

    if(RES)
    {
        if(TEK)
        {
            if (TEK_mu)
                initia_single_particle_MPI_fix_mu(E_0_TEK,P_phi_0_TEK,mu_0_TEK,sigma_TEK);
            else
                initia_single_particle_MPI();
        }
        else
            initia_single_particle();
    }

    if(ONS)
        diagnostics_on_number_of_particle();

    //output Particle_Region in XX-Lambda space
    if(PRL)
        output_particle_region();

	if(myid==0)
		start = (clock_t)MPI_Wtime();
    if(COF && !RES)
    {
        //use small time step for large velocity
        dt = dt/2.0;
        if(OPS == 1)//Oribit in Phase Space (E,P_phi)
            for(int i=0;i<NUM_E;i++)
                for(int j=0;j<NUM_PPHI;j++)
                {	
                    double B_value_par;
                    double E;
                    particle_vector<double> X;
                    OUTPUT_TRACK = false;

                    if(j%(NUM_PPHI/NUM_TRACK_PER_PPHI) == 0 && i%(NUM_E/NUM_TRACK_PER_E) == 0)
                        OUTPUT_TRACK = true;
                
                    SINGLE_R     = R0 + 0.01*a + 0.95*a*j/(NUM_PPHI-1);
                    SINGLE_Z     = 0.00;
                    SINGLE_phi   = 0.00;
                
                    X[0] = SINGLE_R;
                    X[1] = SINGLE_Z;
                    X[2] = SINGLE_phi;

                    B_value_par = get_B_value(X);
 
                    E             =   0.1*0.5*species[0].mass*v_A*v_A + 2.9*0.5*species[0].mass*v_A*v_A*i/(NUM_E-1);
                    SINGLE_mu     =   E*LAMBDA/B0;
                    SINGLE_v_per  =   sqrt(2*SINGLE_mu*B_value_par/species[0].mass);
                    SINGLE_Lambda =   LAMBDA;

                    //true for co-passing while false for counter-passing particles
                    if(CO_PASSING)
                        SINGLE_v_par = + sqrt(2*E/species[0].mass - pow(SINGLE_v_per,2.0));
                    else
                        SINGLE_v_par = - sqrt(2*E/species[0].mass - pow(SINGLE_v_per,2.0));

		            marker[0].X[0]	= SINGLE_R;
		            marker[0].X[1]	= SINGLE_Z;
		            marker[0].X[2]	= SINGLE_phi;
		            marker[0].v_par	= SINGLE_v_par;
		            marker[0].v_per	= SINGLE_v_per;
                    marker[0].mu    = SINGLE_mu;
		            marker[0].id	= 0;
                    
                    LOSS = false;

                    t = 0;
                    stepon_RK4();
                    if(myid == 0)
                        cout << "NUM OF PARTICLES : " << i << " " << j << endl;
                }
        else if(OPS == 2)//Oribit in Phase Space (Lambda,E)
           for(int i=0;i<NUM_LAMBDA;i++)
                for(int j=0;j<NUM_E;j++)
                {	
                    double E_OPS,LAMBDA_OPS;
                    OUTPUT_TRACK = false;

                    if(j%(NUM_E/NUM_TRACK_PER_E) == 0 && i%(NUM_LAMBDA/NUM_TRACK_PER_LAMBDA) == 0)
                        OUTPUT_TRACK = true;


                    //confirm v_|| of particle
                    E_OPS         =   0.1*0.5*species[0].mass*v_A*v_A + 1.0*0.5*species[0].mass*v_A*v_A*j/(NUM_E-1);
                    LAMBDA_OPS    =   LAMBDA_min + (LAMBDA_max - LAMBDA_min)*i/(NUM_LAMBDA-1);
                    SINGLE_mu     =   E_OPS*LAMBDA_OPS;
                    
                    
                    double P_phi_OPS_last = 1e10;
                    int ii;
                    //confirm location of particle in R
                    for(ii=mR/2;ii<mR;ii++)
                    {
                        double P_phi_OPS,R_OPS,psi_OPS,B_phi_OPS,B_value_OPS,SINGLE_v_par_OPS;
                        R_OPS             = R0-a*a_b+ii*dR;
                        psi_OPS           = psi[ii][mZ/2];
                        B_phi_OPS         = B[2][ii][mZ/2];
                        B_value_OPS       = B_value[ii][mZ/2];

                        if(2/species[0].mass*(E_OPS - SINGLE_mu*B_value_OPS) < 0)
                            continue;
                        else
                        {
                            if(CO_PASSING)
                                SINGLE_v_par_OPS  = + sqrt(2/species[0].mass*(E_OPS - SINGLE_mu*B_value_OPS));
                            else
                                SINGLE_v_par_OPS  = - sqrt(2/species[0].mass*(E_OPS - SINGLE_mu*B_value_OPS));


                            P_phi_OPS = species[0].mass*SINGLE_v_par_OPS*R_OPS*(B_phi_OPS/B_value_OPS)+species[0].charge*psi_OPS/alpha;
                        
                            if( ((P_phi_OPS > P_phi_0 && P_phi_OPS_last < P_phi_0) || (P_phi_OPS < P_phi_0 && P_phi_OPS_last > P_phi_0)) && abs(P_phi_OPS - P_phi_OPS_last) < 1e3 )
                            {
                                SINGLE_R     = R_OPS;
                                SINGLE_Z     = 0.00;
                                SINGLE_phi   = 0.00;
                                SINGLE_v_par = SINGLE_v_par_OPS;
                                SINGLE_v_per = sqrt(2*SINGLE_mu*B_value_OPS/species[0].mass);
                                SINGLE_Lambda= LAMBDA_OPS;
                                LOSS = false;
                                break;
                            }
                            P_phi_OPS_last = P_phi_OPS;
                        }
                    }

                    if(ii == mR)
                    {
                        cout << "particle is out of the simulation region!!!!!" << endl;
                        LOSS = true;
                    }

                    

		            marker[0].X[0]	= SINGLE_R;
		            marker[0].X[1]	= SINGLE_Z;
		            marker[0].X[2]	= SINGLE_phi;
		            marker[0].v_par	= SINGLE_v_par;
		            marker[0].v_per	= SINGLE_v_per;
                    marker[0].mu    = SINGLE_mu;
		            marker[0].id	= 0;

                    

                    t = 0;
                    stepon_RK4();
                    if(myid == 0)
                        cout << "NUM OF PARTICLES : " << i << " " << j << endl;
                }
    }
    else
        stepon_RK4();

	if(myid==0)
		finish = (clock_t)MPI_Wtime();

	if(myid==0)
	{	
		totaltime=(double)(finish-start);
		cout << setiosflags(ios::fixed) << setprecision(8) << "totaltime : " << totaltime << endl;
		cout << "complete!  &&  t = " << t << endl;
		cout << "r,theta,phi " << marker[0].X[0] << " " << marker[0].X[1] << "  " << marker[0].X[2] <<  endl;
	}
	MPI::Finalize();
	return 0;

}

