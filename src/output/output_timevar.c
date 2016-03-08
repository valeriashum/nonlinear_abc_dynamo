#include <stdlib.h>
#include <string.h>
#include <complex.h> 
#include "../common.h"
#include "../gfft.h"
#include "../shear.h"
#include "../debug.h"
#include "../particles.h"
#include "../isotropyT.h"
#include "../snoopy.h"
#define MAX_N_BIN					10000
#define OUTPUT_SPECTRUM_K_BIN		(2.0 * M_PI/100.0)

/***********************************************************/
/** 
	find the maximum of a real array of size 2*NTOTAL_COMPLEX
	@param wri array in which we want to know the maximum
*/
/***********************************************************/

double find_max(double *wri) {
	double q0;
	int i,j,k,idx;
	q0=wri[0];
	
	for(i=0 ; i < NX/NPROC ; i++) {
		for(j=0 ; j < NY ; j++) {
			for(k=0 ; k < NZ ; k++) {
#ifdef WITH_2D
				idx = j + i * (NY+2);
#else
				idx = k + j * (NZ + 2) + i * NY * (NZ + 2);
#endif
				if(q0<wri[idx]) q0=wri[idx];
			}
		}
	}
	
	return(q0);
}

/***********************************************************/
/** 
	find the minimum of a real array of size 2*NTOTAL_COMPLEX
	@param wri array in which we want to know the minimum
*/
/***********************************************************/

double find_min(double *wri) {
	double q0;
	int i,j,k,idx;
	q0=wri[0];
	
	for(i=0 ; i < NX/NPROC ; i++) {
		for(j=0 ; j < NY ; j++) {
			for(k=0 ; k < NZ ; k++) {
#ifdef WITH_2D
				idx = j + i * (NY+2);
#else
				idx = k + j * (NZ + 2) + i * NY * (NZ + 2);
#endif
				if(q0>wri[idx]) q0=wri[idx];
			}
		}
	}
	
	return(q0);
}

/***********************************************************/
/** 
	compute the correlation between 2 fields
*/
/***********************************************************/
double compute_2correlation(double *wri1, double *wri2) {
	double q0;
	int i,j,k;
	q0=0.0;
	
	for(i = 0 ; i < 2 * NTOTAL_COMPLEX ; i ++) {
		q0 += wri1[i] * wri2[i] / ((double) NTOTAL);
	}
	
	return(q0);
}

/***********************************************************/
/** 
	compute the correlation between 3 fields
*/
/***********************************************************/
double compute_3correlation(double *wri1, double *wri2, double *wri3) {
	double q0;
	int i,j,k;
	q0=0.0;
	
	for(i = 0 ; i < 2 * NTOTAL_COMPLEX ; i ++) {
		q0 += wri1[i] * wri2[i] * wri3[i] / ((double) NTOTAL);
	}
	
	return(q0);
}

/***********************************************************/
/** 
	Write statistical quantities using text format in the file
	timevar. 
	@param fldi Field structure from which the statistical quantities are derived.
	@param t Current time of the simulation
*/
/***********************************************************/

void output_timevar(const struct Field fldi,
					const double t) {
					
	FILE *ht;
	double output_var;
	static int warning_flag = 0;
    double  threshold, em_curr;
    int i,j,k, ii;
    complex double  B1,B2,B3,B4,B5;
    complex double  r1,r2,r3,r4,r5;
    int             i1,i2,i3,i4,i5;
    int             rank1,rank2,rank3,rank4,rank5;

	DEBUG_START_FUNC;
	

	// Open the timevar file
	if(rank==0) {
		ht=fopen("timevar","a");
	}
    
// Transform a few basic variables.
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = fldi.vx[i];
		w2[i] = fldi.vy[i];
		w3[i] = fldi.vz[i];
        w4[i] = pow(pow(fabs(fldi.vx[i]),2) + pow(fabs(fldi.vy[i]),2) + pow(fabs(fldi.vz[i]),2) , 0.5);

#ifdef BOUSSINESQ
		w4[i] = fldi.th[i];
#endif
#ifdef COMPRESSIBLE
		w4[i] = fldi.d[i];
#endif
#ifdef MHD
		w5[i] = fldi.bx[i];
		w6[i] = fldi.by[i];
		w7[i] = fldi.bz[i];
        w9[i] = pow(pow(fabs(fldi.bx[i]),2) + pow(fabs(fldi.by[i]),2) + pow(fabs(fldi.bz[i]),2) , 0.5);
#endif
	}
    if(param.output_magnetic_field && t>0. ) {
		// Compute the magnetic field B = sqrt BxBx* + ByBy* + BzBz*
       
    B5 = 0.0 + I*0.0;   B4 = 0.0 + I*0.0;   B3 = 0.0 + I*0.0;   B2 = 0.0 + I*0.0;   B1 = 0.0 + I*0.0;
    r5 = 0.0 + I*0.0;   r4 = 0.0 + I*0.0;   r3 = 0.0 + I*0.0;   r2 = 0.0 + I*0.0;   r1 = 0.0 + I*0.0;
    i5 = 0;             i4 = 0;             i3 = 0;             i2 = 0;             i1 = 0;
    rank5=0;            rank4=0;            rank3=0;            rank2=0;            rank1=0;  
    

    for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
			    w16[IDX3D] = pow( pow(fabs(fldi.bx[IDX3D]),2) + pow(fabs(fldi.by[IDX3D]),2) + pow(fabs(fldi.bz[IDX3D]),2)  , 0.5); 
                
                //if(isfinite(wr17[IDX3D])==1){
                  
//                }

                // Find 5 top dominant wavenumbers: 
                //
                if ( creal(B5) < creal(w16[IDX3D])){
                    B1 = B2;         r1 = r2;            i1 = i2;        rank1 = rank2;
                    B2 = B3;         r2 = r3;            i2 = i3;        rank2 = rank3;
                    B3 = B4;         r3 = r4;            i3 = i4;        rank3 = rank4;
                    B4 = B5;         r4 = r5;            i4 = i5;        rank4 = rank5;
                    B5 = w16[IDX3D]; r5 =w17[IDX3D];     i5 = IDX3D;     rank5 = rank;
                }
                else{ 
                    if (creal(B4) < creal(w16[IDX3D])) {
                        B1 = B2;           r1 = r2;             i1 = i2;        rank1 = rank2;
                        B2 = B3;           r2 = r3;             i2 = i3;        rank2 = rank3;
                        B3 = B4;           r3 = r4;             i3 = i4;        rank3 = rank4;
                        B4 = w16[IDX3D];   r4 =w17[IDX3D];      i4 = IDX3D;     rank4 = rank;
                    }
                    else{
                        if( creal(B3) < creal(w16[IDX3D])) {
                            B1 = B2;            r1 = r2;             i1 = i2;        rank1 = rank2;
                            B2 = B3;            r2 = r3;             i2 = i3;        rank2 = rank3;
                            B3 = w16[IDX3D];    r3 = w17[IDX3D];     i3 = IDX3D;     rank3 = rank;
                        } 
                        else{
                            if( creal(B2) < creal(w16[IDX3D])) {
                                B1 = B2;           r1 = r2;             i1 = i2;        rank1 = rank2;  
                                B2 = w16[IDX3D];   r2 = w17[IDX3D];     i2 = IDX3D;     rank2 = rank;
                            } 
                            else{
                                if( creal(B1) < creal(w16[IDX3D])) {
                                   B1 = w16[IDX3D];  r1 = w17[IDX3D]; i1 = IDX3D;  rank1=rank;                           
                                }                         
                            } 
                        } 
                    }   
                }
             }
          //  }
		}
    }
 

    if (rank==rank1) printf("______________________________________________________________________\n");
    if (rank==rank1) printf("r5= %e, B-5 = %e, k=(%f,%f,%f) in rank %d\n",creal(r1), creal(B1), kx[i1],ky[i1],kz[i1],rank);
    if (rank==rank2) printf("r4= %e, B-4 = %e, k=(%f,%f,%f) in rank %d\n",creal(r2), creal(B2), kx[i2],ky[i2],kz[i2],rank);
    if (rank==rank3) printf("r3= %e, B-3 = %e, k=(%f,%f,%f) in rank %d\n",creal(r3), creal(B3), kx[i3],ky[i3],kz[i3],rank);
    if (rank==rank4) printf("r2= %e, B-2 = %e, k=(%f,%f,%f) in rank %d\n",creal(r4), creal(B4), kx[i4],ky[i4],kz[i4],rank);
    if (rank==rank5) printf("r1= %e, B-1 = %e, k=(%f,%f,%f) in rank %d\n",creal(r5), creal(B5), kx[i5],ky[i5],kz[i5],rank);
    if (rank==rank5) printf("______________________________________________________________________\n");

// Normalise by dividing by B_max - range < 1 
       /* if(rank == 0){
        // Use GNUPLOT to plot B vs IDX3D
            FILE *gnuplot_pipe = popen("gnuplot -persist", "w");                      // GNU
            fprintf(gnuplot_pipe, "reset\n");                                         // GNU
            fprintf(gnuplot_pipe, "set term eps\n"); //gnuplot commands
            fprintf(gnuplot_pipe, "set output \"plots/test_%f.eps\"\n", t); //gnuplot commands
            fprintf(gnuplot_pipe, "set title \"<B> vs IDX3D at t=%f\"\n", t);         // GNU
            //fprintf(gnuplot_pipe, "set yrange [0.5:1]\n");                            // GNU
            fprintf(gnuplot_pipe, "set xlabel \"IDX3D = k+jNz+iNzNy \"\n");           // GNU
            fprintf(gnuplot_pipe, "set ylabel \"<B>/<B>max in Fourier Space \"\n");          // GNU
             fprintf(gnuplot_pipe, "set nokey\n");                                    // GNU
            fprintf(gnuplot_pipe, "plot '-' with lines\n");                           // GNU
            for( i = 0; i < NX_COMPLEX/NPROC; i++) {                                  // GNU
                for( j = 0; j < NY_COMPLEX; j++) {                                    // GNU
	        	    for( k = 0; k < NZ_COMPLEX; k++) {                                // GNU
                        fprintf(gnuplot_pipe , "%i %e\n", IDX3D, creal(w7[IDX3D])/creal(B5));   // GNU
                    }
                }
            }

            fprintf(gnuplot_pipe , "e\n");                                            // GNU
            fflush(gnuplot_pipe);                                                     // GNU
        }*/

 }
	/***********************************************************/
	/** 
		 Magnitude of the overlap integral of u with the ABC flow for A=B=C=1
		 @param wri array in which we want to know the maximum

	*/
	/***********************************************************/

	// This variable requires the velocity field in wavenumber space
	// In retrospect I don't think there's any reason that this code needs to
	// be where it is - yes w1, w2 and w3 are transformed to positions space
	// below, but I could still access fldi.vx etc.!
	double ABCmag;
    	int idx3d;
	for( i = 0 ; i < param.timevar_vars.length ; i++ ) {
	   // Code added by Harry Braviner 2015-03-17
	   if(!strcmp(param.timevar_vars.name[i], "ABCmag")){
		  // Magnitude of the overlap integral of u with the ABC flow for A=B=C=1
          ABCmag = 0.0;
		  if (rank == 0){
	        // Add (1/V)*integral[ u_x*sin(2*pi*z) dV ]
            idx3d = 1;
            ABCmag += 2.0*fldi.vx[idx3d]*I;
	        // Add (1/V)*integral[ u_y*cos(2*pi*z) dV ]
            idx3d = 1;
            ABCmag += 2.0*fldi.vy[idx3d];
	        // Add (1/V)*integral[ u_x*cos(2*pi*y) dV ]
            idx3d = (1 * NZ_COMPLEX);
            ABCmag += fldi.vx[idx3d];
            idx3d = ((NY_COMPLEX - 1) * NZ_COMPLEX);
            ABCmag += fldi.vx[idx3d];
	        // Add (1/V)*integral[ u_z*sin(2*pi*y) dV ]
            idx3d = (1 * NZ_COMPLEX);
            ABCmag += fldi.vz[idx3d]*I;
            idx3d = ((NY_COMPLEX - 1) * NZ_COMPLEX);
            ABCmag += fldi.vz[idx3d]*(-I);
	        // Add (1/V)*integral[ u_y*sin(2*pi*x) dV ]
            idx3d = (NZ_COMPLEX * NY_COMPLEX * 1);
            ABCmag += fldi.vy[idx3d]*I;
	        // Add (1/V)*integral[ u_z*cos(2*pi*x) dV ]
            idx3d = (NZ_COMPLEX * NY_COMPLEX * 1);
            ABCmag += fldi.vz[idx3d];
		  }
          if (rank == NPROC-1) {
	        // Add (1/V)*integral[ u_y*sin(2*pi*x) dV ]
            idx3d = (NZ_COMPLEX * NY_COMPLEX * (NX_COMPLEX/NPROC - 1));
            ABCmag += fldi.vy[idx3d]*(-I);
	        // Add (1/V)*integral[ u_z*cos(2*pi*x) dV ]
            idx3d = (NZ_COMPLEX * NY_COMPLEX * (NX_COMPLEX/NPROC - 1));
            ABCmag += fldi.vz[idx3d];
          }
          reduce(&ABCmag,1);
          ABCmag = ABCmag / (6.0 * (double) NTOTAL);
	   }
	   // End of code added by Harry Braviner 2015-03-17
	}

	gfft_c2r(w1);
	gfft_c2r(w2);
	gfft_c2r(w3);
#ifdef BOUSSINESQ	// Boussinesq and compressibility are in principle mutually exclusive.
	gfft_c2r(w4);
#endif
#ifdef COMPRESSIBLE
	gfft_c2r(w4);
	check_positivity(wr4);
#endif
#ifdef MHD
	gfft_c2r(w5);
	gfft_c2r(w6);
	gfft_c2r(w7);
#endif

	
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
#ifndef COMPRESSIBLE
		wr1[i] = wr1[i] / ((double) NTOTAL );
		wr2[i] = wr2[i] / ((double) NTOTAL );
		wr3[i] = wr3[i] / ((double) NTOTAL );
#ifdef BOUSSINESQ
		wr4[i] = wr4[i] / ((double) NTOTAL );
#endif
#else
		wr1[i] = wr1[i] / wr4[i];		// Compute the real velocity field
		wr2[i] = wr2[i] / wr4[i];
		wr3[i] = wr3[i] / wr4[i];
		wr4[i] = wr4[i] / ((double) NTOTAL);
#endif
		
#ifdef MHD
		wr5[i] = wr5[i] / ((double) NTOTAL );
		wr6[i] = wr6[i] / ((double) NTOTAL );
		wr7[i] = wr7[i] / ((double) NTOTAL );
		
#endif
	}
/********************************************/
/** LOOP OVER ALL THE REQUESTED VARIABLES ***/
/** do not use i as a dummy variable      ***/
/********************************************/
    for( i = 0 ; i < param.timevar_vars.length ; i++ ) {
	
		if(!strcmp(param.timevar_vars.name[i],"t")) {
			// current time
			output_var=t;
		}
		else if(!strcmp(param.timevar_vars.name[i],"ev")) {
			// kinetic energy
#ifndef COMPRESSIBLE
			output_var = energy(fldi.vx) + energy(fldi.vy)+energy(fldi.vz);
#else
			output_var = 0.5 * ( compute_3correlation(wr1, wr1, wr4) +
								 compute_3correlation(wr2, wr2, wr4) +
								 compute_3correlation(wr3, wr3, wr4) );
#endif
			reduce(&output_var,1);
        }
		else if(!strcmp(param.timevar_vars.name[i],"disp_rate")){
			output_var = dissipation_rate(fldi.vx, nu) +
			             dissipation_rate(fldi.vy, nu) +
									 dissipation_rate(fldi.vz, nu);
			reduce(&output_var, 1);
		}
		else if(!strcmp(param.timevar_vars.name[i],"Rxy")){
			// x,y component of the box-averaged Reynolds stress
			output_var = correlator(fldi.vx,fldi.vy);
			reduce(&output_var, 1);
		}
		else if(!strcmp(param.timevar_vars.name[i],"RottaEtaSquared")){
			// x,y component of the box-averaged Reynolds stress
			output_var = RottaEtaSquared(fldi);
			// Reduce done in RottaEtaSquared()
		}
		else if(!strcmp(param.timevar_vars.name[i],"RottaXiCubed")){
			// x,y component of the box-averaged Reynolds stress
			output_var = RottaXiSquared(fldi);
			// Reduce done in RottaXiSquared()
		}
        

#ifdef OSCILLATORY_SHEAR
		else if(!strcmp(param.timevar_vars.name[i],"intSRxy")){
			// x,y component of the time integral of sin(wt) * Rxy
			output_var = integrated_sin_Rxy;
			reduce(&output_var, 1);
		}
		else if(!strcmp(param.timevar_vars.name[i],"intCRxy")){
			// x,y component of the time integral of cos(wt) * Rxy
			output_var = integrated_cos_Rxy;
			reduce(&output_var, 1);
		}
#endif
		else if(!strcmp(param.timevar_vars.name[i],"vxmax")) {
			// maximum of vx component
			output_var=find_max(wr1);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[i],"vymax")) {
			// maximum of vy component
			output_var=find_max(wr2);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[i],"vzmax")) {
			// maximum of vz component
			output_var=find_max(wr3);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[i],"vxmin")) {
			// minimum of vx component
			output_var=find_min(wr1);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[i],"vymin")) {
			// minimum of vy component
			output_var=find_min(wr2);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[i],"vzmin")) {
			// minimum of vz component
			output_var=find_min(wr3);
			reduce(&output_var,3);
        }

// k_alfvenic
        else if(!strcmp(param.timevar_vars.name[i],"k_alfvenic")) {
			// characteristic perpendicular wavenumber
            // Dastgeer Shaikh 2009 [<|i k x u|^2>/<u^2>]^1/2

		// nabla x u 
            for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
		        w10[j] = I * mask[j] * (ky[j] * w3[j] - kz[j] * w2[j]);
	            w11[j] = I * mask[j] * (kz[j] * w1[j] - kxt[j]* w3[j]);
	            w12[j] = I * mask[j] * (kxt[j]* w2[j] - ky[j] * w1[j]);
    	    }      
        //  [  <|B x (nabla x B)|^2> / <B^4> ]^0.5
            output_var= pow(( energy(w10)+energy(w11)+energy(w12)) / (energy(w1)+energy(w2)+energy(w3) ), 0.5)/(2*M_PI);
		    reduce(&output_var,1);
        }

		else if(!strcmp(param.timevar_vars.name[i],"vxvy")) {
			// incompressible reynolds stress
			output_var=compute_2correlation(wr1,wr2);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[i],"ABCmag")) {
			// overlap integral with the ABC flow
			output_var=ABCmag;
            // Do NOT perform a reduce here - it is done above.
		}
		else if(!strcmp(param.timevar_vars.name[i],"OL2012A")) {
			// Output the quantity A_1jj1 - 2A_1221 - 2A_2121
			// of equation (31) of Ogilvie & Lesur (2012)
			output_var = OL2012A(fldi.vx, fldi.vy);
			reduce(&output_var,1);
		}
#ifdef WITH_OL2012B
		else if(!strcmp(param.timevar_vars.name[i],"OL2012B")) {
			output_var = OL2012B(fldi);
			reduce(&output_var,1);
		}
#endif
		else if(!strcmp(param.timevar_vars.name[i],"hv")) {
			// kinetic helicity
			// Compute vector potential
			for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
				w10[j] = I * ik2t[j] * (ky[j] * fldi.vz[j] - kz[j] * fldi.vy[j] );
				w11[j] = I * ik2t[j] * (kz[j] * fldi.vx[j] - kxt[j]* fldi.vz[j] );
				w12[j] = I * ik2t[j] * (kxt[j]* fldi.vy[j] - ky[j] * fldi.vx[j] );
			}
	
			gfft_c2r(w10);
			gfft_c2r(w11);
			gfft_c2r(w12);
	
			for( j = 0 ; j < 2*NTOTAL_COMPLEX ; j++) {
				wr10[j] = wr10[j] / ((double) NTOTAL );
				wr11[j] = wr11[j] / ((double) NTOTAL );
				wr12[j] = wr12[j] / ((double) NTOTAL );
			}

			output_var=compute_2correlation(wr10,wr1)+compute_2correlation(wr11,wr2)+compute_2correlation(wr12,wr3);
			reduce(&output_var,1);
		}
		
		else if(!strcmp(param.timevar_vars.name[i],"w2")) {
			// total enstrophy
			// Compute vorticity
			for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
				w10[j] = I * (ky[j]  * fldi.vz[j] - kz[j]  * fldi.vy[j]);
				w11[j] = I * (kz[j]  * fldi.vx[j] - kxt[j] * fldi.vz[j]);
				w12[j] = I * (kxt[j] * fldi.vy[j] - ky[j]  * fldi.vx[j]);
			}
			output_var=energy(w10)+energy(w11)+energy(w12);
			reduce(&output_var,1);
		}

		

#ifdef MHD
		else if(!strcmp(param.timevar_vars.name[i],"em")) {
			// magnetic energy
			output_var = energy(fldi.bx) + energy(fldi.by)+energy(fldi.bz);
            reduce(&output_var,1);   
            //MPI_Printf("From o_timevar.c: t=%f, after energy=%e\n", t, output_var );
            em_curr = output_var;
            
            if (isfinite(em_curr)==0){
                    MPI_Printf("Terminated due to error in EM.\n");
                    if(rank==0) {
		            	fprintf(ht,"%08e\n",output_var);
                		fclose(ht);
                    }
                    dump_immediate(fldi, t);
                    MPI_Printf("Dump file saved.\n");
#ifdef MPI_SUPPORT
		            MPI_Finalize();
#endif
		            exit(0);
            }
		}
     
        else if(!strcmp(param.timevar_vars.name[i],"em_gr")) {
                // growth rate of energy 
                output_var= (log (em_curr/em_prev))/(2*param.toutput_time);
                em_prev = em_curr;
                if (isfinite(em_gr)==1 && isfinite(output_var)==1 && fabs((em_gr - output_var)/output_var) < 5e-4 ){
                    param.toutput_time = 2.0;
                    threshold = 1e-6;  
                }
                if (isfinite(em_gr)==1 && isfinite(output_var)==1 && fabs((em_gr - output_var)/output_var) > 5e-4 ){
                    threshold = 1e-4 ; 
                }

                //MPI_Printf("em_gr_prev=%e, em_gr=%e, dif=%e\n", em_gr,output_var,fabs((em_gr - output_var)/output_var)) ;
                // find the saturated growth rate
                if (isfinite(em_gr)==1 && isfinite(output_var)==1 && fabs((em_gr - output_var)/output_var) < threshold ){
                    MPI_Printf("Terminated due to saturated growth rate.Error=%e\n", fabs((em_gr - output_var)/output_var));
                    if(rank==0) {
		            	fprintf(ht,"%08e\n",output_var);
                		fclose(ht);
                    }
                    dump_immediate(fldi, t);
                    MPI_Printf("Dump file saved.\n");
#ifdef MPI_SUPPORT
		            MPI_Finalize();
#endif
		            exit(0);
                }

                // If the steady state is reached, change output_timevar time 
                // to more frequent outputs 
              
                MPI_Printf("t=%f, energy growth rate=%e, perc error=%e\n",t,output_var, fabs((em_gr -output_var)/output_var));	
                em_gr = output_var;
               
        }
        else if(!strcmp(param.timevar_vars.name[i],"Rm")) {
			// current R_m
			output_var=param.reynolds_m;
		}
		else if(!strcmp(param.timevar_vars.name[i],"bxmax")) {
			// maximum of bx component
			output_var=find_max(wr5);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[i],"bymax")) {
			// maximum of by component
			output_var=find_max(wr6);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[i],"bzmax")) {
			// maximum of bz component
			output_var=find_max(wr7);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[i],"bxmin")) {
			// minimum of bx component
			output_var=find_min(wr5);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[i],"bymin")) {
			// minimum of by component
			output_var=find_min(wr6);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[i],"bzmin")) {
			// minimum of bz component
			output_var=find_min(wr7);
			reduce(&output_var,3);
		}

		
// Characteristic Scales 
// k||
	    else if(!strcmp(param.timevar_vars.name[i],"k_par")) {
			// characteristic wavenumber parallel to the magnetic field
            // Shekochinin 2004 [<|B.del B|^2>/<B^4>]^1/2
            
            for( j = 0 ; j < 2*NTOTAL_COMPLEX ; j++) {
		        wr10[j] = wr5[j] * wr5[j] / ((double) NTOTAL*NTOTAL);
	    	    wr11[j] = wr6[j] * wr6[j] / ((double) NTOTAL*NTOTAL);
		        wr12[j] = wr7[j] * wr7[j] / ((double) NTOTAL*NTOTAL);
	        	wr13[j] = wr5[j] * wr6[j] / ((double) NTOTAL*NTOTAL);
	        	wr14[j] = wr5[j] * wr7[j] / ((double) NTOTAL*NTOTAL);
	        	wr15[j] = wr6[j] * wr7[j] / ((double) NTOTAL*NTOTAL);         
        	}

            gfft_r2c_t(wr10);
            gfft_r2c_t(wr11);
            gfft_r2c_t(wr12);
            gfft_r2c_t(wr13);
    	    gfft_r2c_t(wr14);
	        gfft_r2c_t(wr15);
    	             
            for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
		         w8[j] = I * mask[j] * (kxt[j] * w10[j] + ky[j] * w13[j] + kz[j] * w14[j]);
	    	    w13[j] = I * mask[j] * (kxt[j] * w13[j] + ky[j] * w11[j] + kz[j] * w15[j]);
		        w14[j] = I * mask[j] * (kxt[j] * w14[j] + ky[j] * w15[j] + kz[j] * w12[j]);
            }   
		                 //[  <|B.del B|^2> / <B^4> ]^0.5
            output_var= pow(( energy(w8)+energy(w13)+energy(w14)) / (energy(w10) +energy(w11) +energy(w12) ), 0.5)/(2*M_PI);
		    reduce(&output_var,1);
        	
        }
// k_BxJ
        else if(!strcmp(param.timevar_vars.name[i],"k_BxJ")) {
			// characteristic perpendicular wavenumber
            // Shekochinin 2004 [<|B x J|^2>/<B^4>]^1/2

		// nabla x B 
            for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
		        w10[j] = I * mask[j] * (ky[j] * w7[j] - kz[j] * w6[j]);
	            w11[j] = I * mask[j] * (kz[j] * w5[j] - kxt[j]* w7[j]);
	            w12[j] = I * mask[j] * (kxt[j]* w6[j] - ky[j] * w5[j]);
    	    }
        
            gfft_c2r_t(w5);
    	    gfft_c2r_t(w6);
            gfft_c2r_t(w7);
        	gfft_c2r_t(w10);
        	gfft_c2r_t(w11);
            gfft_c2r_t(w12);
     
        // B x (nabla x B)
            for( j = 0 ; j < 2*NTOTAL_COMPLEX ; j++) {
		        wr13[j] = (wr6[j] * wr12[j] - wr7[j] * wr11[j]) / ((double) NTOTAL*NTOTAL);
		        wr14[j] = (wr7[j] * wr10[j] - wr5[j] * wr12[j]) / ((double) NTOTAL*NTOTAL);
		        wr15[j] = (wr5[j] * wr11[j] - wr6[j] * wr10[j]) / ((double) NTOTAL*NTOTAL);
        // B^2
                wr10[j] = wr5[j]* wr5[j]/ ((double) NTOTAL*NTOTAL); 
                wr11[j] = wr6[j]* wr6[j]/ ((double) NTOTAL*NTOTAL);
                wr12[j] = wr7[j]* wr7[j]/ ((double) NTOTAL*NTOTAL); 
	         }
            gfft_r2c_t(wr10);
            gfft_r2c_t(wr11);
            gfft_r2c_t(wr12);
            gfft_r2c_t(wr13);
	        gfft_r2c_t(wr14);
	        gfft_r2c_t(wr15);
    
               
        //  [  <|B x (nabla x B)|^2> / <B^4> ]^0.5
            output_var= pow(( energy(w13)+energy(w14)+energy(w15)) / (energy(w10) +energy(w11) +energy(w12) ), 0.5)/(2*M_PI);
		    reduce(&output_var,1);
        
    	}

		else if(!strcmp(param.timevar_vars.name[i],"bxby")) {
			// incompressible maxwell stress
			output_var=compute_2correlation(wr5,wr6);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[i],"hc")) {
			// cross helicity (u.b)
			output_var=compute_2correlation(wr1,wr5)+compute_2correlation(wr2,wr6)+compute_2correlation(wr3,wr7);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[i],"hm")) {
			// magnetic helicity
			// Compute vector potential
			for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
				w10[j] = I * ik2t[j] * (ky[j] * fldi.bz[j] - kz[j] * fldi.by[j] );
				w11[j] = I * ik2t[j] * (kz[j] * fldi.bx[j] - kxt[j]* fldi.bz[j] );
				w12[j] = I * ik2t[j] * (kxt[j]* fldi.by[j] - ky[j] * fldi.bx[j] );
			}
	
			gfft_c2r(w10);
			gfft_c2r(w11);
			gfft_c2r(w12);
	
			for( j = 0 ; j < 2*NTOTAL_COMPLEX ; j++) {
				wr10[j] = wr10[j] / ((double) NTOTAL );
				wr11[j] = wr11[j] / ((double) NTOTAL );
				wr12[j] = wr12[j] / ((double) NTOTAL );
			}

			output_var=compute_2correlation(wr10,wr5)+compute_2correlation(wr11,wr6)+compute_2correlation(wr12,wr7);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[i],"j2")) {
			// total current rms
			// Compute current
			for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
				w10[j] = I * (ky[j]  * fldi.bz[j] - kz[j]  * fldi.by[j]);
				w11[j] = I * (kz[j]  * fldi.bx[j] - kxt[j] * fldi.bz[j]);
				w12[j] = I * (kxt[j] * fldi.by[j] - ky[j]  * fldi.bx[j]);
			}
			output_var=energy(w10)+energy(w11)+energy(w12);
			reduce(&output_var,1);
		}
		
		else if(!strcmp(param.timevar_vars.name[i],"az2")) {
			// Square of the vertical component of the vector potential
			for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
				w12[j] = I * ik2t[j] * (kxt[j]* fldi.by[j] - ky[j] * fldi.bx[j] );
			}
			
			output_var=energy(w12);
			
			reduce(&output_var,1);
		}
		
		else if(!strcmp(param.timevar_vars.name[i],"vxaz")) {
			// Source term in the streamfunction equation
			for( j = 0 ; j < NTOTAL_COMPLEX ; j++) {
				w12[j] = I * ik2t[j] * (kxt[j]* fldi.by[j] - ky[j] * fldi.bx[j] );
			}
			
			gfft_c2r(w12);
			
			for( j = 0 ; j < 2*NTOTAL_COMPLEX ; j++) {
				wr12[j] = wr12[j] / ((double) NTOTAL );
			}
			
			output_var=compute_2correlation(wr12,wr1);
			reduce(&output_var,1);
			
		}
#endif
#ifdef BOUSSINESQ
		else if(!strcmp(param.timevar_vars.name[i],"et")) {
			// thermal energy
			output_var = energy(fldi.th);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[i],"thmax")) {
			// maximum of th
			output_var=find_max(wr4);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[i],"thmin")) {
			// minimum of th
			output_var=find_min(wr4);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[i],"thvx")) {
			// turbulent heat flux in x
			output_var=compute_2correlation(wr4,wr1);
			reduce(&output_var,1);
		}
		else if(!strcmp(param.timevar_vars.name[i],"thvz")) {
			// turbulent heat flux in z
			output_var=compute_2correlation(wr4,wr3);
			reduce(&output_var,1);
		}
#endif
#ifdef COMPRESSIBLE
		else if(!strcmp(param.timevar_vars.name[i],"dmax")) {
			// maximum of density
			output_var=find_max(wr4);
			reduce(&output_var,2);
		}
		else if(!strcmp(param.timevar_vars.name[i],"dmin")) {
			// minimum of density
			output_var=find_min(wr4);
			reduce(&output_var,3);
		}
		else if(!strcmp(param.timevar_vars.name[i],"dvxvy")) {
			// compressible reynolds stress
			output_var=compute_3correlation(wr1,wr2,wr4);
			reduce(&output_var,1);
		}
#endif

		else {
			if(!warning_flag) {
				ERROR_HANDLER(ERROR_WARNING,"Unable to produce the requested data\n");
				MPI_Printf("Timevar output string ''%s'' is unknown\n",param.timevar_vars.name[i]);
				warning_flag=1;
			}
			output_var=0.0;
		}
		
		// Write the variable in the timevar file
		if(rank==0) {
			fprintf(ht,"%08e\t",output_var);
			if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing timevar file");
		}
		
#ifdef MPI_SUPPORT
		MPI_Barrier(MPI_COMM_WORLD);
#endif
	}
		
	if(rank==0) {
		fprintf(ht,"\n");
		fclose(ht);
	}

#ifdef OSCILLATORY_SHEAR
#ifdef SHELL_R_XY
	FILE *RxyFile;

	// Open the timevar file
	if(rank==0) {
		RxyFile=fopen("R_xy_Fourier.dat", "a");
		fprintf(RxyFile, "%08e\t", t);
		for(i=0; i < (int) ceil(kmax / OUTPUT_SPECTRUM_K_BIN); i++){
			// N.B.: No reduce() operation is needed here, since the make_spectrum()
			//       function that generates the integrated_..._shell[i] arrays
			//       contains a reduce() step at the end.

			// Write the cumulative fourier integral for each shell
			fprintf(RxyFile, "%08e\t", integrated_sin_Rxy_shell[i]);
			fprintf(RxyFile, "%08e\t", integrated_cos_Rxy_shell[i]);
		}
		fprintf(RxyFile, "\n");
		fclose(RxyFile);
	}
#endif
#endif

	DEBUG_END_FUNC;
	
	return;
}

/**************************************************************************************/
/** 
	Remove the timevar file (if exists) to start from a fresh one.
*/
/**************************************************************************************/
void init_timevar() {
	FILE* ht;
	int i;
	DEBUG_START_FUNC;
	
	if(rank==0) {
		ht=fopen("timevar","w");
		// print a line with the fields we are going to write
		
		for( i = 0 ; i < param.timevar_vars.length ; i++) {
			fprintf(ht,"%s\t\t",param.timevar_vars.name[i]);
		}
		fprintf(ht,"\n");
		
		fclose(ht);
#ifdef WITH_PARTICLES
		ht=fopen("partvar","w");
		fclose(ht);
#endif
	}

#ifdef OSCILLATORY_SHEAR
#ifdef SHELL_R_XY
	FILE *RxyFile;
	if(rank==0){
		RxyFile = fopen("R_xy_Fourier.dat", "w");
		fprintf(RxyFile, "#t\t");
		for(i=0; i < (int) ceil(kmax / OUTPUT_SPECTRUM_K_BIN); i++){
			// Write the wavenumber of each bin
			fprintf(RxyFile, "%08e\t", i * OUTPUT_SPECTRUM_K_BIN);
			fprintf(RxyFile, "%08e\t", i * OUTPUT_SPECTRUM_K_BIN);
		}
		fprintf(RxyFile, "\n");
		fclose(RxyFile);
	}
#endif
#endif

	DEBUG_END_FUNC;
	return;
}
/*********************************************************************************/
/** 
	Kinematic dynamo growth rates & energy at each mode  
	@param fldi Field structure from which the statistical quantities are derived.
	@param t Current time of the simulation
*/
/*********************************************************************************/

void output_rate(const struct Field fldi,
				 const double t) {
					
	FILE *ht;
    int i,j,k, rank_i;
    char  filename[50];
       
	DEBUG_START_FUNC;

    //Print each rank to a different file 
    for(rank_i=0; rank_i < NPROC; rank_i++){
        if(rank==rank_i) {
            sprintf(filename,"rates_r%02i.txt",rank_i);
            ht=fopen(filename,"a");
            for( i = 0; i < NX_COMPLEX/NPROC; i++) {
                for( j = 0; j < NY_COMPLEX; j++) {
                    for( k = 0; k < NZ_COMPLEX; k++) {
                        w16[i] = pow( pow(fabs(fldi.bx[IDX3D]),2) + pow(fabs(fldi.by[IDX3D]),2) + pow(fabs(fldi.bz[IDX3D]),2)  , 0.5);
                        if(isfinite(creal(w17[IDX3D]))==1 && creal(w17[IDX3D])>0.0){
                            fprintf(ht, "%f\t%f\t%f\t%f\t%e\t%e\n",kx[IDX3D],ky[IDX3D],kz[IDX3D],t,creal(w16[IDX3D]), creal(w17[IDX3D]));
                            if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing rate file");
                        }        
                    } 
                } 
            } 
            fclose(ht); 
        }
    }
} 
/**************************************************************************************/
/** 
	Remove the rates file (if exists) to start from a fresh one.
*/
/**************************************************************************************/
void init_rate(){

    FILE* ht;
    int i;
    char  filename[50];
	DEBUG_START_FUNC;

    for (i=0; i<NPROC; i++){
	    if(rank==i) {
            sprintf(filename,"rates_r%02i.txt",i);
		    ht=fopen(filename,"w");
		    // print a line with the fields we are going to write
            fprintf(ht,"KX\tKY\tKZ\tTIME\tEM\tGR\n");
		    fclose(ht);
        }
    }
   
	DEBUG_END_FUNC;
	return;
}

/**************************************************************************************/
/** 
	                            Isotropy T code
*/
/**************************************************************************************/

/* Start of isotropy T code */
void init_isotropyT(){
	FILE *ht;
	if (rank==0){
		ht = fopen("isotropyT.dat", "w");
		fprintf(ht, "#t\tT1111\tT1112\tT1123\tTijkl*Tijkl\tu4\n");
		fclose(ht);
	}
}

void output_isotropyT(const struct Field fld, const double t){
	double *returnT = malloc(sizeof(double)*5);


	isotropyT(fld, returnT);

	if (rank == 0){
		FILE *ht = fopen("isotropyT.dat", "a");

		fprintf(ht, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t, returnT[0], returnT[1], returnT[2], returnT[3], returnT[4]);

		fclose(ht);
		free(returnT);
	}

}
/* End of isotropy T code */

/**************************************************************************************/
/** 
	                            Timeseries code t
    The purpose of this function is to output timeseries 
    of the velocity field at up to two points in the file
*/
/**************************************************************************************/

double last_timeseries_time;

void output_timeseries(double *w1, double *w2, double *w3, const double t){
	
	FILE *ht;
	double ux, uy, uz;
	int idx3d;

	if(((t-last_timeseries_time) >= param.timeseries_time) || (last_timeseries_time == -1.0)) {
		if(param.output_timeseries_1){
#ifdef MPI_SUPPORT
			if(rank == (param.timeseries_1ny*NPROC) / NY){
#endif
				ht = fopen("timeseries1.dat", "a");
#ifdef MPI_SUPPORT
				// If MPI is enabled, and FFTW has MPI support, the X and Y directions
				// get transposed during the call to gfft_c2r_t due to the definition
				// of the c2rfft_mpi_t plan
				idx3d = param.timeseries_1nz + param.timeseries_1nx*(NZ+2)
				       +(param.timeseries_1ny - rank*NY/NPROC)*NX*(NZ+2);
#else
				// The X and Y directions don't get transposed in this case
				idx3d = param.timeseries_1nz + param.timeseries_1ny*(NZ+2)
				       +(param.timeseries_1nx)*NY*(NZ+2);
#endif

				ux = w1[idx3d]/((double)NTOTAL);
				uy = w2[idx3d]/((double)NTOTAL);
				uz = w3[idx3d]/((double)NTOTAL);

				fprintf(ht, "%lf\t%lf\t%lf\t%lf\n", t, ux, uy, uz);
		
				fclose(ht);
#ifdef MPI_SUPPORT
			}
#endif
		}

		if(param.output_timeseries_2){
#ifdef MPI_SUPPORT
			if(rank == (param.timeseries_2ny*NPROC) / NY){
#endif
				ht = fopen("timeseries2.dat", "a");
#ifdef MPI_SUPPORT
				// If MPI is enabled, and FFTW has MPI support, the X and Y directions
				// get transposed during the call to gfft_c2r_t due to the definition
				// of the c2rfft_mpi_t plan
				idx3d = param.timeseries_2nz + param.timeseries_2nx*(NZ+2)
				       +(param.timeseries_2ny - rank*NY/NPROC)*NX*(NZ+2);
#else
				// The X and Y directions don't get transposed in this case
				idx3d = param.timeseries_2nz + param.timeseries_2ny*(NZ+2)
				       +(param.timeseries_2nx)*NY*(NZ+2);
#endif

				ux = w1[idx3d]/((double)NTOTAL);
				uy = w2[idx3d]/((double)NTOTAL);
				uz = w3[idx3d]/((double)NTOTAL);

				fprintf(ht, "%lf\t%lf\t%lf\t%lf\n", t, ux, uy, uz);

				fclose(ht);
#ifdef MPI_SUPPORT
			}
#endif
		}
		last_timeseries_time = t;
	}

}

void init_timeseries(){
	FILE *ht;
	double x,y,z;

	last_timeseries_time= -1.0;

	if(param.output_timeseries_1){
		// Map the x cordinate to the interval [0, Lx)
		while(param.timeseries_1x < 0){
			param.timeseries_1x += param.lx;
		}
		while (param.timeseries_1x >= param.lx){
			param.timeseries_1x -= param.lx;
		}
		// Map the y cordinate to the interval [0, Ly)
		while(param.timeseries_1y < 0){
			param.timeseries_1y += param.ly;
		}
		while (param.timeseries_1y >= param.ly){
			param.timeseries_1y -= param.ly;
		}
		// Map the z cordinate to the interval [0, Lz)
		while(param.timeseries_1z < 0){
			param.timeseries_1z += param.lz;
		}
		while (param.timeseries_1z >= param.lz){
			param.timeseries_1z -= param.lz;
		}

		param.timeseries_1nx = (int)(NX*param.timeseries_1x/param.lx);
		param.timeseries_1ny = (int)(NY*param.timeseries_1y/param.ly);
		param.timeseries_1nz = (int)(NZ*param.timeseries_1z/param.lz);

		x = param.lx*param.timeseries_1nx/NX;
		y = param.ly*param.timeseries_1ny/NY;
		z = param.lz*param.timeseries_1nz/NZ;

		if(rank==0){
			ht = fopen("timeseries1.dat", "w");
			fprintf(ht, "#timeseries from (%lf, %lf, %lf)\n", x, y, z);
			fprintf(ht, "#t\tux\tuy\tuz\n");
			fclose(ht);
		}
	}
	if(param.output_timeseries_2){
		// Map the x cordinate to the interval [0, Lx)
		while(param.timeseries_2x < 0){
			param.timeseries_2x += param.lx;
		}
		while (param.timeseries_2x >= param.lx){
			param.timeseries_2x -= param.lx;
		}
		// Map the y cordinate to the interval [0, Ly)
		while(param.timeseries_2y < 0){
			param.timeseries_2y += param.ly;
		}
		while (param.timeseries_2y >= param.ly){
			param.timeseries_2y -= param.ly;
		}
		// Map the z cordinate to the interval [0, Lz)
		while(param.timeseries_2z < 0){
			param.timeseries_2z += param.lz;
		}
		while (param.timeseries_2z >= param.lz){
			param.timeseries_2z -= param.lz;
		}

		param.timeseries_2nx = (int)(NX*param.timeseries_2x/param.lx);
		param.timeseries_2ny = (int)(NY*param.timeseries_2y/param.ly);
		param.timeseries_2nz = (int)(NZ*param.timeseries_2z/param.lz);

		x = param.lx*param.timeseries_2nx/NX;
		y = param.ly*param.timeseries_2ny/NY;
		z = param.lz*param.timeseries_2nz/NZ;

		if(rank==0){
			ht = fopen("timeseries2.dat", "w");
			fprintf(ht, "#timeseries from (%lf, %lf, %lf)\n", x, y, z);
			fprintf(ht, "#t\tux\tuy\tuz\n");
			fclose(ht);
		}
	}
}
/* End of timeseries code */
