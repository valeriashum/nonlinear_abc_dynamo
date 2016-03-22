/*	
	This file is part of the Snoopy code.

    Snoopy code is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Snoopy code is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Snoopy code.  If not, see <http://www.gnu.org/licenses/>.
*/



#include "common.h"
#include "gfft.h"
#include "output/output_dump.h"
#include "symmetries.h"

#include "debug.h"



/** Allow one to init a structure in real space using ordinary defined x,y,z coordinates */

void init_SpatialStructure(struct Field fldi) {
	double *x,*y,*z;
	int i,j,k;
	
	/*******************************************************************
	** This part does not need to be modified **************************
	********************************************************************/
	// Allocate coordinate arrays
	x = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (x == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for x allocation");
	
	y = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (y == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for y allocation");
	
	z = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (z == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for z allocation");

	// Initialize the arrays
	
	for(i = 0 ; i < NX/NPROC ; i++) {
		for(j = 0 ; j < NY ; j++) {
			for(k = 0 ; k < NZ ; k++) {
				x[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - param.lx / 2 + (param.lx * (i + rank * NX / NPROC)) / NX;
				y[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - param.ly / 2 + (param.ly * j ) / NY;
				z[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - param.lz / 2 + (param.lz * k ) / NZ;
			}
		}
	}
	
	// Initialize the extra points (k=NZ and k=NZ+1) to zero to prevent stupid things from happening...
	for(i = 0 ; i < NX/NPROC ; i++) {
		for(j = 0 ; j < NY ; j++) {
			for(k = NZ ; k < NZ + 2 ; k++) {
				x[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				y[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				z[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
			}
		}
	}
	
	// Init work array to zero
	for(i = 0 ; i < NX/NPROC ; i++) {
		for(j = 0 ; j < NY ; j++) {
			for(k = 0 ; k < NZ + 2 ; k++) {
				wr1[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				wr2[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				wr3[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				wr4[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				wr5[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				wr6[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
			}
		}
	}
	
	/*******************************************************************
	** This part can be modified              **************************
	********************************************************************/
	
	// The velocity field vx,vy,vz is stored in wr1,wr2,wr3
	// The magnetic field bx,by,bz is stored in wr4,wr5,wr6 (ignored if MHD is not set)
	
	for(i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		// Example: init a flux tube in the x direction+a vertical displacement
	  //	wr4[i] = exp(-(y[i]*y[i]+z[i]*z[i])*20.0);
	  //	wr3[i] = 0.5*cos(x[i] * 2.0 * M_PI);

		// Example: twisted flux tube + vertical displacement
	//	wr4[i] = exp(-(y[i]*y[i]+z[i]*z[i])/(0.2*0.2));
	//	wr5[i] = fabs(z[i])*1.0*wr4[i];
	//	wr6[i] = -fabs(y[i])*1.0*wr4[i];
	//	wr3[i] = 0.5*cos(x[i] * 2.0 * M_PI);
	//	if (i==3*NY*(NZ+2)) fprintf(stderr," %d %e",i,x[i]);

	}
	
	/*******************************************************************
	** This part does not need to be modified **************************
	********************************************************************/
	// Fourier transform everything
	gfft_r2c(wr1);
	gfft_r2c(wr2);
	gfft_r2c(wr3);
	gfft_r2c(wr4);
	gfft_r2c(wr5);
	gfft_r2c(wr6);
	
	// Transfer data in the relevant array (including dealiasing mask)
	for(i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		fldi.vx[i] += w1[i] * mask[i];
		fldi.vy[i] += w2[i] * mask[i];
		fldi.vz[i] += w3[i] * mask[i];
#ifdef MHD
		fldi.bx[i] += w4[i] * mask[i];
		fldi.by[i] += w5[i] * mask[i];
		fldi.bz[i] += w6[i] * mask[i];
#endif
	}
	
	// free memory
	fftw_free(x);
	fftw_free(y);
	fftw_free(z);
	
	//done
	return;
}


void init_KidaVortex(struct Field fldi) {
	double a = param.vortex_a;
	double b = param.vortex_b;
	
	int i,j,k;
	
	double w0, x, y;
	double chi;
	
	chi = b / a;
	w0 = 1.0/chi*(chi + 1.0)/(chi-1.0);			// According to Kida!
	
	for(i = 0 ; i < NX/NPROC ; i++) {
		x = - param.lx / 2 + (param.lx * (i + rank * NX / NPROC)) / NX;
		for(j = 0 ; j < NY ; j++) {
			y = - param.ly / 2 + (param.ly * j) / NY;
#ifdef WITH_2D
			if(x * x / (a * a) + y * y / (b * b) < 1) {
					// we are in the vortex
					wr1[j + (NY+2) * i] = -w0;
			}
			else {
				wr1[j + (NY+2) * i] = 0.0;
			}
#else
			for(k = 0 ; k < NZ ; k++) {
				if(x * x / (a * a) + y * y / (b * b) < 1) {
					// we are in the vortex
					wr1[k + j*(NZ+2) + (NZ+2) * NY * i] = -w0;
				}
				else {
					wr1[k + j*(NZ+2) + (NZ+2) * NY * i] = 0.0;
				}
			}
#endif
		MPI_Printf("From initflow.c: check x[i]=%f\n", x);

		}

	}
	
	// transform
	gfft_r2c(wr1);
	
	for(i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		fldi.vx[ i ] +=  I * ky[i] * w1[i] * ik2t[i];
		fldi.vy[ i ] += -I * kxt[i] * w1[i] * ik2t[i];
	}
	
	// done
	return;
}

/************************************/
/** Init some crazy structure involving
/** A kida vortex and a vertical structure
/** for the field */
/***********************************/
void init_Bench(struct Field fldi) {
	const double a = 0.3;
	const double b = 0.4;
	
	int i,j,k;
	
	double w0, x, y;
	double chi;
	
	chi = b / a;
	w0 = 1.0/chi*(chi + 1)/(chi-1.0);			// According to Kida!
	
	for(i = 0 ; i < NX/NPROC ; i++) {
		x = - param.lx / 2. + (param.lx * (i + rank * NX / NPROC)) / NX;
		for(j = 0 ; j < NY ; j++) {
			y = - param.ly / 2. + (param.ly * j) / NY;
			for(k = 0 ; k < NZ ; k++) {
				if(x * x / (a * a) + y * y / (b * b) < 1) {
					// we are in the vortex
					wr1[k + j*(NZ+2) + (NZ+2) * NY * i] = -w0;
				}
				else {
					wr1[k + j*(NZ+2) + (NZ+2) * NY * i] = 0.0;
				}
			}
		}
	}
	
	// transform
	gfft_r2c(wr1);
	
	for(i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		fldi.vx[ i ] +=  I * ky[i] * w1[i] * ik2t[i];
		fldi.vy[ i ] += -I * kxt[i] * w1[i] * ik2t[i];
	}
	
	// Brake vertical symmetry
	if(rank==0) {
		fldi.vx[1] = 1000.0 / NTOTAL;
		fldi.vy[1] = 1000.0 / NTOTAL;
#ifdef MHD
		fldi.bx[1] = 1000.0 / NTOTAL;
		fldi.by[1] = 1000.0 / NTOTAL;
#endif
	}
	// done
	return;
}


void init_LargeScaleNoise(struct Field fldi) {
	int i,j,k;
	int num_force=0;
	int total_num_force;
	double fact;
	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				if(pow(k2t[ IDX3D ], 0.5) / ( 2.0*M_PI ) < 1.0 / param.noise_cut_length) {
					fldi.vx[ IDX3D ] += param.per_amplitude_large * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
					fldi.vy[ IDX3D ] += param.per_amplitude_large * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
					fldi.vz[ IDX3D ] += param.per_amplitude_large * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
#ifdef MHD
					fldi.bx[ IDX3D ] += param.per_amplitude_large * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
					fldi.by[ IDX3D ] += param.per_amplitude_large * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
					fldi.bz[ IDX3D ] += param.per_amplitude_large * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
#endif
					if(mask[IDX3D] > 0) num_force++;
				}
			}
		}
	}
	
	// Get the total number of forced scales.
#ifdef MPI_SUPPORT
	MPI_Allreduce( &num_force, &total_num_force, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
	total_num_force=num_force;
#endif
	
	fact=pow(total_num_force,0.5);
	
	// Divide by the total number of modes
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				fldi.vx[ IDX3D ] = fldi.vx[ IDX3D ] / fact;
				fldi.vy[ IDX3D ] = fldi.vy[ IDX3D ] / fact;
				fldi.vz[ IDX3D ] = fldi.vz[ IDX3D ] / fact;
#ifdef MHD
				fldi.bx[ IDX3D ] = fldi.bx[ IDX3D ] / fact;
				fldi.by[ IDX3D ] = fldi.by[ IDX3D ] / fact;
				fldi.bz[ IDX3D ] = fldi.bz[ IDX3D ] / fact;
#endif
			}
		}
	}
	
  enforce_complex_symm(fldi);  
}

/******************************************
** Large scale 2D (x,y) noise *************
*******************************************/

void init_LargeScale2DNoise(struct Field fldi) {
	int i,j,k;
	int num_force=0;
	int total_num_force;
	double fact;
	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			k=0;
			if(kz[ IDX3D ] == 0.0) {
				if(pow(k2t[ IDX3D ], 0.5) / ( 2.0*M_PI ) < 1.0 / param.noise_cut_length_2D) {
					fldi.vx[ IDX3D ] += param.per_amplitude_large_2D * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
					fldi.vy[ IDX3D ] += param.per_amplitude_large_2D * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
#ifdef MHD
					fldi.bx[ IDX3D ] += param.per_amplitude_large_2D * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
					fldi.by[ IDX3D ] += param.per_amplitude_large_2D * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL;
#endif
					if(mask[IDX3D] > 0) num_force++;
				}
			}
		}
	}
	
	// Get the total number of forced scales.
#ifdef MPI_SUPPORT
	MPI_Allreduce( &num_force, &total_num_force, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
	total_num_force=num_force;
#endif
	
	fact=pow(total_num_force,0.5);
	
	// Divide by the total number of modes
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			k=0;
			if(kz[ IDX3D ] == 0.0) {
				fldi.vx[ IDX3D ] = fldi.vx[ IDX3D ] / fact;
				fldi.vy[ IDX3D ] = fldi.vy[ IDX3D ] / fact;
#ifdef MHD
				fldi.bx[ IDX3D ] = fldi.bx[ IDX3D ] / fact;
				fldi.by[ IDX3D ] = fldi.by[ IDX3D ] / fact;
#endif
			}
		}
	}
	
  enforce_complex_symm(fldi);  
}

//***********************************************************************************************************************************************

void init_White_Noise(struct Field fldi) {
	int i,j,k;
	double fact;
	
	// Excite (2/3)^3*NTOTAL modes
	fact = pow(27.0/8.0*NTOTAL, 0.5);

   	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {                
//				fldi.vx[ IDX3D ] += param.per_amplitude_noise * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * fact;
//				fldi.vy[ IDX3D ] += param.per_amplitude_noise * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * fact;
//				fldi.vz[ IDX3D ] += param.per_amplitude_noise * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * fact;
#ifdef MHD
                fldi.bx[ IDX3D ] += param.per_amplitude_noise * mask[IDX3D] * randm() * cexp( I * 0.25*M_PI ) * fact;
				fldi.by[ IDX3D ] += param.per_amplitude_noise * mask[IDX3D] * randm() * cexp( I * 0.25*M_PI ) * fact;
				fldi.bz[ IDX3D ] += param.per_amplitude_noise * mask[IDX3D] * randm() * cexp( I * 0.25*M_PI ) * fact;
                
               /* if(creal(fldi.bx[IDX3D])!= 0.0 || cimag(fldi.bx[IDX3D])!= 0.0 || creal(fldi.by[IDX3D])!= 0.0 || cimag(fldi.by[IDX3D])!= 0.0 || creal(fldi.bz[IDX3D])!= 0.0 || cimag(fldi.bz[IDX3D])!= 0.0){
                    printf("bx=%e + i%e, i=%d,j=%d,k=%d, rank=%d\n", creal(fldi.bx[IDX3D]), cimag(fldi.bx[IDX3D]),i,j,k, rank);
                    printf("by=%e + i%e, i=%d,j=%d,k=%d, rank=%d\n", creal(fldi.by[IDX3D]), cimag(fldi.by[IDX3D]),i,j,k, rank);
                    printf("bz=%e + i%e, i=%d,j=%d,k=%d, rank=%d\n", creal(fldi.bz[IDX3D]), cimag(fldi.bz[IDX3D]),i,j,k, rank);
                }*/

                //if (i==0 && j==1 && k==1) printf("rank=%d, idx3d=%d, amp=%f, mask[IDX3D]=%d,randm()=%e,cexp( I * 2.0*M_PI*randm() )=%f\n",rank, IDX3D, param.per_amplitude_noise, mask[IDX3D], randm(), cexp( I * 2.0*M_PI*randm() ) );
#endif
			}
		}
	}

	enforce_complex_symm(fldi);
/*    // check non-zero fields
for(i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		if (creal(fldi.bx[i])!= 0.0 || cimag(fldi.bx[i])!= 0.0) {
	        printf("in rank=%d,for kx =%f, ky=%f, kz=%f, bx= %e+ i%e\n",rank, kx[i],ky[i],kz[i], creal(fldi.bx[i]), cimag(fldi.bx[i]));
        }
        if (creal(fldi.by[i])!= 0.0 || cimag(fldi.by[i])!= 0.0) {
	        printf("in rank=%d,for kx =%f, ky=%f, kz=%f, by= %e+ i%e\n",rank, kx[i],ky[i],kz[i], creal(fldi.by[i]), cimag(fldi.by[i]));
        }
        if (creal(fldi.bz[i])!= 0.0 || cimag(fldi.bz[i])!= 0.0) {
	        printf("in rank=%d,for kx =%f, ky=%f, kz=%f, bz= %e+ i%e\n",rank, kx[i],ky[i],kz[i], creal(fldi.bz[i]), cimag(fldi.bz[i]));
        }       
    }*/
  
}//***********************************************************************************************************************************************

void init_No_Noise(struct Field fldi) {
	int i,j,k;
	
	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
                fldi.vx[ IDX3D ] = 0.0;
                fldi.vy[ IDX3D ] = 0.0;
                fldi.vz[ IDX3D ] = 0.0;

#ifdef MHD
				fldi.bx[ IDX3D ] = 0.0;
				fldi.by[ IDX3D ] = 0.0;
				fldi.bz[ IDX3D ] = 0.0;
#endif
			}
		}
	}

    enforce_complex_symm(fldi);
  
}

//***********************************************************************************************************************************************
void init_dominant_scale(struct Field fldi) {
	int i,j,k;
#ifdef MHD	


	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
                if (rank == 0){
				// Set higher amplitude to the modes of i,j,k
		    		//if ((i==1 && j==1 && k==1) || (i==1 && j==0 && k==1) || (i==1 && j==1 && k==0) || (i==1 && j==0 && k==0) || (i==2 && j==1 && k==1) || (i==0 && j==1 && k==1)) {
                    // If both j=k=0, it will set Bx=0 due to solenoidal condition
            	    if (i==1 && j==1 && k==1) {
				        fldi.bx[IDX3D] += param.per_amplitude_noise * cexp( I * 2.1*M_PI ); 
				    	fldi.by[IDX3D] += param.per_amplitude_noise * cexp( I * 1.2*M_PI ); 
					    fldi.bz[IDX3D] += param.per_amplitude_noise * cexp( I * 0.4*M_PI );
                    }
                }
	        }
		}
	}
#endif

enforce_complex_symm(fldi);
}


//***********************************************************************************************************************************************
void init_ABC_flow(struct Field fldi){
	// Initialise the flow to (A sinz + C cosy, B sinx + A cos z, C siny+ B cosx)

	int i,j,k;	// Used to address the fourier modes via the macro IDX3D
       
  
   // SET ABC wavenumber entries
    
	// vx += sin(z), vy += cos(z)
	if(rank == 0){
		i=j=0; k=param.ABC_flow_kz;
		fldi.vx[IDX3D] += ((double)NTOTAL)*0.5*param.ABC_flow_A*(-I); 
		fldi.vy[IDX3D] += ((double)NTOTAL)*0.5*param.ABC_flow_A*(+1);
        k = NZ_COMPLEX - param.ABC_flow_kz; 
	}

	// vx += cos(y), vz += sin(y)
	if(rank == 0){
		i=k=0; j=param.ABC_flow_ky;
		fldi.vx[IDX3D] += ((double)NTOTAL)*0.5*param.ABC_flow_C*(+1);
		fldi.vz[IDX3D] += ((double)NTOTAL)*0.5*param.ABC_flow_C*(-I);
		j = NY_COMPLEX - param.ABC_flow_ky;
		fldi.vx[IDX3D] += ((double)NTOTAL)*0.5*param.ABC_flow_C*(+1);
		fldi.vz[IDX3D] += ((double)NTOTAL)*0.5*param.ABC_flow_C*(+I);
	}

	// vy += sin(x), vz += cos(x)
	if(rank == (param.ABC_flow_kx/(NX_COMPLEX/NPROC))){
		i= (param.ABC_flow_kx%(NX_COMPLEX/NPROC)); j=k=0;
		fldi.vy[IDX3D] += ((double)NTOTAL)*0.5*param.ABC_flow_B*(-I);
		fldi.vz[IDX3D] += ((double)NTOTAL)*0.5*param.ABC_flow_B*(+1);
	}
    
	if(rank == ((NX_COMPLEX - param.ABC_flow_kx)/(NX_COMPLEX/NPROC))){
		i = ((NX_COMPLEX - param.ABC_flow_kx)%(NX_COMPLEX/NPROC)); j=k=0;
		fldi.vy[IDX3D] += ((double)NTOTAL)*0.5*param.ABC_flow_B*(+I);
		fldi.vz[IDX3D] += ((double)NTOTAL)*0.5*param.ABC_flow_B*(+1);
	}
        
	enforce_complex_symm(fldi);
    // I don't really know what this does, but all the other init functions were doing it
}
/***************************************************************************************************************/

void init_MeanField(struct Field fldi) {
    int i,j,k;
#ifdef MHD
	if(rank==0) {
		fldi.bx[0] = param.bx0 * ((double) NTOTAL);
		fldi.by[0] = param.by0 * ((double) NTOTAL);
		fldi.bz[0] = param.bz0 * ((double) NTOTAL);
	}  
#endif
}
/***************************************************************************************************************/

void init_Magnetic_field(struct Field fldi) {
	double *x,*y,*z;
	int i,j,k;
#ifdef MHD

	/*******************************************************************
	** This part does not need to be modified **************************
	********************************************************************/
	// Allocate coordinate arrays
	x = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (x == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for x allocation");
	
	y = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (y == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for y allocation");
	
	z = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (z == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for z allocation");

	// Initialize the arrays
	
	for(i = 0 ; i < NX/NPROC ; i++) {
		for(j = 0 ; j < NY ; j++) {
			for(k = 0 ; k < NZ ; k++) {
				x[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - param.lx / 2 + (param.lx * (i + rank * NX / NPROC)) / NX;
				y[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - param.ly / 2 + (param.ly * j ) / NY;
				z[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - param.lz / 2 + (param.lz * k ) / NZ;
			}
		}
	}
	
	// Initialize the extra points (k=NZ and k=NZ+1) to zero to prevent stupid things from happening...
	for(i = 0 ; i < NX/NPROC ; i++) {
		for(j = 0 ; j < NY ; j++) {
			for(k = NZ ; k < NZ + 2 ; k++) {
				x[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				y[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				z[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
			}
		}
	}
	
	// Init work array to zero
	for(i = 0 ; i < NX/NPROC ; i++) {
		for(j = 0 ; j < NY ; j++) {
			for(k = 0 ; k < NZ + 2 ; k++) {
				wr4[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				wr5[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				wr6[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
			}
		}
	}
	
	/*******************************************************************
	** This part can be modified              **************************
	********************************************************************/
	
	// The magnetic field bx,by,bz is stored in wr4,wr5,wr6
	
	for(i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		// Example: From Bouya (2013) – same periodicity as the flow
		wr4[i] = sin(param.ABC_flow_kz*z[i]) - cos(param.ABC_flow_ky*y[i]);
		wr5[i] = sin(param.ABC_flow_kx*x[i]) - cos(param.ABC_flow_kz*z[i]);
		wr6[i] = sin(param.ABC_flow_ky*y[i]) - cos(param.ABC_flow_kx*x[i]);
	}
	
	/*******************************************************************
	** This part does not need to be modified **************************
	********************************************************************/
	// Fourier transform everything
	gfft_r2c(wr4);
	gfft_r2c(wr5);
	gfft_r2c(wr6);
	// Transfer data in the relevant array (including dealiasing mask)
	for(i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		fldi.bx[i] += w4[i] * mask[i];
		fldi.by[i] += w5[i] * mask[i];
		fldi.bz[i] += w6[i] * mask[i];
	}
	
	// free memory
	fftw_free(x);
	fftw_free(y);
	fftw_free(z);
	
	//done
	return;
#endif
}
void init_modified_ABC_flow(struct Field fldi){
	int ii,i,j,k;	    // Used to address the fourier modes via the macro IDX3D
    double *x,*y,*z;
    
    /*******************************************************************
	** This part does not need to be modified **************************
	********************************************************************/
	// Allocate coordinate arrays
	x = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (x == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for x allocation");
	
	y = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (y == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for y allocation");
	
	z = (double *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (z == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for z allocation");

	// Initialize the arrays
	
	for(i = 0 ; i < NX/NPROC ; i++) {
		for(j = 0 ; j < NY ; j++) {
			for(k = 0 ; k < NZ ; k++) {
				x[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - param.lx / 2 + (param.lx * (i + rank * NX / NPROC)) / NX;
				y[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - param.ly / 2 + (param.ly * j ) / NY;
				z[k + (NZ + 2) * j + (NZ + 2) * NY * i] = - param.lz / 2 + (param.lz * k ) / NZ;
			}
		}
	}
	
	// Initialize the extra points (k=NZ and k=NZ+1) to zero to prevent stupid things from happening...
	for(i = 0 ; i < NX/NPROC ; i++) {
		for(j = 0 ; j < NY ; j++) {
			for(k = NZ ; k < NZ + 2 ; k++) {
				x[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				y[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				z[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;

			}
		}
	}
	
	// Init work array to zero
	for(i = 0 ; i < NX/NPROC ; i++) {
		for(j = 0 ; j < NY ; j++) {
			for(k = 0 ; k < NZ + 2 ; k++) {
				wr1[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				wr2[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
				wr3[k + (NZ + 2) * j + (NZ + 2) * NY * i] = 0.0;
   			}
		}
   }
    /*******************************************************************
	** This part can be modified              **************************
	********************************************************************/
	
	// The velocity field vx,vy,vz is stored in wr1,wr2,wr3
   	for(i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
     
        wr1[i] = param.modified_ABC_flow_D*(
                    cos(x[i]/param.modified_ABC_flow_m)*(
                        param.modified_ABC_flow_A*param.modified_ABC_flow_kz*sin(z[i] /(double) param.modified_ABC_flow_kz ) +
                        param.modified_ABC_flow_C*param.modified_ABC_flow_ky*cos(y[i] /(double) param.modified_ABC_flow_ky)
                    ) + 
                    sin(x[i]/param.modified_ABC_flow_m)/param.modified_ABC_flow_m*(
                        0
                    )                                         
                ); 
        wr2[i] = param.modified_ABC_flow_D*(
                    cos(x[i]/param.modified_ABC_flow_m)*(
                        param.modified_ABC_flow_B*param.modified_ABC_flow_kx*sin(x[i] /(double) param.modified_ABC_flow_kx ) +
                        param.modified_ABC_flow_A*param.modified_ABC_flow_kz*cos(z[i] /(double) param.modified_ABC_flow_kz)
                    ) + 
                    sin(x[i]/param.modified_ABC_flow_m)/param.modified_ABC_flow_m*(
                        param.modified_ABC_flow_B*cos(x[i] /(double) param.modified_ABC_flow_kx) + 
                        param.modified_ABC_flow_C*sin(y[i] /(double) param.modified_ABC_flow_ky)
                    )                                         
                );
        wr3[i] = param.modified_ABC_flow_D*(
                    cos(x[i]/param.modified_ABC_flow_m)*(
                        param.modified_ABC_flow_C*param.modified_ABC_flow_ky*sin(y[i] /(double) param.modified_ABC_flow_ky ) +
                        param.modified_ABC_flow_B*param.modified_ABC_flow_kx*cos(x[i] /(double) param.modified_ABC_flow_kx)
                    ) + 
                    sin(x[i]/param.modified_ABC_flow_m)/param.modified_ABC_flow_m*(
                        -param.modified_ABC_flow_A*cos(z[i] /(double) param.modified_ABC_flow_kz) + 
                        -param.modified_ABC_flow_B*sin(x[i] /(double) param.modified_ABC_flow_kx)
                    )                                         
                );
    }
 


    gfft_r2c(wr1);
	gfft_r2c(wr2);
    gfft_r2c(wr3);
	
    for(i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		fldi.vx[i] = w1[i] * mask[i];
		fldi.vy[i] = w2[i] * mask[i];
		fldi.vz[i] = w3[i] * mask[i];
    }   

    for(i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		if (fabs(creal(fldi.vx[i])) > 1.0 || fabs(cimag(fldi.vx[i]))>1.0) {
	        printf("in rank=%d,for kx =%f, ky=%f, kz=%f, vx= %e+ i%e\n",rank, kx[i],ky[i],kz[i], creal(fldi.vx[i]), cimag(fldi.vx[i]));
        }
        if (fabs(creal(fldi.vy[i]))  > 1.0 || fabs(cimag(fldi.vy[i]))> 1.0) {
	        printf("in rank=%d,for kx =%f, ky=%f, kz=%f, vy= %e+ i%e\n",rank, kx[i],ky[i],kz[i], creal(fldi.vy[i]), cimag(fldi.vy[i]));
        }
        if (fabs(creal(fldi.vz[i]) ) > 1.0 || fabs(cimag(fldi.vz[i]))> 1.0) {
	        printf("in rank=%d,for kx =%f, ky=%f, kz=%f, vz= %e+ i%e\n",rank, kx[i],ky[i],kz[i], creal(fldi.vz[i]), cimag(fldi.vz[i]));
        }
    }
    return; 
}
/*****************************************************************************************************************/

/** Init the flow arrays... */	
void init_flow(struct Field fldi) {
	int i,n;
	int j,k;
	
	double dummy_var;
	
	DEBUG_START_FUNC;
	// Initialise vectors to 0
	
	for( n = 0 ; n < fldi.nfield ; n++) {
		for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
			fldi.farray[n][i] = 0.0;
		}
	}


	
#ifdef COMPRESSIBLE
	// Initialise the density to 1...
	if(rank==0) {
		fldi.d[0] = (double) NTOTAL;
	}
#endif
	if(param.init_large_scale_noise) init_LargeScaleNoise(fldi);
	
	if(param.init_large_scale_2D_noise) init_LargeScale2DNoise(fldi);

	if(param.init_vortex) init_KidaVortex(fldi);

	if(param.init_spatial_structure) init_SpatialStructure(fldi);

	if(param.init_white_noise) init_White_Noise(fldi);

	if(param.init_ABC_flow) init_ABC_flow(fldi);

	if(param.init_bench) init_Bench(fldi);

    if(param.init_modified_ABC_flow) init_modified_ABC_flow(fldi);
	
    if(param.init_mean_field) init_MeanField(fldi);

	if(param.init_dominant_scale) init_dominant_scale(fldi);

	if(param.init_Magnetic_field) init_Magnetic_field(fldi);

    if(param.init_No_Noise) init_No_Noise(fldi);

    if(param.init_dump) {
		read_dump(fldi, &dummy_var,"init.dmp");
		MPI_Printf("Initial conditions read successfully from the restart dump\n");
	}

#ifdef BOUNDARY_C
	boundary_c(fldi);
#endif

	projector(fldi.vx,fldi.vy,fldi.vz);

#ifdef MHD
	projector(fldi.bx,fldi.by,fldi.bz);
#endif

#ifdef WITH_PARTICLES
#ifdef WITH_ROTATION
		if(rank==0) {
			kappa_tau2 = 2.0*param.omega*(2.0*param.omega-param.shear) * param.particles_stime * param.particles_stime + (param.particles_dg_ratio + 1.0) * (param.particles_dg_ratio + 1.0);

	// This is a non trivial equilibrium for the particles+gas system
			fldi.vx[0] = param.particles_epsilon*param.particles_stime*param.particles_dg_ratio / kappa_tau2 * ( (double) NTOTAL);
			fldi.vy[0] = param.particles_epsilon*param.particles_dg_ratio*(1.0+param.particles_dg_ratio)/(2.0*param.omega*kappa_tau2) * ( (double) NTOTAL);
		}
#endif
#endif

// check non-zero fields
for(i = 0 ; i < NTOTAL_COMPLEX ; i++) {
	/*	if (creal(fldi.vx[i])!= 0.0 || cimag(fldi.vx[i])!= 0.0) {
	        printf("in rank=%d,for kx =%f, ky=%f, kz=%f, vx= %e+ i%e\n",rank, kx[i],ky[i],kz[i], creal(fldi.vx[i]), cimag(fldi.vx[i]));
        }
        if (creal(fldi.vy[i])!= 0.0 || cimag(fldi.vy[i])!= 0.0) {
	        printf("in rank=%d,for kx =%f, ky=%f, kz=%f, vy= %e+ i%e\n",rank, kx[i],ky[i],kz[i], creal(fldi.vy[i]), cimag(fldi.vy[i]));
        }
        if (creal(fldi.vz[i])!= 0.0 || cimag(fldi.vz[i])!= 0.0) {
	        printf("in rank=%d,for kx =%f, ky=%f, kz=%f, vz= %e+ i%e\n",rank, kx[i],ky[i],kz[i], creal(fldi.vz[i]), cimag(fldi.vz[i]));
        }

       /* if (creal(fldi.bx[i])!= 0.0 || cimag(fldi.bx[i])!= 0.0) {
	        printf("in rank=%d,for kx =%f, ky=%f, kz=%f, bx= %e+ i%e\n",rank, kx[i],ky[i],kz[i], creal(fldi.bx[i]), cimag(fldi.bx[i]));
        }
        if (creal(fldi.by[i])!= 0.0 || cimag(fldi.by[i])!= 0.0) {
	        printf("in rank=%d,for kx =%f, ky=%f, kz=%f, by= %e+ i%e\n",rank, kx[i],ky[i],kz[i], creal(fldi.by[i]), cimag(fldi.by[i]));
        }
        if (creal(fldi.bz[i])!= 0.0 || cimag(fldi.bz[i])!= 0.0) {
	        printf("in rank=%d,for kx =%f, ky=%f, kz=%f, bz= %e+ i%e\n",rank, kx[i],ky[i],kz[i], creal(fldi.bz[i]), cimag(fldi.bz[i]));
        }*/

        
    }

#ifdef DEBUG
	MPI_Printf("Initflow:\n");
	D_show_all(fldi);
	MPI_Printf("**************************************************************************************\n");
#endif	
	
	DEBUG_END_FUNC;
	
	return;
}
	
	
