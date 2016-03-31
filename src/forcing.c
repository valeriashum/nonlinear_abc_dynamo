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

#include <math.h>
#include <complex.h>

#include "common.h"

#ifdef FORCING

/*****************************************************
** Here are several forcing possibilities which can **
** be used. Please remove comments to use a specific**
** forcing                                          **
******************************************************/



/*****************************************************
** Random noise forcing ******************************
******************************************************/

void forcing(struct Field fldi,
			 double dt) {
			 
// Force random velocity field
	const double kf = 3.0 * M_PI * 2.0;
	const double deltakf = kf * 0.2;
	const double amplitude_forcing = 0.1;

	// Force all the vector
	
	
	int i,j,k;
	int num_force=0;
	int total_num_force;
	double fact;
	double q0;
	
	q0=pow(dt,0.5);
	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				if( (k2t[ IDX3D ]>(kf-deltakf)*(kf-deltakf)) && (k2t[ IDX3D ]<(kf+deltakf)*(kf+deltakf))) {
					w4[ IDX3D ] = amplitude_forcing * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL*q0;
					w5[ IDX3D ] = amplitude_forcing * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL*q0;
					w6[ IDX3D ] = amplitude_forcing * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL*q0;

					if(mask[IDX3D] > 0) num_force++;
				}
				else {
					w4[ IDX3D ] = 0.0;
					w5[ IDX3D ] = 0.0;
					w6[ IDX3D ] = 0.0;
				}
			}
		}
	}
	
  symmetrize_complex(w4);
  if(rank==0) w4[0]=0.0;
  symmetrize_complex(w5);
  if(rank==0) w5[0]=0.0;
  symmetrize_complex(w6);
  if(rank==0) w6[0]=0.0;
  
	// Get the total number of forced scales.
#ifdef MPI_SUPPORT
	MPI_Allreduce( &num_force, &total_num_force, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
	total_num_force=num_force;
#endif
	
	fact=pow(total_num_force,0.5);

	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				fldi.vx[ IDX3D ] += w4[ IDX3D ] / fact;
				fldi.vy[ IDX3D ] += w5[ IDX3D ] / fact;
				fldi.vz[ IDX3D ] += w6[ IDX3D ] / fact;
			}
		}
	}
	
	projector(fldi.vx,fldi.vy,fldi.vz);
	
	return;
}

#endif

/*****************************************************
********* U_III Forcing ******************************
******************************************************/
#ifdef U_III_FORCING

void u_iii_forcing(struct Field fldi, double dt) {
	double *x, *y, *z;		 	
	int i,j,k;

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
	
	// The velocity field vx,vy,vz is stored in wr1,wr2,wr3
   	for(i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
     
        wr4[i] = param.modified_ABC_flow_D*(
                    cos(x[i]/param.modified_ABC_flow_m)*(
                        param.modified_ABC_flow_A*param.modified_ABC_flow_kz*sin(z[i] /(double) param.modified_ABC_flow_kz ) +
                        param.modified_ABC_flow_C*param.modified_ABC_flow_ky*cos(y[i] /(double) param.modified_ABC_flow_ky)
                    ) + 
                    sin(x[i]/param.modified_ABC_flow_m)/param.modified_ABC_flow_m*(
                        0
                    )                                         
                ); 
        wr5[i] = param.modified_ABC_flow_D*(
                    cos(x[i]/param.modified_ABC_flow_m)*(
                        param.modified_ABC_flow_B*param.modified_ABC_flow_kx*sin(x[i] /(double) param.modified_ABC_flow_kx ) +
                        param.modified_ABC_flow_A*param.modified_ABC_flow_kz*cos(z[i] /(double) param.modified_ABC_flow_kz)
                    ) + 
                    sin(x[i]/param.modified_ABC_flow_m)/param.modified_ABC_flow_m*(
                        param.modified_ABC_flow_B*cos(x[i] /(double) param.modified_ABC_flow_kx) + 
                        param.modified_ABC_flow_C*sin(y[i] /(double) param.modified_ABC_flow_ky)
                    )                                         
                );
        wr6[i] = param.modified_ABC_flow_D*(
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
    gfft_r2c(wr4);
	gfft_r2c(wr5);
    gfft_r2c(wr6);
	
    for(i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w4[i] = w4[i] * mask[i];
		w5[i] = w5[i] * mask[i];
		w6[i] = w6[i] * mask[i];
    }   

	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				fldi.vx[ IDX3D ] +=  w4[ IDX3D ] * (1 - exp(-nu*k2t[IDX3D]*dt)) ;
				fldi.vy[ IDX3D ] +=  w5[ IDX3D ] * (1 - exp(-nu*k2t[IDX3D]*dt)) ;
				fldi.vz[ IDX3D ] +=  w6[ IDX3D ] * (1 - exp(-nu*k2t[IDX3D]*dt)) ;
			}
		}
	}

#ifdef U_III_FORCING_EXTRA
    gfft_c2r_t(w4);
	gfft_c2r_t(w5);
	gfft_c2r_t(w6);

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr10[i] = wr4[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr11[i] = wr5[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
#ifndef WITH_2D
		wr12[i] = wr6[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
#endif
		wr7[i] = wr4[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr8[i] = wr4[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr9[i] = wr5[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
	}
	
	gfft_r2c_t(wr10);
	gfft_r2c_t(wr11);
#ifndef WITH_2D
	gfft_r2c_t(wr12);
#endif
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		fldi.vx[i] +=  I * mask[i] / (nu * k2t) * (1 - exp(-nu*k2t[IDX3D]*dt)) * (
					kxt[i] * w10[i] + ky[i] * w7[i] + kz[i] * w8[i] );
		fldi.vy[i] +=  I * mask[i] / (nu * k2t) * (1 - exp(-nu*k2t[IDX3D]*dt)) * (
					kxt[i] * w7[i] + ky[i] * w11[i] + kz[i] * w9[i] );
		fldi.vz[i] +=  I * mask[i] / (nu * k2t) * (1 - exp(-nu*k2t[IDX3D]*dt)) * (
					kxt[i] * w8[i] + ky[i] * w9[i] + kz[i] * w12[i] );	// since kz=0 in 2D, kz*w6 gives 0, even if w6 is some random array
	}

#endif

	projector(fldi.vx,fldi.vy,fldi.vz);	
	return;
}

#endif

/*****************************************************
** ABC FLOW forcing                       ************
** (Thanks to E. Rempel)                  ************
******************************************************/

#ifdef ABC_FORCING

void ABC_forcing(struct Field fldi,
			 double dt) {

  const double A=NTOTAL*(2*M_PI)*(2*M_PI)*nu, B=NTOTAL*(2*M_PI)*(2*M_PI)*nu, C=NTOTAL*(2*M_PI)*(2*M_PI)*nu;

  const double kf=1.0;
  int i, j, k;
  int divisor,quociente,resto;

	
// ABC forcing

// Na direcao z nao precisa setar valores negativos de k, pois espelhamento e' feito automaticamente       
	if ((rank==0)){
	i=0;
        j=0;
        k=1;
	fldi.vx[IDX3D] += -0.5*A*I * dt;
	fldi.vy[IDX3D] += 0.5*A * dt;
	i=0;
	j=1;
        k=0;
	fldi.vx[IDX3D] += 0.5*C * dt;
	fldi.vz[IDX3D] += -0.5*C*I * dt;
	j=NY_COMPLEX - 1;
	fldi.vx[IDX3D] += 0.5*C * dt;
	fldi.vz[IDX3D] += 0.5*C*I * dt;
	}

	divisor=NX_COMPLEX/NPROC;
	quociente=(int)1/divisor;
	resto=fmod(1,divisor);

	// Somente kx eh paralelizado:
	if (rank==quociente){
	i=resto;
	j=0;
	k=0;
	fldi.vy[IDX3D] += -0.5*B*I * dt;
	fldi.vz[IDX3D] += 0.5*B * dt;
	}
	if (rank==(NPROC-1-quociente)){
	i = NX_COMPLEX/NPROC-resto;
	j=0;
	k=0;
	fldi.vy[IDX3D] += 0.5*B*I * dt;
	fldi.vz[IDX3D] += 0.5*B * dt;
	}
	
	projector(fldi.vx,fldi.vy,fldi.vz);
	
	return;
}


#endif

#ifdef IMPLICIT_ABC_FORCING

// Snoopy handles forcing in the same section of code in which
// the implicit dissipation is handled. Therefore it seems consistent
// to implement the forcing implicitly.

void implicit_ABC_forcing(struct Field fldi, double dt){
	int i,j,k;	// Used to address the fourier modes via the macro IDX3D

	// The amplitude is chosen so as to give us velocities matching Podvigina & Pouquet 1994, A=B=C=1

	// F_x += sin(z), F_y += cos(z)
	if(rank == 0){
		i=j=0; k=param.ABC_flow_kz;
		fldi.vx[IDX3D] += ((double)NTOTAL)*0.5*param.ABC_flow_A*(-I)*(1 - exp(-nu*k2t[IDX3D]*dt)); 
		fldi.vy[IDX3D] += ((double)NTOTAL)*0.5*param.ABC_flow_A*(+1)*(1 - exp(-nu*k2t[IDX3D]*dt)); 
	}
	// F_x += cos(y), F_z += sin(y)
	if(rank ==0 ){
		i=k=0; j=param.ABC_flow_ky;
		fldi.vx[IDX3D] += ((double)NTOTAL)*0.5*param.ABC_flow_C*(+1)*(1 - exp(-nu*k2t[IDX3D]*dt));
		fldi.vz[IDX3D] += ((double)NTOTAL)*0.5*param.ABC_flow_C*(-I)*(1 - exp(-nu*k2t[IDX3D]*dt));
		j = NY_COMPLEX - param.ABC_flow_ky;
		fldi.vx[IDX3D] += ((double)NTOTAL)*0.5*param.ABC_flow_C*(+1)*(1 - exp(-nu*k2t[IDX3D]*dt));
		fldi.vz[IDX3D] += ((double)NTOTAL)*0.5*param.ABC_flow_C*(+I)*(1 - exp(-nu*k2t[IDX3D]*dt));
	}
	// F_y += sin(x), F_z += cos(x)
	if(rank == (param.ABC_flow_kx/(NX_COMPLEX/NPROC))){
		i= (param.ABC_flow_kx%(NX_COMPLEX/NPROC)); j=k=0;
		fldi.vy[IDX3D] += ((double)NTOTAL)*0.5*param.ABC_flow_B*(-I)*(1 - exp(-nu*k2t[IDX3D]*dt));
		fldi.vz[IDX3D] += ((double)NTOTAL)*0.5*param.ABC_flow_B*(+1)*(1 - exp(-nu*k2t[IDX3D]*dt));
	}
	if(rank == ((NX_COMPLEX - param.ABC_flow_kx)/(NX_COMPLEX/NPROC))){
		i = ((NX_COMPLEX - param.ABC_flow_kx)%(NX_COMPLEX/NPROC)); j=k=0;
		fldi.vy[IDX3D] += ((double)NTOTAL)*0.5*param.ABC_flow_B*(+I)*(1 - exp(-nu*k2t[IDX3D]*dt));
		fldi.vz[IDX3D] += ((double)NTOTAL)*0.5*param.ABC_flow_B*(+1)*(1 - exp(-nu*k2t[IDX3D]*dt));
	}
	return;
}

#endif
