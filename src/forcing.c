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
