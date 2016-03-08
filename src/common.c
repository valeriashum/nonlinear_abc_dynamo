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



#include "snoopy.h"
#include "common.h"

#include <stdlib.h>

// Timing support
#ifndef MPI_SUPPORT
#ifndef _OPENMP
#include <stdio.h>
#include <time.h>
#endif
#endif

#include "error.h"
#include "gfft.h"
#include "debug.h"

#define MAX_N_BIN					10000
#define OUTPUT_SPECTRUM_K_BIN		(2.0 * M_PI/100.0)

// This are global variables used throughout the code
// Wave number pointers
double	*kx;	/**< x Wavevector */
double	*ky;	/**< y Wavevector */
double	*kz;	/**< z Wavevector */
double	*kxt;	/**< Time dependant x Wavevector. Different from kx only when SHEAR is present.*/
double	*k2t;	/**<  k squared Wavevector, function of time when SHEAR is present.*/
double	*ik2t;  /**< inverse of k2t Wavevector, function of time when SHEAR is present. set to 0 wheh k2t=0 to avoid singularity*/
double	kxmax,	kymax,  kzmax,	kmax;	/**< Maximum wavevectors */

#ifdef OSCILLATORY_SHEAR
// Fourier-integrated components of the Reynolds stress.
// Need to be global so that output() can see it. Bit of a kludge to be honest!
double integrated_cos_Rxy, integrated_sin_Rxy;
// The two variables below are to allow me to use RK3 integration
double integrated_cos_Rxy1, integrated_sin_Rxy1;
#ifdef SHELL_R_XY
double integrated_sin_Rxy_shell[ MAX_N_BIN ];
double integrated_cos_Rxy_shell[ MAX_N_BIN ];
double integrated_sin_Rxy_shell1[ MAX_N_BIN ];
double integrated_cos_Rxy_shell1[ MAX_N_BIN ];
#endif
#endif


// Mask for dealiasing
double   *mask;	/**< Deasliasing Mask*/

double	*wr1,	*wr2,	*wr3;		/** Temporary real array (alias of complex w**) */
double	*wr4,	*wr5,	*wr6;		/** Temporary real array (alias of complex w**) */
double	*wr7,	*wr8,	*wr9;		/** Temporary real array (alias of complex w**) */
double	*wr10,	*wr11,	*wr12;
double	*wr13,	*wr14,	*wr15;
double  *wr16,  *wr17;

double complex		*w1,	*w2,	*w3;	/**< Temporary complex array (alias of real wr**) */
double complex		*w4,	*w5,	*w6;	/**< Temporary complex array (alias of real wr**) */
double complex		*w7,	*w8,	*w9;	/**< Temporary complex array (alias of real wr**) */
double complex		*w10,	*w11,	*w12;
double complex		*w13,	*w14,	*w15;
double complex		*w16,   *w17;


// Parameters
struct Parameters			param;

// Physics variables 
double	nu;
#ifdef BOUSSINESQ
double	nu_th;
#endif
#ifdef MHD
double	eta;
double  em_prev;
double  em_100_prev;
double  em_110_prev;
double  em_111_prev;
double  em_k00_prev;
double  em_kk0_prev;
double  em_kkk_prev;
double  em_gr;
#endif

#ifdef MPI_SUPPORT
int		NPROC;									/**< NPROC is a variable when MPI is on. Otherwise, it is preprocessor macro in gvars.h */
#endif

int		rank;
int		nthreads;								/**< Number of OpenMP threads */


/* Function prototypes */
void allocate_field(struct Field *fldi);
void deallocate_field(struct Field *fldi);
void init_N2_profile();
void init_real_mask();

/***************************************************************/
/**
	Init all global variables, aligning them in memory
*/
/***************************************************************/
void init_common(void) {
	/* This routine will initialize everything */
	int i,j,k;
	
	DEBUG_START_FUNC;
 printf("rank=%d\n", rank);	
#ifdef MPI_SUPPORT
#ifdef FFTW3_MPI_SUPPORT	
	fftw_mpi_init();
#endif
#endif
#ifdef _OPENMP
	if( !(fftw_init_threads()) ) ERROR_HANDLER( ERROR_CRITICAL, "Threads initialisation failed");
#endif
	
	/* We start with the coordinate system */
	kx = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (kx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for kx allocation");
	
	ky = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (ky == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for ky allocation");
	
	kz = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (kz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for kz allocation");
	
	kxt = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (kxt == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for kxt allocation");
	
	k2t = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (k2t == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for k2t allocation");
	
	ik2t = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (ik2t == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for ik2t allocation");


	for( i = 0; i < NX_COMPLEX / NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				kx[ IDX3D ] = (2.0 * M_PI) / param.lx *
						(fmod( NX_COMPLEX * rank / NPROC  + i + (NX_COMPLEX / 2) ,  NX_COMPLEX ) - NX_COMPLEX / 2 );
#ifdef WITH_2D
				ky[ IDX3D ] = (2.0 * M_PI) / param.ly * j;
					 
				kz[ IDX3D ] = 0.0;
#else	
				ky[ IDX3D ] = (2.0 * M_PI) / param.ly *
						(fmod( j + (NY_COMPLEX / 2) ,  NY_COMPLEX ) - NY_COMPLEX / 2 );
				//if( j < 5 && i<5 && k<5) MPI_Printf("From common.c: for i = %i, j =%i, k=%i, ky = %f\n", i,j,k, ky[IDX3D]);
	 
				kz[ IDX3D ] = (2.0 * M_PI) / param.lz * k;
#endif

				kxt[ IDX3D ]= kx[IDX3D];
			
				k2t[ IDX3D ] = kxt[IDX3D] * kxt[IDX3D] +
								ky[IDX3D] * ky[IDX3D] +
								kz[IDX3D] * kz[IDX3D];
							  
				if ( k2t[IDX3D] == 0.0 ) ik2t[IDX3D] = 1.0;
				else	ik2t[IDX3D] = 1.0 / k2t[IDX3D];
			}
		}
	}
	
	kxmax = 2.0 * M_PI/ param.lx * ( (NX / 2) - 1);
	kymax = 2.0 * M_PI/ param.ly * ( (NY / 2) - 1);
	kzmax = 2.0 * M_PI/ param.lz * ( (NZ / 2) - 1);

				/*	MPI_Printf("From common.c: check kxmax = %f\n", kxmax);
					MPI_Printf("From common.c: check kymax = %f\n", kymax);
					MPI_Printf("From common.c: check kzmax = %f\n", kzmax);
					MPI_Printf("From common.c: lx = %f, ly =  %f,lz =  %f\n ", param.lx, param.ly,param.lz);*/

		
#ifdef WITH_2D
	kzmax = 0.0;
#endif
	
	kmax=pow(kxmax*kxmax+kymax*kymax+kzmax*kzmax,0.5);


	/* Initialize the dealiazing mask Or the nyquist frequency mask (in case dealiasing is not required) */
	
	mask = (double *) fftw_malloc( sizeof(double) * NTOTAL_COMPLEX );
	if (mask == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for mask allocation");
	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {

				mask[ IDX3D ] = 1.0;
				if(param.antialiasing) {
					if( fabs( kx[ IDX3D ] ) > 2.0/3.0 * kxmax)
						mask[ IDX3D ] = 0.0;
                    				
					if( fabs( ky[ IDX3D ] ) > 2.0/3.0 * kymax)
						mask[ IDX3D ] = 0.0;
#ifndef WITH_2D
					if( fabs( kz[ IDX3D ] ) > 2.0/3.0 * kzmax)
						mask[ IDX3D ] = 0.0;
#endif
				}
				else {			
					if (  NX_COMPLEX / NPROC * rank + i == NX_COMPLEX / 2 ) 
						mask[ IDX3D ] = 0.0;
					if ( j == NY_COMPLEX / 2 )  
						mask[ IDX3D ] = 0.0;
#ifndef WITH_2D
					if ( k == NZ_COMPLEX ) 
						mask[ IDX3D ] = 0.0;
#endif
				}
			}
		}
	}

	if(param.antialiasing) {
		kxmax = kxmax * 2.0 / 3.0;
		kymax = kymax * 2.0 / 3.0;
		kzmax = kzmax * 2.0 / 3.0;
		kmax = kmax * 2.0 / 3.0;
	}
	

// Allocate fields
// Complex fields

	w1 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w1 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w1 allocation");
	
	w2 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w2 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w2 allocation");
	
	w3 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w3 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w3 allocation");
	
	w4 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w4 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w4 allocation");
	
	w5 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w5 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w5 allocation");
	
	w6 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w6 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w6 allocation");
	
	w7 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w7 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w7 allocation");
	
	w8 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w8 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w8 allocation");
	
	w9 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w9 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w9 allocation");
	
	w10 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w10 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w10 allocation");
	
	w11 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w11 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w11 allocation");
	
	w12 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w12 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w12 allocation");
	
	w13 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w13 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w13 allocation");
	
	w14 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w14 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w14 allocation");
	
	w15 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w15 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w15 allocation");
		
    w16 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w16 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w16 allocation");

    w17 =(double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (w17 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w17 allocation");


	/* Will use the same memory space for real and complex fields */
	
	wr1 = (double *) w1;
	wr2 = (double *) w2;
	wr3 = (double *) w3;
	wr4 = (double *) w4;
	wr5 = (double *) w5;
	wr6 = (double *) w6;
	wr7 = (double *) w7;
	wr8 = (double *) w8;
	wr9 = (double *) w9;
	wr10 = (double *) w10;
	wr11 = (double *) w11;
	wr12 = (double *) w12;
	wr13 = (double *) w13;
	wr14 = (double *) w14;
	wr15 = (double *) w15;
    wr16 = (double *) w16;
    wr17 = (double *) w17;
	
// Physic initialisation
//	init_real_mask();
	
	nu = 1.0 / param.reynolds;
#ifdef BOUSSINESQ	
	nu_th = 1.0 / param.reynolds_th;
#endif
#ifdef MHD
    //eta = param.etta;
	eta = 1.0 / param.reynolds_m;
    em_prev = 0.0;
    em_100_prev = 0.0;
    em_110_prev = 0.0;
    em_111_prev = 0.0;
    em_k00_prev = 0.0;
    em_kk0_prev = 0.0;
    em_kkk_prev = 0.0;
    em_gr = 10.0;
#endif
	DEBUG_END_FUNC;
	return;
}

void finish_common(void) {
	free(kx);
	free(ky);
	free(kz);
	free(kxt);
	free(k2t);
	free(ik2t);
	free(mask);
	
	free(w1);
	free(w2);
	free(w3);
	free(w4);
	free(w5);
	free(w6);
	free(w7);
	free(w8);
	free(w9);
	free(w10);
	free(w11);
	free(w12);
	free(w13);
	free(w14);
	free(w15);
    free(w16);
    free(w17);

	return;
}


/*********************************************/
/**
	Allocate a field structure according to the code
	current configuration
	This routine allows one to add extra fields
	to the code very easily.

	@param *fldi: pointer to a field structure to initialize
**/
/*********************************************/
void allocate_field(struct Field *fldi) {
	int current_field, i;
	
	DEBUG_START_FUNC;
	
	// We want to allocate a field structure
	fldi->nfield = 3;	
		
#ifdef BOUSSINESQ
	fldi->nfield++;
#endif

#ifdef MHD
	fldi->nfield=fldi->nfield+3;
#endif

#ifdef COMPRESSIBLE
	fldi->nfield=fldi->nfield+1;
#endif

	// Now we want to initialize the pointers of the field structure
	
	// farray will point to each of the array previously allocated
	fldi->farray = (double complex **) fftw_malloc( sizeof(double complex *) * fldi->nfield);
	if (fldi->farray == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fldi->farray allocation");
	
	// fname will point to the name of each field
	fldi->fname = (char **) fftw_malloc(sizeof(char *) * fldi->nfield);
	if (fldi->fname == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fldi->fname allocation");
	
	// Initialise the pointers
	for(i=0 ; i < fldi->nfield ; i++) {
		fldi->fname[i] = (char *) fftw_malloc(sizeof(char) * 10); // 10 character to describe each field
		if (fldi->fname[i] == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fldi->fname[i] allocation");
	}
	
	// Allocate the arrays and put the right value in each pointer
	
	current_field = 0;

	fldi->vx = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fldi->vx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fldi->vx allocation");
#ifdef WITH_OL2012B
	fldi->vdotx = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fldi->vdotx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fldi->vdotx allocation");
#endif
	fldi->farray[current_field] = fldi->vx;
#ifndef COMPRESSIBLE
	sprintf(fldi->fname[current_field],"vx");
#else
	sprintf(fldi->fname[current_field],"px");
#endif
	current_field++;
	
	fldi->vy = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fldi->vy == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fldi->vy allocation");
#ifdef WITH_OL2012B
	fldi->vdoty = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fldi->vdoty == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fldi->vdoty allocation");
#endif
	fldi->farray[current_field] = fldi->vy;
#ifndef COMPRESSIBLE
	sprintf(fldi->fname[current_field],"vy");
#else
	sprintf(fldi->fname[current_field],"py");
#endif
	current_field++;
	
	fldi->vz = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fldi->vz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fldi->vz allocation");
#ifdef WITH_OL2012B
	fldi->vdotz = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fldi->vdotz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fldi->vdotz allocation");
#endif
	fldi->farray[current_field] = fldi->vz;
#ifndef COMPRESSIBLE
	sprintf(fldi->fname[current_field],"vz");
#else
	sprintf(fldi->fname[current_field],"pz");
#endif
	current_field++;
	
#ifdef BOUSSINESQ
	fldi->th = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fldi->th == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fldi->th allocation");
	fldi->farray[current_field] = fldi->th;
	sprintf(fldi->fname[current_field],"th");
	current_field++;
#endif
#ifdef MHD
	fldi->bx = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fldi->bx == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fldi->bx allocation");
#ifdef KINEMATIC_GROWTH
    fldi->bx0 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fldi->bx0 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fldi->bx0 allocation");
#endif
	fldi->farray[current_field] = fldi->bx;
	sprintf(fldi->fname[current_field],"bx");
	current_field++;
	
	fldi->by = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fldi->by == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fldi->by allocation");
#ifdef KINEMATIC_GROWTH
    fldi->by0 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fldi->by0 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fldi->by0 allocation");
#endif
	fldi->farray[current_field] = fldi->by;
	sprintf(fldi->fname[current_field],"by");
	current_field++;
	
	fldi->bz = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fldi->bz == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fldi->bz allocation");
#ifdef KINEMATIC_GROWTH
    fldi->bz0 = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fldi->bz0 == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fldi->bz0 allocation");
#endif
	fldi->farray[current_field] = fldi->bz;
	sprintf(fldi->fname[current_field],"bz");
	current_field++;
#endif
#ifdef COMPRESSIBLE
	fldi->d = (double complex *) fftw_malloc( sizeof(double complex) * NTOTAL_COMPLEX);
	if (fldi->d == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for fldi->d allocation");
	fldi->farray[current_field] = fldi->d;
	sprintf(fldi->fname[current_field]," d");
	current_field++;
#endif
	// Add a field here if you need one... (don't forget to ajust fldi.nfield accordingly)
	// *


#ifdef WITH_PARTICLES	
	// Init space for particle storage
	fldi->part = (struct Particle *) malloc(sizeof(struct Particle) * param.particles_n);
	
#endif
	// Ok, all done...
	
	DEBUG_END_FUNC;
	
	return;
}
/*********************************************/
/**
Deallocate a field structure created by
allocate_field
**/
/*********************************************/
void deallocate_field(struct Field *fldi) {
	int i;
	// Free a field structure
	
	DEBUG_START_FUNC;
	
	for(i=0 ; i < fldi->nfield ; i++) {
		fftw_free(fldi->fname[i]);
		fftw_free(fldi->farray[i]);
	}
	
	fftw_free(fldi->farray);
	fftw_free(fldi->fname);
	
#ifdef WITH_PARTICLES
	free(fldi->part);
#endif
	
	// Done
	
	DEBUG_END_FUNC;
	
	return;
}
	

/*********************************************/
/**
Customized random number generator
Allow one to have consistant random numbers
generators on different architectures.
**/
/*********************************************/
double randm(void) {
	const int a = 16807;
	const int m = 2147483647;
	static int in0 = 13763;
	int q;

#ifdef RAND_OFFSET
	// When using mpi, this allows us to have different number series in each process...
	if(in0 == 13763) in0 += 2543 * (rank + RAND_OFFSET);
#else
	// When using mpi, this allows us to have different number series in each process...
	if(in0 == 13763) in0 += 2543 * rank;
#endif
	
	/* find random number  */
	q= (int) fmod((double) a * in0, m);
	in0=q;
	
	return((double)q/(double)m);
}
/*********************************************/
/**
	 * Normal distribution
	 * Algorithm by D.E. Knut, 1997, The Art of Computer Programmin, Addison-Wesley. 
	 */
/*********************************************/
	 
double randm_normal(void) {
	double v1, v2;
	double rsq=1.0;
	
	while(rsq>=1. || rsq==0.0) {
		v1=2.*randm()-1.0;
		v2=2.*randm()-1.0;
		rsq=v1*v1+v2*v2;
	}
	
	return( v1*sqrt(-2.0 * log(rsq) / rsq));
}

/****************************************************/
/**
	Remove the divergence from a 3D field using
	the projector operator:
	q=q-k.q/k^2
	
	@param qx: x component of the field
	@param qy: y component of the field
	@param qz: z component of the field
*/
/****************************************************/
	
void projector( double complex qx[],
			    double complex qy[],
			    double complex qz[]) {
				
	int i;
	double complex q0;
	
	DEBUG_START_FUNC;
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		q0 = kxt[i] * qx[i] + ky[i] * qy[i] + kz[i] * qz[i];
		qx[i] = qx[i] - kxt[i] * q0 * ik2t[i];
		qy[i] = qy[i] - ky[i] * q0 * ik2t[i];
		qz[i] = qz[i] - kz[i] * q0 * ik2t[i];
	}
	
	DEBUG_END_FUNC;
	
	return;
}

/*********************************************/
/** Compute the energy of a given field.
	@param q complex array containing the field for which
				  we want the total energy
*/
/*********************************************/

double energy(const double complex q[]) {
	
	int i,j,k;
	double energ_tot;
	
	energ_tot=0.0;
	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k=0; k < NZ_COMPLEX; k++) {
#ifdef WITH_2D
				if( j == 0)
#else
				if( k == 0)
#endif
					// k=0, we have all the modes.
					energ_tot = energ_tot + creal( 0.5 * q[ IDX3D ] * conj( q[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
				else
					// k>0, only half of the complex plane is represented.
					energ_tot = energ_tot + creal( q[ IDX3D ] * conj( q[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
			}
		}
	}
//	energ_tot = 0;
	return(energ_tot);
}

/*********************************************/
/* Compute the 'dissipation' of a particular
	field, i.e. the sum over wavenumbers of
	- visc * k^2 * |q_k|^2
	Based on Geoffroy's energy() function above.
*/
/*********************************************/
double dissipation_rate(const double complex q[], const double visc){
	
	int i,j,k;
	double disp_tot;
	
	disp_tot=0.0;
	
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k=0; k < NZ_COMPLEX; k++) {
#ifdef WITH_2D
				if( j == 0)
#else
				if( k == 0)
#endif
					// k=0, we have all the modes.
					disp_tot = disp_tot + visc * k2t[IDX3D] * creal( q[ IDX3D ] * conj( q[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
				else
					// k>0, only half of the complex plane is represented.
					disp_tot = disp_tot + visc * k2t[IDX3D] * creal( 2.0 * q[ IDX3D ] * conj( q[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
			}
		}
	}
	return(disp_tot);
}


/*********************************************/
/*  Compute the A_1jj1 - 2A_1221 - 2A_2121 tensor 
    of equation (31) of Ogilvie & Lesur (2012)
    Rather, it computes the contribution to this quantity
    from the wavenumbers held by this process
*/
/*********************************************/
double OL2012A(const double complex vx[], const double complex vy[]){

    int i,j,k;
	double A_output;

	A_output=0;

	for(i=0; i < NX_COMPLEX/NPROC; i++){
		for(j=0; j < NY_COMPLEX; j++){
			for(k=0; k < NZ_COMPLEX; k++){
				if(!(i==0 && j==0 && k==0 && rank==0)){
#ifdef WITH_2D
					if(j==0)
#else
					if(k==0)
#endif
						// k=0, we have all of the modes
						A_output = A_output + creal((conj(vx[IDX3D])*(kxt[IDX3D]*kxt[IDX3D] - ky[IDX3D]*ky[IDX3D] + kz[IDX3D]*kz[IDX3D])
						                             -2*conj(vy[IDX3D])*kx[IDX3D]*ky[IDX3D])*vx[IDX3D] / (((double) NTOTAL*NTOTAL)*k2t[IDX3D]));
					else
						// k>0, one half of the complex plane is represented
						A_output = A_output + 2.0*creal((conj(vx[IDX3D])*(kxt[IDX3D]*kxt[IDX3D] - ky[IDX3D]*ky[IDX3D] + kz[IDX3D]*kz[IDX3D])
						                                -2*conj(vy[IDX3D])*kx[IDX3D]*ky[IDX3D])*vx[IDX3D] / (((double) NTOTAL*NTOTAL)*k2t[IDX3D]));
				}
			}
		}
	}
	return A_output;
}


/*********************************************/
/*  Compute and return the tensor at the top of
    page 5 of OL2012 that controls the dissipation
*/
/*********************************************/


#ifdef WITH_OL2012B
double OL2012B(const struct Field fldi){
	
	int i,j,k;
	double B_output, C_output, D_output, E_output;
	
	// Compute B_1jj1 - B_1221 - B_1122
	B_output=0;
	for(i=0; i < NX_COMPLEX/NPROC; i++){
		for(j=0; j < NY_COMPLEX; j++){
			for(k=0; k < NZ_COMPLEX; k++){
				if(!(i==0 && j==0 && k==0 && rank==0)){
#ifdef WITH_2D
					if(j==0)
#else
					if(k==0)
#endif
						// k=0, we have all of the modes
						B_output = B_output + creal((conj(fldi.vdotx[IDX3D])*((kxt[IDX3D]*kxt[IDX3D] + kz[IDX3D]*kz[IDX3D])*fldi.vx[IDX3D]
						                                                       - kxt[IDX3D]*ky[IDX3D]*fldi.vy[IDX3D])) / (((double) NTOTAL*NTOTAL)*k2t[IDX3D]));
					else
						// k>0, one half of the complex plane is represented
						B_output = B_output + 2.0*creal((conj(fldi.vdotx[IDX3D])*((kxt[IDX3D]*kxt[IDX3D] + kz[IDX3D]*kz[IDX3D])*fldi.vx[IDX3D]
						                                                           - kxt[IDX3D]*ky[IDX3D]*fldi.vy[IDX3D])) / (((double) NTOTAL*NTOTAL)*k2t[IDX3D]));
				}
			}
		}
	}

	// Compute -C_1221 + C_1jj1 + 3C_1122
	C_output=0;
	for(i=0; i < NX_COMPLEX/NPROC; i++){
		for(j=0; j < NY_COMPLEX; j++){
			for(k=0; k < NZ_COMPLEX; k++){
				if(!(i==0 && j==0 && k==0 && rank==0)){
#ifdef WITH_2D
					if(j==0)
#else
					if(k==0)
#endif
						// k=0, we have all of the modes
						C_output = C_output + (-nu)*(-1.0)*creal((conj(fldi.vx[IDX3D])*((kxt[IDX3D]*kxt[IDX3D] + kz[IDX3D]*kz[IDX3D])*fldi.vx[IDX3D]
						                                                         + 3.0*kxt[IDX3D]*ky[IDX3D]*fldi.vy[IDX3D])) / ((double) NTOTAL*NTOTAL));
					else
						// k>0, one half of the complex plane is represented
						C_output = C_output + 2.0*(-nu)*(-1.0)*creal((conj(fldi.vx[IDX3D])*((kxt[IDX3D]*kxt[IDX3D] + kz[IDX3D]*kz[IDX3D])*fldi.vx[IDX3D]
						                                                             + 3.0*kxt[IDX3D]*ky[IDX3D]*fldi.vy[IDX3D])) / ((double) NTOTAL*NTOTAL));
				}
			}
		}
	}

	D_output=0;
	/* Compute v_i v_j by convolution (lifted from timestep.c) */
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fldi.vx[i];
		w2[i] =  fldi.vy[i];
		w3[i] =  fldi.vz[i];
	}

	projector(w1,w2,w3);
	
	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr4[i] = wr1[i] * wr1[i] / ((double) NTOTAL*NTOTAL);
		wr5[i] = wr2[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
#ifndef WITH_2D
		wr6[i] = wr3[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
#endif
		wr7[i] = wr1[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr8[i] = wr1[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr9[i] = wr2[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
	}
	
	gfft_r2c_t(wr4);
	gfft_r2c_t(wr5);
#ifndef WITH_2D
	gfft_r2c_t(wr6);
#endif
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);

	/* Done with Convolution for v_i v_j computation */

// Commented out since D_1jj1 vanishes
//	// Compute 2 * D_1jj1 = 2 * < v_1 v_j d_j v_1 >
//	for(i=0; i < NX_COMPLEX/NPROC; i++){
//		for(j=0; j < NY_COMPLEX; j++){
//			for(k=0; k < NZ_COMPLEX; k++){
//				if(!(i==0 && j==0 && k==0 && rank==0)){
//#ifdef WITH_2D
//					if(j==0)
//#else
//					if(k==0)
//#endif
//						// k=0, we have all of the modes
//						D_output = D_output + 2.0*creal((conj(w4[IDX3D])*I*kxt[IDX3D] + conj(w7[IDX3D])*I*ky[IDX3D]
//						                                 + conj(w8[IDX3D])*I*kz[IDX3D])*fldi.vx[IDX3D] / ((double) NTOTAL*NTOTAL));
//					else
//						// k>0, one half of the complex plane is represented
//						D_output = D_output + 2.0*2.0*creal((conj(w4[IDX3D])*I*kxt[IDX3D] + conj(w7[IDX3D])*I*ky[IDX3D]
//						                                     + conj(w8[IDX3D])*I*kz[IDX3D])*fldi.vx[IDX3D] / ((double) NTOTAL*NTOTAL));
//				}
//			}
//		}
//	}

	// Compute -2 * D_1jj221 = -2 * < v_1 v_j d_j d_2 d_2 d^-2 v_1>
	for(i=0; i < NX_COMPLEX/NPROC; i++){
		for(j=0; j < NY_COMPLEX; j++){
			for(k=0; k < NZ_COMPLEX; k++){
				if(!(i==0 && j==0 && k==0 && rank==0)){
#ifdef WITH_2D
					if(j==0)
#else
					if(k==0)
#endif
						// k=0, we have all of the modes
						D_output = D_output + -2.0*creal((conj(w4[IDX3D])*I*kxt[IDX3D]*ky[IDX3D]*ky[IDX3D] + conj(w7[IDX3D])*I*ky[IDX3D]*ky[IDX3D]*ky[IDX3D]
						                                  + conj(w8[IDX3D])*I*kz[IDX3D]*ky[IDX3D]*ky[IDX3D])*(fldi.vx[IDX3D]/k2t[IDX3D]) / ((double) NTOTAL*NTOTAL));
					else
						// k>0, one half of the complex plane is represented
						D_output = D_output + 2.0*(-2.0)*creal((conj(w4[IDX3D])*I*kxt[IDX3D]*ky[IDX3D]*ky[IDX3D] + conj(w7[IDX3D])*I*ky[IDX3D]*ky[IDX3D]*ky[IDX3D]
						                                        + conj(w8[IDX3D])*I*kz[IDX3D]*ky[IDX3D]*ky[IDX3D])*(fldi.vx[IDX3D]/k2t[IDX3D]) / ((double) NTOTAL*NTOTAL));
				}
			}
		}
	}

	// Compute -2 * D_2jj121 = -2 * < v_2 v_j d_j d_1 d_2 d^-2 v_1>
	for(i=0; i < NX_COMPLEX/NPROC; i++){
		for(j=0; j < NY_COMPLEX; j++){
			for(k=0; k < NZ_COMPLEX; k++){
				if(!(i==0 && j==0 && k==0 && rank==0)){
#ifdef WITH_2D
					if(j==0)
#else
					if(k==0)
#endif
						// k=0, we have all of the modes
						D_output = D_output + -2.0*creal((conj(w7[IDX3D])*I*kxt[IDX3D]*kxt[IDX3D]*ky[IDX3D] + conj(w5[IDX3D])*I*ky[IDX3D]*kxt[IDX3D]*ky[IDX3D]
						                                  + conj(w9[IDX3D])*I*kz[IDX3D]*kxt[IDX3D]*ky[IDX3D])*(fldi.vx[IDX3D]/k2t[IDX3D]) / ((double) NTOTAL*NTOTAL));
					else
						// k>0, one half of the complex plane is represented
						D_output = D_output + 2.0*(-2.0)*creal((conj(w7[IDX3D])*I*kxt[IDX3D]*kxt[IDX3D]*ky[IDX3D] + conj(w5[IDX3D])*I*ky[IDX3D]*kxt[IDX3D]*ky[IDX3D]
						                                        + conj(w9[IDX3D])*I*kz[IDX3D]*kxt[IDX3D]*ky[IDX3D])*(fldi.vx[IDX3D]/k2t[IDX3D]) / ((double) NTOTAL*NTOTAL));
				}
			}
		}
	}

	// Compute +3 * D_1jj221 = +3 * < v_1 v_j d_j d_2 d_2 d^-2 v_1 >
	for(i=0; i < NX_COMPLEX/NPROC; i++){
		for(j=0; j < NY_COMPLEX; j++){
			for(k=0; k < NZ_COMPLEX; k++){
				if(!(i==0 && j==0 && k==0 && rank==0)){
#ifdef WITH_2D
					if(j==0)
#else
					if(k==0)
#endif
						// k=0, we have all of the modes
						D_output = D_output + 3.0*creal((conj(w4[IDX3D])*I*kxt[IDX3D]*ky[IDX3D]*ky[IDX3D] + conj(w7[IDX3D])*I*ky[IDX3D]*ky[IDX3D]*ky[IDX3D]
						                                 + conj(w8[IDX3D])*I*kz[IDX3D]*ky[IDX3D]*ky[IDX3D])*(fldi.vx[IDX3D]/k2t[IDX3D]) / ((double) NTOTAL*NTOTAL));
					else
						// k>0, one half of the complex plane is represented
						D_output = D_output + 2.0*(3.0)*creal((conj(w4[IDX3D])*I*kxt[IDX3D]*ky[IDX3D]*ky[IDX3D] + conj(w7[IDX3D])*I*ky[IDX3D]*ky[IDX3D]*ky[IDX3D]
						                                       + conj(w8[IDX3D])*I*kz[IDX3D]*ky[IDX3D]*ky[IDX3D])*(fldi.vx[IDX3D]/k2t[IDX3D]) / ((double) NTOTAL*NTOTAL));
				}
			}
		}
	}

	// Compute +3 * D_1jj212 = +3 * < v_1 v_j d_j d_2 d_1 d^-2 v_2 >
	for(i=0; i < NX_COMPLEX/NPROC; i++){
		for(j=0; j < NY_COMPLEX; j++){
			for(k=0; k < NZ_COMPLEX; k++){
				if(!(i==0 && j==0 && k==0 && rank==0)){
#ifdef WITH_2D
					if(j==0)
#else
					if(k==0)
#endif
						// k=0, we have all of the modes
						D_output = D_output + 3.0*creal((conj(w4[IDX3D])*I*kxt[IDX3D]*ky[IDX3D]*kxt[IDX3D] + conj(w7[IDX3D])*I*ky[IDX3D]*ky[IDX3D]*kxt[IDX3D]
						                                 + conj(w8[IDX3D])*I*kz[IDX3D]*ky[IDX3D]*kxt[IDX3D])*(fldi.vy[IDX3D]/k2t[IDX3D]) / ((double) NTOTAL*NTOTAL));
					else
						// k>0, one half of the complex plane is represented
						D_output = D_output + 2.0*(3.0)*creal((conj(w4[IDX3D])*I*kxt[IDX3D]*ky[IDX3D]*kxt[IDX3D] + conj(w7[IDX3D])*I*ky[IDX3D]*ky[IDX3D]*kxt[IDX3D]
						                                       + conj(w8[IDX3D])*I*kz[IDX3D]*ky[IDX3D]*kxt[IDX3D])*(fldi.vy[IDX3D]/k2t[IDX3D]) / ((double) NTOTAL*NTOTAL));
				}
			}
		}
	}

	// Compute -D_ijij1221 = - < v_i v_j d_i d_j d_1 d_2 d_2 d^-4 v_1 >
	for(i=0; i < NX_COMPLEX/NPROC; i++){
		for(j=0; j < NY_COMPLEX; j++){
			for(k=0; k < NZ_COMPLEX; k++){
				if(!(i==0 && j==0 && k==0 && rank==0)){
#ifdef WITH_2D
					if(j==0)
#else
					if(k==0)
#endif
						// k=0, we have all of the modes
						D_output = D_output + (-1.0)*creal((     conj(w4[IDX3D])*kxt[IDX3D]*kxt[IDX3D]
						                                    +    conj(w5[IDX3D])*ky[IDX3D] *ky[IDX3D]
						                                    +    conj(w6[IDX3D])*kz[IDX3D] *kz[IDX3D]
						                                    +2.0*conj(w7[IDX3D])*kxt[IDX3D]*ky[IDX3D]
						                                    +2.0*conj(w8[IDX3D])*kxt[IDX3D]*kz[IDX3D]
						                                    +2.0*conj(w9[IDX3D])*ky[IDX3D] *kz[IDX3D])
						                                   *(I*kxt[IDX3D]*ky[IDX3D]*ky[IDX3D])*(fldi.vx[IDX3D]/(k2t[IDX3D]*k2t[IDX3D])) / ((double) NTOTAL*NTOTAL));
					else
						// k>0, one half of the complex plane is represented
						D_output = D_output + 2.0*(-1.0)*creal((     conj(w4[IDX3D])*kxt[IDX3D]*kxt[IDX3D]
						                                        +    conj(w5[IDX3D])*ky[IDX3D] *ky[IDX3D]
						                                        +    conj(w6[IDX3D])*kz[IDX3D] *kz[IDX3D]
						                                        +2.0*conj(w7[IDX3D])*kxt[IDX3D]*ky[IDX3D]
						                                        +2.0*conj(w8[IDX3D])*kxt[IDX3D]*kz[IDX3D]
						                                        +2.0*conj(w9[IDX3D])*ky[IDX3D] *kz[IDX3D])
						                                       *(I*kxt[IDX3D]*ky[IDX3D]*ky[IDX3D])*(fldi.vx[IDX3D]/(k2t[IDX3D]*k2t[IDX3D])) / ((double) NTOTAL*NTOTAL));
				}
			}
		}
	}

	// Compute -D_ijij1212 = - < v_i v_j d_i d_j d_1 d_2 d_1 d^-4 v_2 >
	for(i=0; i < NX_COMPLEX/NPROC; i++){
		for(j=0; j < NY_COMPLEX; j++){
			for(k=0; k < NZ_COMPLEX; k++){
				if(!(i==0 && j==0 && k==0 && rank==0)){
#ifdef WITH_2D
					if(j==0)
#else
					if(k==0)
#endif
						// k=0, we have all of the modes
						D_output = D_output + (-1.0)*creal((     conj(w4[IDX3D])*kxt[IDX3D]*kxt[IDX3D]
						                                    +    conj(w5[IDX3D])*ky[IDX3D] *ky[IDX3D]
						                                    +    conj(w6[IDX3D])*kz[IDX3D] *kz[IDX3D]
						                                    +2.0*conj(w7[IDX3D])*kxt[IDX3D]*ky[IDX3D]
						                                    +2.0*conj(w8[IDX3D])*kxt[IDX3D]*kz[IDX3D]
						                                    +2.0*conj(w9[IDX3D])*ky[IDX3D] *kz[IDX3D])
						                                   *(I*kxt[IDX3D]*ky[IDX3D]*kxt[IDX3D])*(fldi.vy[IDX3D]/(k2t[IDX3D]*k2t[IDX3D])) / ((double) NTOTAL*NTOTAL));
					else
						// k>0, one half of the complex plane is represented
						D_output = D_output + 2.0*(-1.0)*creal((     conj(w4[IDX3D])*kxt[IDX3D]*kxt[IDX3D]
						                                        +    conj(w5[IDX3D])*ky[IDX3D] *ky[IDX3D]
						                                        +    conj(w6[IDX3D])*kz[IDX3D] *kz[IDX3D]
						                                        +2.0*conj(w7[IDX3D])*kxt[IDX3D]*ky[IDX3D]
						                                        +2.0*conj(w8[IDX3D])*kxt[IDX3D]*kz[IDX3D]
						                                        +2.0*conj(w9[IDX3D])*ky[IDX3D] *kz[IDX3D])
						                                       *(I*kxt[IDX3D]*ky[IDX3D]*kxt[IDX3D])*(fldi.vy[IDX3D]/(k2t[IDX3D]*k2t[IDX3D])) / ((double) NTOTAL*NTOTAL));
				}
			}
		}
	}

	// Compute D_ijij22 = < v_i v_j d_i d_j d_2 d^-2 v_2 >
	for(i=0; i < NX_COMPLEX/NPROC; i++){
		for(j=0; j < NY_COMPLEX; j++){
			for(k=0; k < NZ_COMPLEX; k++){
				if(!(i==0 && j==0 && k==0 && rank==0)){
#ifdef WITH_2D
					if(j==0)
#else
					if(k==0)
#endif
						// k=0, we have all of the modes
						D_output = D_output + (+1.0)*creal((     conj(w4[IDX3D])*kxt[IDX3D]*kxt[IDX3D]
						                                    +    conj(w5[IDX3D])*ky[IDX3D] *ky[IDX3D]
						                                    +    conj(w6[IDX3D])*kz[IDX3D] *kz[IDX3D]
						                                    +2.0*conj(w7[IDX3D])*kxt[IDX3D]*ky[IDX3D]
						                                    +2.0*conj(w8[IDX3D])*kxt[IDX3D]*kz[IDX3D]
						                                    +2.0*conj(w9[IDX3D])*ky[IDX3D] *kz[IDX3D])
						                                   *(+I*ky[IDX3D])*(fldi.vy[IDX3D]/k2t[IDX3D]) / ((double) NTOTAL*NTOTAL));
					else
						// k>0, one half of the complex plane is represented
						D_output = D_output + 2.0*(+1.0)*creal((     conj(w4[IDX3D])*kxt[IDX3D]*kxt[IDX3D]
						                                        +    conj(w5[IDX3D])*ky[IDX3D] *ky[IDX3D]
						                                        +    conj(w6[IDX3D])*kz[IDX3D] *kz[IDX3D]
						                                        +2.0*conj(w7[IDX3D])*kxt[IDX3D]*ky[IDX3D]
						                                        +2.0*conj(w8[IDX3D])*kxt[IDX3D]*kz[IDX3D]
						                                        +2.0*conj(w9[IDX3D])*ky[IDX3D] *kz[IDX3D])
						                                       *(+I*ky[IDX3D])*(fldi.vy[IDX3D]/k2t[IDX3D]) / ((double) NTOTAL*NTOTAL));
				}
			}
		}
	}
	
	/* Begin computation of 4*E_2121 + 4*E_2112 */
	E_output = 0.0;
	// E_ijkl = <v_m (d_m d_n d^-2 d_i v_j) (d_n d^-2 d_k v_l)>
	// The derivatives are not all on the same term, making this somewhat harder

	// Convolutions first
	for(i=0; i < NX_COMPLEX/NPROC; i++){
		for(j=0; j < NY_COMPLEX; j++){
			for(k=0; k < NZ_COMPLEX; k++){
				if(!(i==0 && j==0 && k==0 && rank==0)){
					w1[IDX3D] = I*kxt[IDX3D]*kxt[IDX3D]*ky[IDX3D]*fldi.vx[IDX3D] / k2t[IDX3D];	// m=1, n=1
					w2[IDX3D] = I*kxt[IDX3D]*ky[IDX3D]*ky[IDX3D]*fldi.vx[IDX3D] / k2t[IDX3D];	// m=1, n=2
					w3[IDX3D] = I*kxt[IDX3D]*kz[IDX3D]*ky[IDX3D]*fldi.vx[IDX3D] / k2t[IDX3D];	// m=1, n=3
					//w4[IDX3D] = I*ky[IDX3D]*kxt[IDX3D]*ky[IDX3D]*fldi.vx[IDX3D] / k2t[IDX3D];	// m=2, n=1
					w5[IDX3D] = I*ky[IDX3D]*ky[IDX3D]*ky[IDX3D]*fldi.vx[IDX3D] / k2t[IDX3D];	// m=2, n=2
					w6[IDX3D] = I*ky[IDX3D]*kz[IDX3D]*ky[IDX3D]*fldi.vx[IDX3D] / k2t[IDX3D];	// m=2, n=3
					//w7[IDX3D] = I*kz[IDX3D]*kxt[IDX3D]*ky[IDX3D]*fldi.vx[IDX3D] / k2t[IDX3D];	// m=3, n=1
					//w8[IDX3D] = I*kz[IDX3D]*ky[IDX3D]*ky[IDX3D]*fldi.vx[IDX3D] / k2t[IDX3D];	// m=3, n=2
					w9[IDX3D] = I*kz[IDX3D]*kz[IDX3D]*ky[IDX3D]*fldi.vx[IDX3D] / k2t[IDX3D];	// m=3, n=3
				} else {
					w1[IDX3D] = 0.0;
					w2[IDX3D] = 0.0;
					w3[IDX3D] = 0.0;
					//w4[IDX3D] = 0.0;
					w5[IDX3D] = 0.0;
					w6[IDX3D] = 0.0;
					//w7[IDX3D] = 0.0;
					//w8[IDX3D] = 0.0;
					w9[IDX3D] = 0.0;
				}
					w10[IDX3D] = fldi.vx[IDX3D];
					w11[IDX3D] = fldi.vy[IDX3D];
					w12[IDX3D] = fldi.vz[IDX3D];
			}
		}
	}
	gfft_c2r(w1);
	gfft_c2r(w2);
	gfft_c2r(w3);
	//gfft_c2r(w4);
	gfft_c2r(w5);
	gfft_c2r(w6);
	//gfft_c2r(w7);
	//gfft_c2r(w8);
	gfft_c2r(w9);
	gfft_c2r(w10);
	gfft_c2r(w11);
	gfft_c2r(w12);

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr1[i] = (wr1[i]*wr10[i] + wr2[i]*wr11[i] + wr3[i]*wr12[i]) / ((double) NTOTAL*NTOTAL);	// n=1
		wr2[i] = (wr2[i]*wr10[i] + wr5[i]*wr11[i] + wr6[i]*wr12[i]) / ((double) NTOTAL*NTOTAL);	// n=2
		wr3[i] = (wr3[i]*wr10[i] + wr6[i]*wr11[i] + wr9[i]*wr12[i]) / ((double) NTOTAL*NTOTAL);	// n=3
	}
	gfft_r2c(wr1);
	gfft_r2c(wr2);
	gfft_r2c(wr3);
	// Now acutally compute 4*< v_m (d_m d_n d^-2 d_2 v_1) (d_n d^-2 d_2 v_1)>
	//                     +4*< v_m (d_m d_n d^-2 d_2 v_1) (d_n d^-2 d_1 v_2)>
	for(i=0; i < NX_COMPLEX/NPROC; i++){
		for(j=0; j < NY_COMPLEX; j++){
			for(k=0; k < NZ_COMPLEX; k++){
				if(!(i==0 && j==0 && k==0 && rank==0)){
#ifdef WITH_2D
					if(j==0)
#else
					if(k==0){
#endif
						// k=0, we have all of the modes
						E_output = E_output + 4.0*creal((conj(w1[IDX3D])*kxt[IDX3D] + conj(w2[IDX3D])*ky[IDX3D]
						                                 + conj(w3[IDX3D])*kz[IDX3D])*(ky[IDX3D]*fldi.vx[IDX3D]/k2t[IDX3D]) / ((double) NTOTAL*NTOTAL));
						E_output = E_output + 4.0*creal((conj(w1[IDX3D])*kxt[IDX3D] + conj(w2[IDX3D])*ky[IDX3D]
						                                 + conj(w3[IDX3D])*kz[IDX3D])*(kxt[IDX3D]*fldi.vy[IDX3D]/k2t[IDX3D]) / ((double) NTOTAL*NTOTAL));
					} else {
						// k>0, one half of the complex plane is represented
						E_output = E_output + 2.0*4.0*creal((conj(w1[IDX3D])*kxt[IDX3D] + conj(w2[IDX3D])*ky[IDX3D]
						                                     + conj(w3[IDX3D])*kz[IDX3D])*(ky[IDX3D]*fldi.vx[IDX3D]/k2t[IDX3D]) / ((double) NTOTAL*NTOTAL));
						E_output = E_output + 2.0*4.0*creal((conj(w1[IDX3D])*kxt[IDX3D] + conj(w2[IDX3D])*ky[IDX3D]
						                                     + conj(w3[IDX3D])*kz[IDX3D])*(kxt[IDX3D]*fldi.vy[IDX3D]/k2t[IDX3D]) / ((double) NTOTAL*NTOTAL));
					}
				}
			}
		}
	}
	/* End computation of 4*E_2121 + 4*E_2112 */

	/* Restore the w* and wr* memory */
	// Transform a few basic variables.
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = fldi.vx[i];
		w2[i] = fldi.vy[i];
		w3[i] = fldi.vz[i];
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
#endif
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
	/* Finished restoring the w* and wr* memory */

	return B_output + C_output + D_output + E_output;

}
#endif



/*********************************************/
/*  Compute components of 
    b_ij = (<u_i u_j>/<u_k u_k>) - (delat_ij / 3)
	Output b_ij b_ji
*/
/*********************************************/

double RottaEtaSquared(const struct Field fldi){
	
	int i,j,k;	// For addressing wavenumbers via the IDX3D macro
	double ev = 0.0;
	double b11 = 0.0, b12 = 0.0, b13 = 0.0, b22 = 0.0, b23 = 0.0, b33 = 0.0;

	// Compute the energy first
	for(i=0; i< NX_COMPLEX/NPROC; i++){
		for(j=0; j< NY_COMPLEX; j++){
			for(k=0; k< NZ_COMPLEX; k++){
#ifdef WITH_2D
				if( j == 0){
#else
				if( k == 0){
#endif
					// k=0, we have all the modes.
					ev = ev + creal(fldi.vx[IDX3D]*conj(fldi.vx[IDX3D]) + fldi.vy[IDX3D]*conj(fldi.vy[IDX3D])
					                + fldi.vz[IDX3D]*conj(fldi.vz[IDX3D])) / ((double) NTOTAL*NTOTAL);
				} else {
					// k>0, only half of the complex plane is represented.
					ev = ev + 2.0*creal(fldi.vx[IDX3D]*conj(fldi.vx[IDX3D]) + fldi.vy[IDX3D]*conj(fldi.vy[IDX3D])
					                    + fldi.vz[IDX3D]*conj(fldi.vz[IDX3D])) / ((double) NTOTAL*NTOTAL);
				}
			}
		}
	}

	// Get the actual energy, not just that of the wavemodes on this process
	reduce(&ev, 1);

	// Now compute the b_ij components
	for(i=0; i< NX_COMPLEX/NPROC; i++){
		for(j=0; j< NY_COMPLEX; j++){
			for(k=0; k< NZ_COMPLEX; k++){
#ifdef WITH_2D
				if( j == 0){
#else
				if( k == 0){
#endif
					// k=0, we have all the modes.
					b11 = b11 + creal(fldi.vx[IDX3D]*conj(fldi.vx[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
					b22 = b22 + creal(fldi.vy[IDX3D]*conj(fldi.vy[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
					b33 = b33 + creal(fldi.vz[IDX3D]*conj(fldi.vz[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
					b12 = b12 + creal(fldi.vx[IDX3D]*conj(fldi.vy[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
					b13 = b13 + creal(fldi.vx[IDX3D]*conj(fldi.vz[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
					b23 = b23 + creal(fldi.vy[IDX3D]*conj(fldi.vz[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
				} else {
					// k>0, only half of the complex plane is represented.
					b11 = b11 + 2.0*creal(fldi.vx[IDX3D]*conj(fldi.vx[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
					b22 = b22 + 2.0*creal(fldi.vy[IDX3D]*conj(fldi.vy[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
					b33 = b33 + 2.0*creal(fldi.vz[IDX3D]*conj(fldi.vz[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
					b12 = b12 + 2.0*creal(fldi.vx[IDX3D]*conj(fldi.vy[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
					b13 = b13 + 2.0*creal(fldi.vx[IDX3D]*conj(fldi.vz[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
					b23 = b23 + 2.0*creal(fldi.vy[IDX3D]*conj(fldi.vz[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
				}
			}
		}
	}

	// reduce() step has to come HERE, while still linear!!!!
	reduce(&b11, 1);
	reduce(&b22, 1);
	reduce(&b33, 1);
	reduce(&b12, 1);
	reduce(&b13, 1);
	reduce(&b23, 1);

	// Make the tensor traceless
	b11 = b11 - 1.0/3.0;
	b22 = b22 - 1.0/3.0;
	b33 = b33 - 1.0/3.0;

	return (b11*b11 + b22*b22 + b33*b33 + 2.0*(b12*b12 + b13*b13 + b23*b23))/6.0;

}


/*********************************************/
/*  Compute components of 
    b_ij = (<u_i u_j>/<u_k u_k>) - (delat_ij / 3)
	Output b_ij b_jik b_ki
*/
/*********************************************/
double RottaXiSquared(const struct Field fldi){

	int i,j,k;	// For addressing wavenumbers via the IDX3D macro
	double ev = 0.0;
	double b11 = 0.0, b12 = 0.0, b13 = 0.0, b22 = 0.0, b23 = 0.0, b33 = 0.0;

	// Compute the energy first
	for(i=0; i< NX_COMPLEX/NPROC; i++){
		for(j=0; j< NY_COMPLEX; j++){
			for(k=0; k< NZ_COMPLEX; k++){
#ifdef WITH_2D
				if( j == 0){
#else
				if( k == 0){
#endif
					// k=0, we have all the modes.
					ev = ev + creal(fldi.vx[IDX3D]*conj(fldi.vx[IDX3D]) + fldi.vy[IDX3D]*conj(fldi.vy[IDX3D])
					                + fldi.vz[IDX3D]*conj(fldi.vz[IDX3D])) / ((double) NTOTAL*NTOTAL);
				} else {
					// k>0, only half of the complex plane is represented.
					ev = ev + 2.0*creal(fldi.vx[IDX3D]*conj(fldi.vx[IDX3D]) + fldi.vy[IDX3D]*conj(fldi.vy[IDX3D])
					                    + fldi.vz[IDX3D]*conj(fldi.vz[IDX3D])) / ((double) NTOTAL*NTOTAL);
				}
			}
		}
	}

	// Get the actual energy, not just that of the wavemodes on this process
	reduce(&ev, 1);

	// Now compute the b_ij components
	for(i=0; i< NX_COMPLEX/NPROC; i++){
		for(j=0; j< NY_COMPLEX; j++){
			for(k=0; k< NZ_COMPLEX; k++){
#ifdef WITH_2D
				if( j == 0){
#else
				if( k == 0){
#endif
					// k=0, we have all the modes.
					b11 = b11 + creal(fldi.vx[IDX3D]*conj(fldi.vx[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
					b22 = b22 + creal(fldi.vy[IDX3D]*conj(fldi.vy[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
					b33 = b33 + creal(fldi.vz[IDX3D]*conj(fldi.vz[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
					b12 = b12 + creal(fldi.vx[IDX3D]*conj(fldi.vy[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
					b13 = b13 + creal(fldi.vx[IDX3D]*conj(fldi.vz[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
					b23 = b23 + creal(fldi.vy[IDX3D]*conj(fldi.vz[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
				} else {
					// k>0, only half of the complex plane is represented.
					b11 = b11 + 2.0*creal(fldi.vx[IDX3D]*conj(fldi.vx[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
					b22 = b22 + 2.0*creal(fldi.vy[IDX3D]*conj(fldi.vy[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
					b33 = b33 + 2.0*creal(fldi.vz[IDX3D]*conj(fldi.vz[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
					b12 = b12 + 2.0*creal(fldi.vx[IDX3D]*conj(fldi.vy[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
					b13 = b13 + 2.0*creal(fldi.vx[IDX3D]*conj(fldi.vz[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
					b23 = b23 + 2.0*creal(fldi.vy[IDX3D]*conj(fldi.vz[IDX3D])/ev) / ((double) NTOTAL*NTOTAL);
				}
			}
		}
	}

	// reduce() step has to come HERE, while still linear!!!!
	reduce(&b11, 1);
	reduce(&b22, 1);
	reduce(&b33, 1);
	reduce(&b12, 1);
	reduce(&b13, 1);
	reduce(&b23, 1);

	// Make the tensor traceless
	b11 = b11 - 1.0/3.0;
	b22 = b22 - 1.0/3.0;
	b33 = b33 - 1.0/3.0;

	return (b11*b11*b11 + b22*b22*b22 + b33*b33*b33 + 6.0*b12*b23*b13
	        + 3.0*b11*(b12*b12 + b13*b13) + 3.0*b22*(b12*b12 + b23*b23) + 3.0*b33*(b13*b13 + b23*b23))/6.0;

}
/*********************************************/
/*  Added by Harry Braviner as part of the 
 *  OSCILLATORY_SHEAR additions
	Note the factor of two difference from 
    the energy function above correlator computes 
    < q1 q2 > whereas energy computes 0.5*<q q>
*/
/*********************************************/

double correlator(const double complex q1[], const double complex q2[]){
	
    int i,j,k;
	double c_tot = 0.0;

	for(i=0; i< NX_COMPLEX/NPROC; i++){
		for(j=0; j< NY_COMPLEX; j++){
			for(k=0; k< NZ_COMPLEX; k++){
#ifdef WITH_2D
				if( j == 0)
#else
				if( k == 0)
#endif
					// k=0, we have all the modes.
					c_tot = c_tot + creal(q1[ IDX3D ] * conj( q2[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
				else
					// k>0, only half of the complex plane is represented.
					c_tot = c_tot + creal(2.0 * q1[ IDX3D ] * conj( q2[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
			}
		}
	}
	return c_tot;
}

#ifdef SHELL_R_XY
// Essentially the same as the write_spectrum function of output_spectrum.c, but this saves the spectrum to an array rather than writing it to a file
void make_spectrum(const double complex wi[], const double complex wj[], double shell_spectrum[]){
	DEBUG_START_FUNC;
	int i,j,k,m;
	int nbin;
	
	nbin = (int) ceil(kmax / OUTPUT_SPECTRUM_K_BIN);
	
	for( i = 0; i < MAX_N_BIN; i++ )
		shell_spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor( pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if ( m < nbin) {
#ifdef WITH_2D
					if( j == 0)
#else
					if( k == 0) 
#endif
						// k=0, we have all the modes.
						shell_spectrum[ m ] = shell_spectrum[ m ] + creal( wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
					else
						// k>0, only half of the complex plane is represented.
						shell_spectrum[ m ] = shell_spectrum[ m ] + creal( 2.0 * wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
				}
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin; m++)
		reduce(&shell_spectrum[m], 1);
#endif

	DEBUG_END_FUNC;
	return;
}
#endif

/********************************************/
/*
    Return the localtime in seconds. Use different
    implementation depending on the avaiable
    libraries  
*/
/********************************************/
double get_c_time(void) {
#ifdef MPI_SUPPORT
	// We have MPI
	return(MPI_Wtime());
#else
#ifdef _OPENMP
	// We don't have MPI, but we have OpenMP
	return(omp_get_wtime());
#else
	// We really have nothing...
	clock_t now;
	now = clock();
	
	return( (double) now / ( (double) CLOCKS_PER_SEC));

#endif
#endif
}

/******************************************/
/**
	Reduce a variable over all the avaiable processes
	Can add a value on all the process, find a maximum
	or a minimum.
	
	NOTE: This routine makes sense only when MPI_SUPPORT
	is set. If not, this routine does nothing.
	
	@param *var: variable to be reduced
	@param op: operation needed to be done. Can be set to:
		1= Sum over all the processes
		2= Find the maximum over all the processes
		3= Find the minimum over all the processes
*/
/*******************************************/
	
void reduce(double *var, const int op) {
	// op=1 ADD
	// op=2 Max
	// op=3 Min
	
#ifdef MPI_SUPPORT
	double mpi_temp;
	
	mpi_temp=*var;
	
	if(op==1) MPI_Allreduce( &mpi_temp, var, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	if(op==2) MPI_Allreduce( &mpi_temp, var, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	if(op==3) MPI_Allreduce( &mpi_temp, var, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	
#endif	

	// If no MPI, then this routine does nothing...
	return;
}

/* ****************************************************************************/
/** Determines if the machine is little-endian.  If so, 
    it will force the data to be big-endian. 
	@param in_number floating point number to be converted in big endian */
/* *************************************************************************** */

float big_endian(float in_number)
{
    static int doneTest = 0;
    static int shouldSwap = 0;
	
    if (!doneTest)
    {
        int tmp1 = 1;
        unsigned char *tmp2 = (unsigned char *) &tmp1;
        if (*tmp2 != 0)
            shouldSwap = 1;
        doneTest = 1;
    }

    if (shouldSwap)
    {
		unsigned char *bytes = (unsigned char*) &in_number;
        unsigned char tmp = bytes[0];
        bytes[0] = bytes[3];
        bytes[3] = tmp;
        tmp = bytes[1];
        bytes[1] = bytes[2];
        bytes[2] = tmp;
    }
	return(in_number);
}

/******************************************************************************/
/** Test a double for Not a Number error 
	@param xi double to be checked */
/* ************************************************************************** */

void c_nan(double xi, const char ErrorRoutine[], const int line, const char Filename[]) {
	if(isnan(xi)) {
		error_h( 3, "Not a number detected", ErrorRoutine, line, Filename);
	}
	return;
}

#ifdef COMPRESSIBLE
/***************************************************************/
/** Check that the field is definite positive *******************/
/****************************************************************/

void check_positivity(double *wri) {
	int i;
#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for(i=0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		if(wri[i] < 1.0e-2*NTOTAL) wri[i] = 1.0e-2*NTOTAL;
	}
}
#endif
