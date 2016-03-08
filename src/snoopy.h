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



#ifndef __SNOOPY_H__
#define __SNOOPY_H__

#include <math.h>
#include <complex.h>
#include <fftw3.h>

#ifdef MPI_SUPPORT
#include <mpi.h>
#ifdef FFTW3_MPI_SUPPORT
#include <fftw3-mpi.h>
#endif
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "gvars.h"

#ifdef WITH_2D
#undef NZ
#define NZ		1		// enforce NZ when in 2D
#endif

#include "error.h"

// Enforce No elssaser when not using MHD
#ifdef ELSASSER_FORMULATION
#ifndef MHD
#undef ELSASSER_FORMULATION
#endif
#endif

#ifdef MPI_SUPPORT
#define		MPI_Printf			if (rank==0) printf
#else
#define		MPI_Printf			printf
#endif


// fix NPROC if MPI_SUPPORT is disabled
#ifndef	MPI_SUPPORT
#ifndef		NPROC
#define		NPROC			1
#endif
#endif

#define		NTOTAL			NX * NY * NZ	/**< Total number of grid points over all the processors */

#ifdef WITH_2D
#define		NX_COMPLEX		NX				/**< Number of complex point in X over all the processes (Fourier space) */
#define		NY_COMPLEX		(NY / 2 + 1)	/**< Number of complex point in Y over all the processes (Fourier space) */
#define		NZ_COMPLEX		1				/**< Number of complex point in Z over all the processes (Fourier space) */
#else
#define		NX_COMPLEX		NX				/**< Number of complex point in X over all the processes (Fourier space) */
#define		NY_COMPLEX		NY				/**< Number of complex point in Y over all the processes (Fourier space) */
#define		NZ_COMPLEX		(NZ / 2 + 1)	/**< Number of complex point in Z over all the processes (Fourier space) */
#endif

#define		NTOTAL_COMPLEX	(NX_COMPLEX * NY_COMPLEX * NZ_COMPLEX / NPROC)	/**< Number of complex points in one MPI process */

#define		IDX3D			(k + j * NZ_COMPLEX + NZ_COMPLEX * NY_COMPLEX * i)  /**< General wrapper for 3D Arrays */

// Filenames

#define	OUTPUT_DUMP					"dump.dmp"			/**< Dump files filename. */
#define OUTPUT_DUMP_SAV				"dump_sav.dmp"      /**< Previous (saved) output dump. */
#define OUTPUT_DUMP_WRITE			"dump_write.dmp"	/**< dump currently written. */


// Structures

#ifdef WITH_PARTICLES
struct Particle	{
	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	double mass;	// unused for the moment
	double stime;	// Stopping time
};
#endif

struct Field {
	int	   nfield;
	double complex **farray;
	char   **fname;
	
	 // These are the actual fields. 
	double complex *vx;
	double complex *vy;
	double complex *vz;
#ifdef WITH_OL2012B
	double complex *vdotx;
	double complex *vdoty;
	double complex *vdotz;
#endif
#ifdef KINEMATIC_GROWTH
    double complex *bx0;
    double complex *by0;
    double complex *bz0;
#endif   
#ifdef BOUSSINESQ
	double complex *th;
#endif
#ifdef MHD
	double complex *bx;
	double complex *by;
	double complex *bz;
#endif
#ifdef COMPRESSIBLE
	double complex *d;
#endif
#ifdef WITH_PARTICLES
	struct Particle *part;
#endif
};

struct VarName {
	int		length;		// Number of varnames
	char**	name;		// Varnames (need to be allocated properly)
};

struct Parameters {
	// Physics Parameters
	double lx;				/**< Box length in X*/
	double ly;				/**< Box length in Y*/
	double lz;				/**< Box length in Z*/
	
	double reynolds;		/**< Reynolds number (actully the inverse of the viscosity) */
	double reynolds_m;		/**< Magnetic Reynolds number (actully the inverse of the resistivity)  Used only when MHD is on*/
	double etta;            /**< Can define resistivity itself. set in common.c*/ 

	double reynolds_th;		/**< Thermal Reynolds number (actully the inverse of the thermal diffusivity)  Used only when Boussinesq is on*/
	double reynolds_B;		/**< Reynolds number based on Braginskii viscosity */
	
	double x_hall;			/**< Hall parameter */
	
	double N2;				/**< Brunt Vaissala frequency squared */
	
	double omega;			/**< Vertical rotation rate (if Shear=1, Keplerian if found for (2.0/3.0). Only when WITH_ROTATION is on. */
	
	double shear;			/**< Shear rate (only when WITH_SHEAR is on) */
	
	double omega_shear;		/**< Pulsation of the time dependant shear (only when WITH_SHEAR and TIME_DEPENDANT_SHEAR is on, or alternatively WITH_LINEAR_TIDE) */

	// The following should only be used if OSCILLATORY_SHEAR is defined
	double oscillatory_shear_freq;	/**< Frequency of the oscillatory shear */
	double oscillatory_shear_amp;	/**< Amplitude of the oscillatory shear */
	
	double cs;				/**< Sound speed (only used when compressible is on) */
	
	// Particles parameters
	int    particles_n;		/**< Number of particles */
	double particles_mass;  /**< Mass of the particles */
	double particles_stime;	/**< Stopping time of the particles */
	double particles_dg_ratio;	/**< Dust to gas mass ratio for the particles feedback */
	double particles_epsilon;	/**< Pressure gradient epsilon */
	
	// Code parameters
	
	double cfl;				/**< CFL safety factor. Should be smaller than sqrt(3) for RK3 to be stable.*/
    double cfl_diss;        /**< CFL safety factor for the dissipative effects. Should be smaller than cfl factor above. */
	double safety_source;	/**< Safety factor for SHEAR, Coriolis and Boussinesq terms (should be ~0.2 for accuracy) */
	double steps_per_shear;	// Number of timesteps to take per shear period if oscillatory shear is enabled
	
	double t_initial;		/**< Initial time of the simulation */
	double t_final;			/**< Simulation will stop if it reaches this time */
	double max_t_elapsed;	/**< Maximum elapsed time (in hours). Will stop after this elapsed time */
	
	int    interface_check;	/**< Number of loops between two checks for a user input. On slow filesystems, increase this number */
	int    interface_output_file;	/**< Set this option to force code outputs to a file instead of the screen */
	
	int    force_symmetries;	/**< set to enforce spectral symmetries and mean flow to zero. Useful when N^2 or kappa^2 < 0. (see enforce_symm() )*/
	int    symmetries_step;		/**< Number of loops between which the symmetries are enforced. Should be around ~20 for Boussinesq convection*/
	
	int    antialiasing;		/**< 2/3 Antialisaing rule. Could be removed if you assume is unimportant (untested feature). */
	
	int    restart;
	
	// Output parameters
	double toutput_time;		/**< Time between two outputs in the timevar file */
	double toutput_flow;		/**< Time between two snapshot outputs */
	double toutput_dump;		/**< Time between two restart dump outputs (restart dump are erased) */
	
	int		output_vorticity;	    /**< Output the vorticity field in the 3D snapshots */
    int		output_magnetic_field;	/**< Output the <B> field in the 3D snapshots */

/* Start of timeseries code */
	int output_timeseries_1;	// 1 if I'm to output a timeseries of u(x1), 0 otherwise
	int output_timeseries_2;	// 1 if I'm to output a timeseries of u(x2), 0 otherwise
	double timeseries_1x;
	double timeseries_1y;
	double timeseries_1z;
	double timeseries_2x;
	double timeseries_2y;
	double timeseries_2z;
	double timeseries_time;		// Time between two velocity timeseries outputs
	int timeseries_1nx;	// Indices of the collocation points that we want to
	int timeseries_1ny;	// take a timeseries at.
	int timeseries_1nz;
	int timeseries_2nx;
	int timeseries_2ny;
	int timeseries_2nz;
/* End of timeseries code */
	
	struct VarName timevar_vars;	/**< Name of the variables needed in the timevar file */
	
	// initial conditions
	int	   init_vortex;			/**< Add a 2D Kida vortex in the box. Assumes S=1. Requires b>a*/
	double vortex_a;			/**< x dimension of the vortex */
	double vortex_b;			/**< y dimension of the vortex */
	
	int    init_spatial_structure;	/**< Init a user-defined spatial structure (see initflow.c) */
	
	int	   init_large_scale_noise;	/**< Init a large scale random noise */
	double per_amplitude_large;		/**< Amplitude of the large scale random noise */
	double noise_cut_length;		/**< Wavelength over which the noise is applied */
	
	int	   init_large_scale_2D_noise;	/**< Init a large scale 2D (x,y) random noise  */
	double per_amplitude_large_2D;		/**< Amplitude of the 2D large scale random noise */
	double noise_cut_length_2D;		    /**< Wavelength over which the 2D noise is applied */
	
	int    init_white_noise;		/**< Init a random white noise on all the modes */
	double per_amplitude_noise;		/**< total amplitude of the perturbation */
	
	int    init_ABC_flow;		/**< Init an ABC flow with wavenumber 1 */
    double ABC_flow_A;		/**< Initial amplitude of the z-dependent flow */
	double ABC_flow_B;		/**< Initial amplitude of the x-dependent flow */
	double ABC_flow_C;		/**< Initial amplitude of the y-dependent flow */
	long    ABC_flow_kx;
	long    ABC_flow_ky;
	long    ABC_flow_kz;

    
    int    init_modified_ABC_flow;	/**< Init a modified ABC flow with wavenumber 1 */
    double modified_ABC_flow_A;		
    double modified_ABC_flow_B;		
    double modified_ABC_flow_C;		
    double modified_ABC_flow_D;		
    long   modified_ABC_flow_kx;
	long   modified_ABC_flow_ky;
	long   modified_ABC_flow_kz;
    long   modified_ABC_flow_m;

	int    init_Magnetic_field;		/**< Init a magnetic field */

	int	init_dominant_scale;		/**< Assign kx,ky,kz = 1 more energy on the onset */

    int    init_No_Noise;           /**< Set every wavenumber B to zero */

	int    init_mean_field;			/**< Force the mean magnetic field to a given value. */
	double bx0;						/**< Mean magnetic field in the x direction */
	double by0;						/**< Mean magnetic field in the y direction */
	double bz0;						/**< Mean magnetic field in the z direction */
	
	int    init_dump;				/**< Use a dump file as an initial condition (everything else, including t, noutput (...) is reinitialised) */
	
	int	   init_bench;				/**< Init the Benchmark initial conditions */
};

#endif

