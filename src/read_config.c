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

#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "debug.h"
#include "libconfig/libconfig.h"

#define CONFIG_FILENAME		"snoopy.cfg"

void read_config() {
	// Read the config file and initialize everyting
	config_t	config;		// Initialize the structure
	config_setting_t * setting;	// a setting structure
	long tmp_v;
	int i,n;
	
	const char * configname;
	
	const char * temp_string;
	
	
	DEBUG_START_FUNC;
	
	if(rank==0) {
		config_init(&config);
	
		if(!config_read_file(&config, CONFIG_FILENAME)) {
			MPI_Printf("Error reading configuration file in line %d: %s\n", config_error_line(&config), config_error_text(&config));
			ERROR_HANDLER(ERROR_CRITICAL, "Failed to read the configuration file");
		}

		if(config_lookup_string(&config,"configname",&configname)) {
			MPI_Printf("Using config file: %s.\n",configname);
		}
		// read physics parameters-------------------------------------------------------------------------------
		if(!config_lookup_float(&config, "physics.boxsize.[0]",&param.lx)) {
			param.lx = 1.0;
		}
		if(!config_lookup_float(&config, "physics.boxsize.[1]",&param.ly)) {
		param.ly = 1.0;
		}
		if(!config_lookup_float(&config, "physics.boxsize.[2]",&param.lz)) {
			param.lz = 1.0;
		}
		
		if(!config_lookup_float(&config, "physics.reynolds",&param.reynolds)) {
			param.reynolds = 1.0;
		}	
	
		if(!config_lookup_float(&config, "physics.reynolds_magnetic",&param.reynolds_m)) {
			param.reynolds_m = 1.0;
		}
        if(!config_lookup_float(&config, "physics.etta",&param.etta)) {
			param.etta = 1.0;
		}
		if(!config_lookup_float(&config, "physics.reynolds_thermic",&param.reynolds_th)) {
			param.reynolds_th = 1.0;
		}
	
		if(!config_lookup_float(&config, "physics.reynolds_Braginskii",&param.reynolds_B)) {
			param.reynolds_B = 1.0;
		}
		
		if(!config_lookup_float(&config, "physics.x_Hall",&param.x_hall)) {
			param.x_hall = 1.0;
		}
		
		if(!config_lookup_float(&config, "physics.brunt_vaissala_squared",&param.N2)) {
			param.N2 = 0.0;
		}
	
		if(!config_lookup_float(&config, "physics.omega",&param.omega)) {
			param.omega = 0.0;
		}
#ifndef WITH_ROTATION
		// Omega should be forced to zero in order to be fool-proof
		param.omega = 0.0;
#endif
		if(!config_lookup_float(&config, "physics.shear",&param.shear)) {
			param.shear = 0.0;
		}
#ifndef WITH_SHEAR
		// same for the shear
		param.shear = 0.0;
#endif	
		if(!config_lookup_float(&config, "physics.omega_shear",&param.omega_shear)) {
			param.omega_shear = 0.0;
		}
		
		if(!config_lookup_float(&config, "physics.sound_speed",&param.cs)) {
			param.cs = 1.0;
		}

		// The following should only be used if OSCILLATORY_SHEAR is defined
		if(!config_lookup_float(&config, "physics.oscillatory_shear_freq",&param.oscillatory_shear_freq)){
			param.oscillatory_shear_freq = 0.0;
		}
		if(!config_lookup_float(&config, "physics.oscillatory_shear_amp",&param.oscillatory_shear_amp)){
			param.oscillatory_shear_amp = 0.0;
		}
	
		// Particles parameters-------------------------------------------------------------------------------------
		if(!config_lookup_int(&config, "particles.n",&tmp_v)) {
			param.particles_n = 1000;
		}
		else {
			param.particles_n = (int) tmp_v;
		}
		
		if(!config_lookup_float(&config, "particles.mass",&param.particles_mass)) {
			param.particles_mass = 1.0;
		}
		
		if(!config_lookup_float(&config, "particles.stime",&param.particles_stime)) {
			param.particles_stime = 1.0;
		}
		
		if(!config_lookup_float(&config, "particles.dg_ratio",&param.particles_dg_ratio)) {
			param.particles_dg_ratio = 0.01;
		}
		
		if(!config_lookup_float(&config, "particles.epsilon",&param.particles_epsilon)) {
			param.particles_epsilon = 0.1;
		}
		
		// Code parameters-------------------------------------------------------------------------------------
	
		if(!config_lookup_float(&config, "code.cfl",&param.cfl)) {
			param.cfl = 1.5;
		}
        if(!config_lookup_float(&config, "code.cfl_diss",&param.cfl_diss)) {
			param.cfl_diss = 1.5;
		}
        if(!config_lookup_float(&config, "code.safety_source",&param.safety_source)) {
			param.safety_source = 0.2;
		}
		if(!config_lookup_float(&config, "code.steps_per_shear",&param.steps_per_shear)) {
			param.steps_per_shear = 10.0;
		}
		if(!config_lookup_float(&config, "code.t_initial",&param.t_initial)) {
			param.t_initial = 0.0;
		}
		if(!config_lookup_float(&config, "code.t_final",&param.t_final)) {
			param.t_final = 1.0;
		}
		if(!config_lookup_float(&config, "code.max_t_elapsed",&param.max_t_elapsed)) {
			param.max_t_elapsed = 1e30;
		}
		if(!config_lookup_int(&config, "code.interface_check",&tmp_v)) {
			param.interface_check = 5;
		}
		else {
			param.interface_check = (int) tmp_v;
		}
		if(!config_lookup_bool(&config, "code.interface_output_file",&param.interface_output_file)) {
			param.interface_output_file = 0;
		}
		if(!config_lookup_bool(&config, "code.force_symmetries",&param.force_symmetries)) {
			param.force_symmetries = 0;
		}
		if(!config_lookup_int(&config, "code.symmetries_step",&tmp_v)) {
			param.symmetries_step = 20;
		}
		else {
			param.symmetries_step = (int) tmp_v;
		}
		if(!config_lookup_bool(&config, "code.antialiasing",&param.antialiasing)) {
			param.antialiasing = 1;
		}
		if(!config_lookup_bool(&config, "code.restart",&param.restart)) {
			param.restart = 0;
		}

		// Output parameters-------------------------------------------------------------------------------------
		if(!config_lookup_float(&config, "output.timevar_step",&param.toutput_time)) {
			param.toutput_time = 1.0;
		}
		if(!config_lookup_float(&config, "output.snapshot_step",&param.toutput_flow)) {
			param.toutput_flow = 1.0;
		}
		if(!config_lookup_float(&config, "output.dump_step",&param.toutput_dump)) {
			param.toutput_dump = 1.0;
		}
		if(!config_lookup_bool(&config, "output.vorticity",&param.output_vorticity)) {
			param.output_vorticity = 0;
		}
        if(!config_lookup_bool(&config, "output.magnetic_field",&param.output_magnetic_field)) {
			param.output_magnetic_field = 0;
		}

		// find which parameters are requested in the timevar file
		setting = config_lookup(&config, "output.timevar_vars");
		
		if(setting == NULL) {
			ERROR_HANDLER(ERROR_WARNING, "You did not provide any variable in timevar outputs");
		}
		else {
			param.timevar_vars.length = config_setting_length( setting );
		
			// Allocate output_vars
			param.timevar_vars.name = malloc( param.timevar_vars.length * sizeof(char*) );
		
			for(i = 0 ; i < param.timevar_vars.length ; i++) {
				temp_string = config_setting_get_string_elem( setting, i);
			
				// Allocate the string
				param.timevar_vars.name[i] = malloc( sizeof(char) * (strlen(temp_string) + 1));
			
				// Copy the string in the right location
				strcpy(param.timevar_vars.name[i], temp_string);
			}
		}
/* Start of timeseries code */
		param.output_timeseries_1 = 1;
		param.output_timeseries_2 = 1;
		if(!(  config_lookup_float(&config, "output.timeseries1x", &param.timeseries_1x)
		     &&config_lookup_float(&config, "output.timeseries1y", &param.timeseries_1y)
		     &&config_lookup_float(&config, "output.timeseries1z", &param.timeseries_1z))){
			param.output_timeseries_1 = 0;
		}
		if(!(  config_lookup_float(&config, "output.timeseries2x", &param.timeseries_2x)
		     &&config_lookup_float(&config, "output.timeseries2y", &param.timeseries_2y)
		     &&config_lookup_float(&config, "output.timeseries2z", &param.timeseries_2z))){
			param.output_timeseries_2 = 0;
		}
		if(!config_lookup_float(&config, "output.timeseries_step", &param.timeseries_time)){
			param.output_timeseries_1 = 0;
			param.output_timeseries_2 = 0;
		}
/* End of timeseries code */
		
		
		// Initial conditions parameters-------------------------------------------------------------------------
		if(!config_lookup_bool(&config, "init.vortex.enable",&param.init_vortex)) {
			param.init_vortex = 0;
		}
		if(!config_lookup_float(&config, "init.vortex.a",&param.vortex_a)) {
			param.vortex_a = 1.0;
		}
		if(!config_lookup_float(&config, "init.vortex.b",&param.vortex_b)) {
			param.vortex_b = 2.0;
		}
		if(!config_lookup_bool(&config, "init.spatial_structure",&param.init_spatial_structure)) {
			param.init_spatial_structure = 0;
		}
		if(!config_lookup_bool(&config, "init.large_scale_noise.enable",&param.init_large_scale_noise)) {
			param.init_large_scale_noise = 0;
		}	
		if(!config_lookup_float(&config, "init.large_scale_noise.amplitude",&param.per_amplitude_large)) {
			param.per_amplitude_large = 0.0;
		}	
		if(!config_lookup_float(&config, "init.large_scale_noise.cut_length",&param.noise_cut_length)) {
			param.noise_cut_length = 0.0;
		}
		if(!config_lookup_bool(&config, "init.large_scale_2D_noise.enable",&param.init_large_scale_2D_noise)) {
			param.init_large_scale_2D_noise = 0;
		}	
		if(!config_lookup_float(&config, "init.large_scale_2D_noise.amplitude",&param.per_amplitude_large_2D)) {
			param.per_amplitude_large_2D = 0.0;
		}
		if(!config_lookup_float(&config, "init.large_scale_2D_noise.cut_length",&param.noise_cut_length_2D)) {
			param.noise_cut_length_2D = 0.0;
		}
		if(!config_lookup_bool(&config, "init.white_noise.enable",&param.init_white_noise)) {
			param.init_white_noise = 0;
		}
		if(!config_lookup_float(&config, "init.white_noise.amplitude",&param.per_amplitude_noise)) {
			param.per_amplitude_noise = 0.0;
		}
		if(!config_lookup_bool(&config, "init.ABC_flow.enable",&param.init_ABC_flow)) {
			param.init_ABC_flow = 0;
		}
        if(!config_lookup_bool(&config, "init.modified_ABC_flow.enable",&param.init_modified_ABC_flow)) {
			param.init_modified_ABC_flow = 0;
		}
		if(!config_lookup_bool(&config, "init.dominant_scale.enable",&param.init_dominant_scale)) {
			param.init_dominant_scale = 0;
		}
		if(!config_lookup_float(&config, "init.dominant_scale.amplitude",&param.per_amplitude_noise)) {
			param.per_amplitude_noise = 0.0;
		}
		if(!config_lookup_bool(&config, "init.Magnetic_field.enable",&param.init_Magnetic_field)) {
			param.init_Magnetic_field = 0;
		}
		if(!config_lookup_float(&config, "init.ABC_flow.A",&param.ABC_flow_A)) {
			param.ABC_flow_A = 0.0;
		}
		if(!config_lookup_float(&config, "init.ABC_flow.B",&param.ABC_flow_B)) {
			param.ABC_flow_B = 0.0;
		}
		if(!config_lookup_float(&config, "init.ABC_flow.C",&param.ABC_flow_C)) {
			param.ABC_flow_C = 0.0;
		}
		if(!config_lookup_int(&config, "init.ABC_flow.kx", &param.ABC_flow_kx)) {
			param.ABC_flow_kx = 1;
		}
		if(!config_lookup_int(&config, "init.ABC_flow.ky", &param.ABC_flow_ky)) {
			param.ABC_flow_ky = 1;
		}
		if(!config_lookup_int(&config, "init.ABC_flow.kz", &param.ABC_flow_kz)) {
			param.ABC_flow_kz = 1;
		}
        if(!config_lookup_float(&config, "init.modified_ABC_flow.A",&param.modified_ABC_flow_A)) {
			param.modified_ABC_flow_A = 0.0;
		}
		if(!config_lookup_float(&config, "init.modified_ABC_flow.B",&param.modified_ABC_flow_B)) {
			param.modified_ABC_flow_B = 0.0;
		}
		if(!config_lookup_float(&config, "init.modified_ABC_flow.C",&param.modified_ABC_flow_C)) {
    	    param.modified_ABC_flow_C = 0.0;
		}
        if(!config_lookup_float(&config, "init.modified_ABC_flow.D",&param.modified_ABC_flow_D)) {
			param.modified_ABC_flow_D = 0.0;
		}
		if(!config_lookup_int(&config, "init.modified_ABC_flow.kx", &param.modified_ABC_flow_kx)) {
			param.modified_ABC_flow_kx = 1;
		}
		if(!config_lookup_int(&config, "init.modified_ABC_flow.ky", &param.modified_ABC_flow_ky)) {
			param.modified_ABC_flow_ky = 1;
		}
		if(!config_lookup_int(&config, "init.modified_ABC_flow.kz", &param.modified_ABC_flow_kz)) {
			param.modified_ABC_flow_kz = 1;
        }
        if(!config_lookup_int(&config, "init.modified_ABC_flow.m", &param.modified_ABC_flow_m)) {
			param.modified_ABC_flow_m = 1;
        }
		if(!config_lookup_bool(&config, "init.mean_field.enable",&param.init_mean_field)) {
			param.init_mean_field = 0;
		}
		if(!config_lookup_float(&config, "init.mean_field.bx0",&param.bx0)) {
			param.bx0 = 0.0;
		}
		if(!config_lookup_float(&config, "init.mean_field.by0",&param.by0)) {
			param.by0 = 0.0;
		}
		if(!config_lookup_float(&config, "init.mean_field.bz0",&param.bz0)) {
			param.bz0 = 0.0;
		}
        if(!config_lookup_bool(&config, "init.No_Noise",&param.init_No_Noise)) {
			param.init_No_Noise = 0;
		}     
        if(!config_lookup_bool(&config, "init.dump",&param.init_dump)) {
			param.init_dump = 0;
		}
		if(!config_lookup_bool(&config, "init.bench",&param.init_bench)) {
			param.init_bench = 0;
		}
		config_destroy(&config);
	}
#ifdef MPI_SUPPORT
	MPI_Bcast( &param, sizeof(struct Parameters), MPI_CHAR, 0, MPI_COMM_WORLD);
	
	// Copy varname structures properly (Broadcast does not work because of the allocation structure we use)
	if(rank !=0 ) {
		// Allocate the name list
		param.timevar_vars.name = malloc( param.timevar_vars.length * sizeof(char*) );
	}
	
	// Next, allocate each name and copy it
	for(i = 0 ; i < param.timevar_vars.length ; i++) {
		if(rank==0) n = strlen(param.timevar_vars.name[i]);
		
		// Broadcast the string length and allocate it
		MPI_Bcast( &n, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		if(rank != 0) param.timevar_vars.name[i] = malloc( sizeof(char) * (n + 1));
		
		// Broadcast the string itself
		MPI_Bcast( param.timevar_vars.name[i], n+1, MPI_CHAR, 0, MPI_COMM_WORLD);
		
	}
	
#endif
		
	DEBUG_END_FUNC;
	
	return;
}
	
