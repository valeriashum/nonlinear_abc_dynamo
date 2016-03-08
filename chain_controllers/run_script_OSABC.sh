#! /bin/bash

# Parameters for list of runs, may change
num_runs=5
osc_shear_amp=0.02	# Amplitude to use for the oscillatory shear
osc_shear_freq=0.05	# Frequency to use for the oscillatory shear
first_seed=1		# The 'offset' random seed to use
seed_inc=4		# How many to change the seed by between runs
resolution=64		# NX, NY, NZ
reynolds=1000.0		# Reynolds number
ABCA=1.0	# Velocity magnitudes
ABCB=1.0	# Velocity magnitudes
ABCC=1.0	# Velocity magnitudes
t_final=2000.00		# Time to run to, in program units
first_label=161		# Used to give the filenames a sensible labelling for storage
# End of parameters, do not change below here!

stress_list=""
energy_list=""

for i in `seq 1 $num_runs`; do
	# Clean up the diretory from any previous runs
	rm snoopy spectrum.dat timevar R_xy_Fourier.dat
	# Build a diretory containing the modified gvars.h and snoopy.cfg
	rm -rf src/problem/OSABC_temp
	mkdir -p src/problem/OSABC_temp
	cp src/problem/OSABC/gvars.h    src/problem/OSABC_temp/
	cp src/problem/OSABC/snoopy.cfg src/problem/OSABC_temp/
	sed -i -r "s/([ \t]*oscillatory_shear_amp[ \t]*=[ \t]*)[0-9,.]+/\1$osc_shear_amp/" src/problem/OSABC_temp/snoopy.cfg
	sed -i -r "s/([ \t]*oscillatory_shear_freq[ \t]*=[ \t]*)[0-9,.]+/\1$osc_shear_freq/" src/problem/OSABC_temp/snoopy.cfg
	sed -i -r "s/([ \t]*reynolds[ \t]*=[ \t]*)[0-9,.]+/\1$reynolds/" src/problem/OSABC_temp/snoopy.cfg
	sed -i -r "s/([ \t]*t_final[ \t]*=[ \t]*)[0-9,.]+/\1$t_final/" src/problem/OSABC_temp/snoopy.cfg
	sed -i -r "s/([ \t]*A[ \t]*=[ \t]*)[0-9,.]+/\1$ABCA/" src/problem/OSABC_temp/snoopy.cfg
	sed -i -r "s/([ \t]*B[ \t]*=[ \t]*)[0-9,.]+/\1$ABCB/" src/problem/OSABC_temp/snoopy.cfg
	sed -i -r "s/([ \t]*C[ \t]*=[ \t]*)[0-9,.]+/\1$ABCC/" src/problem/OSABC_temp/snoopy.cfg
	let "this_offset = first_seed + seed_inc * i"
	sed -i -r "s/(#define[ \t]*RAND_OFFSET[ \t]*)[0-9]+/\1$this_offset/" src/problem/OSABC_temp/gvars.h
	sed -i -r "s/(#define[ \t]*N[X,Y,Z][ \t]*)[0-9]+/\1$resolution/" src/problem/OSABC_temp/gvars.h
	#sed -i -r "s/(#define[ \t]*N[X,Y,Z][ \t]*)[0-9]+/\18/" src/problem/OSABC_temp/gvars.h
	./configure --with-problem=OSABC_temp --enable-mpi
	#./configure --with-problem=OSABC_temp
	make
	let "this_label = first_label + i - 1"
	mpirun -n 4 ./snoopy >> output$this_label.txt
	#./snoopy >> output$this_label.txt
	cp timevar timevar$this_label.dat
	# Check that the energy_spectrum program is compiled
	if [ ! -e ./spectrum_reader/energy_spectrum ]
	then
		gcc -Wall ./spectrum_reader/energy_spectrum.c -o ./spectrum_reader/energy_spectrum
	fi

	# Take a mean energy spectrum
	./spectrum_reader/energy_spectrum spectrum.dat 50.0 $t_final 500 time_av_spectrum$this_label.dat
	energy_list="$energy_list time_av_spectrum$this_label.dat"
	# Fit regression lines to the fourier integrals of Reynolds stress
	./spectrum_reader/Rxy_shell_spectrum.py R_xy_Fourier.dat --start 50 --output shell_stress_$this_label.dat
	stress_list="$stress_list shell_stress_$this_label.dat"
done

let "this_label = first_label"
# Compute means of reynolds stress gradient spectra
./spectrum_reader/error_calc.py $stress_list --output averaged_shell_$this_label.dat
rm $stress_list

# Compute means of energy spectra
./spectrum_reader/error_calc.py $energy_list --output averaged_energy_$this_label.dat
rm $energy_list
