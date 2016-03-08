#! /bin/bash

for i in `seq 1 5`; do
	rm snoopy spectrum.dat timevar nohup.out R_xy_Fourier.dat
	./configure --with-problem=OSABC$i --enable-mpi
	make
	mpirun -n 4 ./snoopy >> output$i.txt
	cp timevar timevar$i.dat
	# Take a mean energy spectrum
	./spectrum_reader/energy_spectrum spectrum.dat 50.0 2000.0 500 time_av_spectrum$i.dat
	# Fit regression lines to the fourier integrals of Reynolds stress
	./spectrum_reader/Rxy_shell_spectrum.py R_xy_Fourier.dat --start 50 --output shell_stress_$i.dat
done

# Compute means of reynolds stress gradient spectra
./spectrum_reader/error_calc.py shell_stress_1.dat shell_stress_2.dat shell_stress_3.dat shell_stress_4.dat shell_stress_5.dat --output averaged_shell.dat
rm shell_stress_1.dat shell_stress_2.dat shell_stress_3.dat shell_stress_4.dat shell_stress_5.dat

# Compute means of energy spectra
./spectrum_reader/error_calc.py time_av_spectrum1.dat time_av_spectrum2.dat time_av_spectrum3.dat time_av_spectrum4.dat time_av_spectrum5.dat --output energy_spectrum.dat
rm time_av_spectrum1.dat time_av_spectrum2.dat time_av_spectrum3.dat time_av_spectrum4.dat time_av_spectrum5.dat
