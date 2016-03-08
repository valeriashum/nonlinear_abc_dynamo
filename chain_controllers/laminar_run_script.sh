#! /bin/bash

this_label=1	# Used to index runs

## Default run parameters
t_final=1000.0
osc_shear_amp=0.02
osc_shear_freq=30.0
resolution=16
offset=1	# Part of the random seed, doesn't really matter for the laminar flows
reynolds=10.0
ABCA=1.0
ABCB=1.0
ABCC=1.0

set_temp_file_and_run () {
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
	sed -i -r "s/(#define[ \t]*RAND_OFFSET[ \t]*)[0-9]+/\1$offset/" src/problem/OSABC_temp/gvars.h
	sed -i -r "s/(#define[ \t]*N[X,Y,Z][ \t]*)[0-9]+/\1$resolution/" src/problem/OSABC_temp/gvars.h
	# Write data about the run to a file so I know what we actually ran!
	echo "t_final:    $t_final"         >> laminarOutput/run$this_label
	echo "a0:         $osc_shear_amp"   >> laminarOutput/run$this_label
	echo "omega:      $osc_shear_freq"  >> laminarOutput/run$this_label
	echo "resolution: $resolution"      >> laminarOutput/run$this_label
	echo "Re:         $reynolds"        >> laminarOutput/run$this_label
	echo "A:          $ABCA"        >> laminarOutput/run$this_label
	echo "B:          $ABCB"        >> laminarOutput/run$this_label
	echo "C:          $ABCC"        >> laminarOutput/run$this_label
	# Configure, compile, and run
	echo "Configuring..."
	./configure --with-problem=OSABC_temp --enable-mpi 1> /dev/null
	make 1> /dev/null 2> /dev/null
	echo "Running..."
	mpirun -n 4 ./snoopy 1> /dev/null
	# Save the output
	mv timevar laminarOutput/timevar$this_label.dat
	# Fit gradients to the Rxy curves in timevar, echo output to a file
	Cbeta1=$(run_stats/single_run_stats.py laminarOutput/timevar$this_label.dat --start 50 --C)
	Sbeta1=$(run_stats/single_run_stats.py laminarOutput/timevar$this_label.dat --start 50 --S)
	echo "$osc_shear_amp   $osc_shear_freq   $resolution   $reynolds   $Cbeta1   $Sbeta1   $ABCA   $ABCB   $ABCC" >> laminarOutput/gradList.dat
	let "this_label = this_label + 1"
}

## Clear the directory that we'll use for storing output
mkdir laminarOutput
rm laminarOutput/*
## Wrtie header lines for gradient file
echo "#a0   omega   resolution   Re   Cbeta1   Sbeta1   ABCA   ABCB   ABCC" >> laminarOutput/gradList.dat

## Run code for a0=0.2
osc_shear_amp=0.02
resolution=16
offset=1	# Part of the random seed, doesn't really matter for the laminar flows
reynolds=10.0

osc_shear_freq=5.0
for i in `seq 1 20`; do
	set_temp_file_and_run
	osc_shear_freq=$(echo "$osc_shear_freq + 5.0" | bc)
done
osc_shear_freq=110.0
for i in `seq 1 10`; do
	set_temp_file_and_run
	osc_shear_freq=$(echo "$osc_shear_freq + 10.0" | bc)
done


## Run code for a0=0.4 (fewer runs)
osc_shear_amp=0.04
resolution=16
offset=1	# Part of the random seed, doesn't really matter for the laminar flows
reynolds=10.0

osc_shear_freq=10.0
for i in `seq 1 10`; do
	set_temp_file_and_run
	osc_shear_freq=$(echo "$osc_shear_freq + 10.0" | bc)
done
osc_shear_freq=130.0
for i in `seq 1 3`; do
	set_temp_file_and_run
	osc_shear_freq=$(echo "$osc_shear_freq + 30.0" | bc)
done


## Run code for Re = 20 (fewer runs)
osc_shear_amp=0.02
resolution=16
offset=1	# Part of the random seed, doesn't really matter for the laminar flows
reynolds=20.0

osc_shear_freq=5.0
for i in `seq 1 10`; do
	set_temp_file_and_run
	osc_shear_freq=$(echo "$osc_shear_freq + 10.0" | bc)
done
osc_shear_freq=120.0
for i in `seq 1 3`; do
	set_temp_file_and_run
	osc_shear_freq=$(echo "$osc_shear_freq + 30.0" | bc)
done


## Run code for resolution of 32
osc_shear_amp=0.02
resolution=32
offset=1	# Part of the random seed, doesn't really matter for the laminar flows
reynolds=20.0

osc_shear_freq=15.0
for i in `seq 1 4`; do
	set_temp_file_and_run
	osc_shear_freq=$(echo "$osc_shear_freq + 15.0" | bc)
done
osc_shear_freq=140.0
set_temp_file_and_run

## Run code for C = 0.5
osc_shear_amp=0.02
ABCC=0.5
resolution=16
offset=1	# Part of the random seed, doesn't really matter for the laminar flows
reynolds=10.0

osc_shear_freq=5.0
for i in `seq 1 20`; do
	set_temp_file_and_run
	osc_shear_freq=$(echo "$osc_shear_freq + 5.0" | bc)
done
osc_shear_freq=110.0
for i in `seq 1 10`; do
	set_temp_file_and_run
	osc_shear_freq=$(echo "$osc_shear_freq + 10.0" | bc)
done

## Run code for C = 0.0
osc_shear_amp=0.02
ABCC=0.0
resolution=16
offset=1	# Part of the random seed, doesn't really matter for the laminar flows
reynolds=10.0

osc_shear_freq=5.0
for i in `seq 1 20`; do
	set_temp_file_and_run
	osc_shear_freq=$(echo "$osc_shear_freq + 5.0" | bc)
done
osc_shear_freq=110.0
for i in `seq 1 10`; do
	set_temp_file_and_run
	osc_shear_freq=$(echo "$osc_shear_freq + 10.0" | bc)
done

## Run code for C = 2.0
osc_shear_amp=0.02
ABCC=2.0
resolution=16
offset=1	# Part of the random seed, doesn't really matter for the laminar flows
reynolds=10.0

osc_shear_freq=5.0
for i in `seq 1 20`; do
	set_temp_file_and_run
	osc_shear_freq=$(echo "$osc_shear_freq + 5.0" | bc)
done
osc_shear_freq=110.0
for i in `seq 1 10`; do
	set_temp_file_and_run
	osc_shear_freq=$(echo "$osc_shear_freq + 10.0" | bc)
done

echo "Finished!"
