#! /bin/bash
# Constants
pi=3.14159265359 
number_of_cores=8 #(min=4 for the ABC flow if resolution=number_of_cores)
advection_time=14.0

# Set in gvarsh.h
N_x=32          #Number of 2pi boxes in direction (x)
N_y=2           #Number of 2pi boxes in direction (y)
N_z=1           #Number of 2pi boxes in direction (y)
Lx=$(echo "$N_x*2*$pi" |bc)
Ly=$(echo "$N_y*2*$pi" |bc)
Lz=$(echo "$N_z*2*$pi" |bc)

# Parameters for list of runs, may change
this_label=1	                                    # Used to index runs
t_final=1000000.0                                      # Time to run to, in program units
ABCA=1.0                                            # Velocity magnitudes
ABCB=$(echo "$ABCA" |bc)                            # Velocity magnitudes
ABCC=$(echo "$ABCA" |bc)                            # Velocity magnitudes
ABCD=1.0                                            # Velocity magnitudes
kux=1                                               #  velocity k_x
kuy=1                                               #  velocity k_y
kuz=1                                               #  velocity k_z
m=$(echo "$N_x" |bc)                                #  modified m
timevarstep=1.0		                            # Time between two outputs in the timevar file
snapshotstep=1.0		                            # Time between two snapshot outputs
dumpstep=20.0			                            # Time between two restart dump outputs (restart dump are erased)
restar="false"			                            # set to true to restart from a dump file. 

# Define the subroutine to create a problem folder and run 
# the code for various parameters 
set_temp_file_and_run () {
    # Clean up the diretory from any previous runs
	rm -rf src/problem/kinematic_ABC_flow_temp

    # Build a diretory containing the modified gvars.h and snoopy.cfg
	mkdir -p src/problem/kinematic_ABC_flow_temp
	cp src/problem/kinematic_dynamo_box_10x2x1/gvars.h    src/problem/kinematic_ABC_flow_temp/
	cp src/problem/kinematic_dynamo_box_10x2x1/snoopy.cfg src/problem/kinematic_ABC_flow_temp/
    
    # Modify the simulation settings in the temp folder
	sed -i -r "s/([ \t]*timevar_step[ \t]*=[ \t]*)[0-9,.]+/\1$timevarstep/"         src/problem/kinematic_ABC_flow_temp/snoopy.cfg
	sed -i -r "s/([ \t]*snapshot_step[ \t]*=[ \t]*)[0-9,.]+/\1$snapshotstep/"       src/problem/kinematic_ABC_flow_temp/snoopy.cfg
    sed -i -r "s/([ \t]*dump_step[ \t]*=[ \t]*)[0-9,.]+/\1$dumpstep/"               src/problem/kinematic_ABC_flow_temp/snoopy.cfg
    sed -i -r "s/([ \t]*restart[ \t]*=)[ \t]*[a-z]+/\1$restar/"                     src/problem/kinematic_ABC_flow_temp/snoopy.cfg
    sed -i -r "s/([ \t]*reynolds_magnetic[ \t]*=[ \t]*)[0-9,.]+/\1$reynolds_m/"     src/problem/kinematic_ABC_flow_temp/snoopy.cfg
	sed -i -r "s/([ \t]*t_final[ \t]*=[ \t]*)[0-9,.]+/\1$t_final/"                  src/problem/kinematic_ABC_flow_temp/snoopy.cfg
	sed -i -r "s/([ \t]*A[ \t]*=[ \t]*)[0-9,.]+/\1$ABCA/"                           src/problem/kinematic_ABC_flow_temp/snoopy.cfg
	sed -i -r "s/([ \t]*B[ \t]*=[ \t]*)[0-9,.]+/\1$ABCB/"                           src/problem/kinematic_ABC_flow_temp/snoopy.cfg
	sed -i -r "s/([ \t]*C[ \t]*=[ \t]*)[0-9,.]+/\1$ABCC/"                           src/problem/kinematic_ABC_flow_temp/snoopy.cfg
    sed -i -r "s/([ \t]*D[ \t]*=[ \t]*)[0-9,.]+/\1$ABCD/"                           src/problem/kinematic_ABC_flow_temp/snoopy.cfg
    sed -i -r "s/([ \t]*m[ \t]*=[ \t]*)[0-9,.]+/\1$m/"                              src/problem/kinematic_ABC_flow_temp/snoopy.cfg
    sed -i -r "s/([ \t]*kx[ \t]*=[ \t]*)[0-9,.]+/\1$kux/"                           src/problem/kinematic_ABC_flow_temp/snoopy.cfg
    sed -i -r "s/([ \t]*ky[ \t]*=[ \t]*)[0-9,.]+/\1$kuy/"                           src/problem/kinematic_ABC_flow_temp/snoopy.cfg
    sed -i -r "s/([ \t]*kz[ \t]*=[ \t]*)[0-9,.]+/\1$kuz/"                           src/problem/kinematic_ABC_flow_temp/snoopy.cfg
	sed -i -r "s/(#define[ \t]*NX[ \t]*)[0-9]+/\1$resolution_x/"                    src/problem/kinematic_ABC_flow_temp/gvars.h
	sed -i -r "s/(#define[ \t]*NY[ \t]*)[0-9]+/\1$resolution_y/"                    src/problem/kinematic_ABC_flow_temp/gvars.h
	sed -i -r "s/(#define[ \t]*NZ[ \t]*)[0-9]+/\1$resolution_z/"                    src/problem/kinematic_ABC_flow_temp/gvars.h


	# Write data about the run to a file so I know what we actually ran!
	echo "t_final:    $t_final"         >> kinematicOutput/run$this_label
	echo "dt_timevar: $timevarstep"    >> kinematicOutput/run$this_label
	echo "dt_dump:    $dumpstep"        >> kinematicOutput/run$this_label    
	echo "Rm:         $reynolds_m"      >> kinematicOutput/run$this_label
    echo "Restart:    $restar"          >> kinematicOutput/run$this_label
	echo "A:          $ABCA"            >> kinematicOutput/run$this_label
	echo "B:          $ABCB"            >> kinematicOutput/run$this_label
	echo "C:          $ABCC"            >> kinematicOutput/run$this_label
    echo "NX:         $resolution_x"    >> kinematicOutput/run$this_label
    echo "NY:         $resolution_y"    >> kinematicOutput/run$this_label
    echo "NZ:         $resolution_z"    >> kinematicOutput/run$this_label

	# Configure, compile, and run
	echo "Configuring..."
    make clean                                                      1> /dev/null 2> /dev/null
	./configure --with-problem=kinematic_ABC_flow_temp --enable-mpi >> kinematicOutput/make$this_label      #1> /dev/null
	make                                                            1> kinematicOutput/make$this_label  2> kinematicOutput/make$this_label #1> /dev/null 2> /dev/null
    echo "Running..."
	
    mpirun -n $number_of_cores ./snoopy                                            1> /dev/null  >> kinematicOutput/summary$this_label
	# Save the output
    mkdir -p kinematicOutput/rates$this_label
    cp rate* kinematicOutput/rates$this_label
    cp -r data kinematicOutput/data$this_label
    rm data/*
	cp timevar kinematicOutput/timevar$this_label.dat
    cp spectrum.dat kinematicOutput/spectrum$this_label.dat
    cp dump.dmp kinematicOutput/dump$this_label.dmp
    cp dump_sav.dmp kinematicOutput/dump_sav$this_label.dmp
	# Fit gradients to the Rxy curves in timevar, echo output to a file
	#Cbeta1=$(run_stats/single_run_stats.py kinematicOutput/timevar$this_label.dat --start 50 --C)
	#Sbeta1=$(run_stats/single_run_stats.py kinematicOutput/timevar$this_label.dat --start 50 --S)
	#echo "$osc_shear_amp   $osc_shear_freq   $resolution   $reynolds   $Cbeta1   $Sbeta1   $ABCA   $ABCB   $ABCC" >> kinematicOutput/gradList.dat
	let "this_label = this_label + 1"
}

## Clear the directory that we'll use for storing output
##mkdir kinematicOutput
rm -rf kinematicOutput/*
rm kinematicOutput/*
rm plots/*
rm data/*
## Write header lines for gradient file
#echo "#a0   omega   resolution   Re   Cbeta1   Sbeta1   ABCA   ABCB   ABCC" >> kinematicOutput/gradList.dat


resolution_x=$(echo "$N_x*16" |bc)   # NX
resolution_y=$(echo "$N_y*16" |bc)   # NX
resolution_z=$(echo "$N_z*16" |bc)   # NZ


reynolds_m=0.153

for i in `seq 1 3`; do
    echo "Rm = $reynolds_m"
    set_temp_file_and_run
    reynolds_m=$(echo "$reynolds_m + 0.05" | bc)
done





echo "$resolution run done!"
