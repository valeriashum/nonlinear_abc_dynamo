#! /usr/bin/python

###########################################################
# Python script written by Harry Braviner in April 2013
#
# This script is intended to automate the configuration
# and running of a large number of runs of the 'SNOOPY'
# code (due to Geoffry Lesur)
#
# You may freely modify this script as necessary
#
###########################################################

# Start of script parameters
# These may need to be changed to match your local
# environment

# Location of the 'snoopy' and 'configure' executables:
executable_dir = "../"
# Location of the 'Rxy_analyser' executable:
RxyExecDir = "../Rxy_analyser"
# Configuration file that 'configure' will use:
comp_cfg_filename = "../src/problem/OSABC/snoopy.cfg"
# Schedule of runs to make:
runlist_filename = "./runlist.dat"
# Base directory for output (subdirectory will be made):
output_directory = "../../snoopy_output"
# File to which log of this script will be written:
log_filename = "runlog.dat"
# Base version of snoopy.cfg:
local_cfg_filename = './snoopy.cfg'

# End of script parameters

from datetime import datetime
from os import path
import os, re, subprocess, shutil

# Start of set-up

# Check whether the necessary paths exist, create them if
# the do not

configure_executable = path.abspath(path.join(executable_dir, "configure"))
snoopy_executable = path.abspath(path.join(executable_dir, "snoopy"))
Rxy_analyser_executable = path.abspath(path.join(RxyExecDir, "Rxy_analyser.py"))

if not os.access(snoopy_executable, os.X_OK) and os.access(configure_executable, os.X_OK):
    print "At least one of %s and %s either do not exist or are not executable." % (snoopy_executable, configure_executable)
    exit(1)

if not os.access(Rxy_analyser_executable, os.X_OK):
    print "%s either does not exist or is not executable." % Rxy_analyser_executable
    exit(1)

if not os.access(output_directory, os.W_OK):
    print "The output base directory %s either does not exist or you do not have write permissions for it." % output_directory
    exit(1)

outdir = path.join(output_directory, "%d_%d_%d" % (datetime.now().year, datetime.now().month, datetime.now().day))

if not path.exists(outdir):
    print "%s does not exist, creating it." % outdir
    try: os.mkdir(outdir)
    except OSError:
        print "Failed to make the directory %s, exiting." % outdir
	exit(1)

# Make a numbered sub-directory within the dated directory
i = 1
while path.exists(path.join(outdir, "automated_runs_%d" % i)) : i += 1
outdir = path.join(outdir, "automated_runs_%d" % i)
try: os.mkdir(outdir)
except OSError:
        print "Failed to make the directory %s, exiting." % path.abspath(outdir)
	exit()

# Open log-file, with no buffer
try: logfile = open(path.join(outdir, log_filename), 'w', 0)
except IOError:
    print "Failed to open log file %s for writing, exiting." % path.abspath(logfile)
    exit()

# Get the schedule of runs to perform from the file
try: runlist_file = open(runlist_filename, 'r')
except IOError:
    print "Failed to open %s for reading, exiting." % runlist_filename
    exit()
runlist = []
for line in runlist_file:
    if line.strip()[0] != '#': # Is line not a comment?
        try: runlist.append([float(line.strip().split()[0]), float(line.strip().split()[1])])
	except:
	    print "It seems that the file %s contains the line:\n%swhich is not of the form\nfloat float\nIgnoring that line...\n" % (runlist_filename, line)
runlist_file.close()
# End of set-up

###########################################################
# We should now be holding the following variables:
#
# snoopy_executable, config_executable - we know that we
#   have execution permissions on these
# outdir - a dated and numbered output directory on which
#   we have write permissions
# runlist - a list of two-entry lists FIXME - what do these hold?
# logfile - a file open for reading, to which we will write
#   updates regarding the status of the runs
###########################################################

###########################################################
# Start of main loop (over runs)
###########################################################
for shear_params in runlist:
    # Log that we are beginning the run
    logfile.write("Beginning run with a_0 = %f, w = %f on %d/%d/%d at %d:%d:%d\n" % (shear_params[0], shear_params[1], datetime.now().year, datetime.now().month, datetime.now().day, datetime.now().hour, datetime.now().minute, datetime.now().second))
    
    # Start of compiling and running snoopy

    # Update the compiling copy of snoopy.cfg based on the
    # local copy
    try:
        local_cfg_file = open(local_cfg_filename, 'r')
	comp_cfg_file = open(comp_cfg_filename, 'w')
    except IOError:
        logfile.write("There was a problem opening one of %s or %s for reading and writing respectively. Check existence and permissions. Exiting!\n" % (path.abspat(local_cfg_filename), path.abspath(comp_cfg_filename)))
	exit()
    # Line-by-line replacement, changing shear parameters
    for line in local_cfg_file:
        comp_cfg_file.write(re.sub( r"(oscillatory_shear_amp)\s*=.*", (r"\1 = %f;" % shear_params[0]), re.sub(r"(oscillatory_shear_freq)\s*=.*", (r"\1 = %f;" % shear_params[1]), line)))
    local_cfg_file.close()
    comp_cfg_file.close()

    # Run the 'configuration' script
    startdir = os.getcwd()
    os.chdir(executable_dir)
    try:
        if subprocess.call([configure_executable, "--with-problem=OSABC"]):
            logfile.write("Configuration script returned non-zero return code, configuration process may have failed.\n")
        else: logfile.write("Configuration appeared to succeed.\n")
    except OSError:
        logfile.write("Attempting to run configuration script threw an exception!\n")
    os.chdir(startdir)

    # Run 'make'
    startdir = os.getcwd()
    os.chdir(executable_dir)
    try:
        if subprocess.call("make"):
            logfile.write("make returned non-zero return code, make process may have failed.\n")
        else: logfile.write("make appeared to succeed.\n")
    except OSError:
        logfile.write("Attempting to run make threw an exception!\n")
    os.chdir(startdir)
    
    # Run snoopy
    startdir = os.getcwd()
    os.chdir(executable_dir)
    try:
        if subprocess.call(snoopy_executable):
            logfile.write("Snoopy returned a non-zero return code.\n")
	else:
	    logfile.write("SNOOPY execution appeared to succeed.\n")
    except OSError:
        logfile.write("Attempting to run snoopy threw an exception!\n")
    os.chdir(startdir)

    # End of compiling and running snoopy

    # Start of copying output files
    timevar_dest = path.abspath(path.join(outdir, "timevar_amp_%f_freq_%f.dat" % (shear_params[0], shear_params[1])))
    try:
        shutil.copyfile(path.join(executable_dir, "timevar"), path.join(outdir, timevar_dest))
    except IOError:
        logfile.write("From the run with parameters (%f, %f), copying the timevar file failed\n" % (shear_params[0], shear_params[1]))
    # End of copying output files

    # Start of running analysis programs on output files
    try:
        # FIXME - find a more intelligent way of choosing t1 and t2
        G_vals = subprocess.Popen([Rxy_analyser_executable, "%f" % shear_params[1], "%f" % shear_params[0], "%f" % 30.0, "%f" % 90.0, timevar_dest],stderr=None, stdout=subprocess.PIPE).communicate()[0]
        G_vals = [float(g) for g in G_vals.strip().split()]
    except:
        logfile.write("From the run with parameters (%f, %f), running %s failed\n" % (shear_params[0], shear_params[1], Rxy_analyser_executable))
    # If a G-file does not exist for this amp, create one
    GFileName = path.join(outdir, "Rxy_amp_%f.dat" % shear_params[0])
    if not path.exists(GFileName):
        GFile = open(GFileName, 'w')
	GFile.write("#a_0\tw\tG_0\tG_1\n")
	GFile.close()
    GFile = open(GFileName, 'a')
    GFile.write("%f\t%f\t%f\t%f\n" % (shear_params[0], shear_params[1], G_vals[0], G_vals[1]))
    GFile.close()

    # Record in the logfile that we've finished this run
    logfile.write("Finished run with a_0 = %f, w= %f on %d/%d/%d at %d:%d:%d\n" % (shear_params[0], shear_params[1], datetime.now().year, datetime.now().month, datetime.now().day, datetime.now().hour, datetime.now().minute, datetime.now().second))
    # End of running analysis programs on output files

###########################################################
# End of main loop (over runs)
###########################################################

# Start of clean-up
logfile.close()
# End of clean-up
