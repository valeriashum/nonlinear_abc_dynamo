#! /usr/bin/python2.7

###### Shear / No Shear Comparison Runscript ######
##  Written by Harry Braviner 2015-02-02
## The purpose of this script is to automate succesive
## runnings of the snoopy code, alternating between shearing
## and no-shear in the OSABC problem.
##
## The reason for doing so is to try and beat down the noise
## in the xy component of the Reynolds stress. We shall do so
## by:
## 1) Running the code to time=50 with no shear, to allow decay
##    from the ABC flow to a turbulent state.
## 2) Setting t0 <- time. Save the dump file as dump1.
##    Run the code with no shear up to t = t0 + t_switch.
## 3) Save the dump file as dump2. Restore dump1 as the dump file.
## 4) Re-run the code from t0 to t0 + t_switch with shear switched on.
## 5) Calculate the difference of
##    int[t0, t0+t_switch] of R_xy cos(omega*T) between the two runs.

import re
import sys, os
from subprocess import call
from math import pi

omega = 50.0     # Shear angular frequency
a0    = 0.04    # Shear amplitude
t_turb = 50.0		# Time to run the simulation for to let turbulence fully develop
diff_from_shear = False # If true, we restart from sheared simulation, else from unsheared
t_switch = 15.0      # Time to run for before switching between sheared and unsheared
t_final = 500.0     # Time to finish at
delta_diff = 0.05   # Resolution for the differences file

## Convert t_turb, t_switch and t_final to be an integer number of oscillation periods
t_turb = t_turb - (t_turb%(2.0*pi/omega))
t_switch = t_switch - (t_switch%(2.0*pi/omega))
t_final = t_final - (t_final%(2.0*pi/omega))

## Name of config files
snoopyDir = "/home/hb295/snoopy6/"
scriptDir = os.path.dirname(os.path.realpath(__file__))
configFilename = snoopyDir + "src/problem/OSABC/snoopy.cfg"
#configFilename = snoopyDir + "./snoopy.cfg"
gvarsFilename = snoopyDir + "src/problem/OSABC/gvars.h"
savedDumpFilename = scriptDir + "/" + "saved_dump.dmp"
tempDumpFilename = scriptDir + "/" + "temp_dump.dmp"
shearRxyFilename = scriptDir + "/" + "shearRxy.dat"
diffRxyFilename = scriptDir + "/" + "diffRxy.dat"
noshearRxyFilename = scriptDir + "/" + "noshearRxy.dat"
timevarFilename = snoopyDir + "timevar"

def enableShear():
    regex1 = r'^\s*oscillatory_shear_amp\s*=\s*[0-9.]+\s*;$'
    repl1  = '    oscillatory_shear_amp = ' + str(a0) + ';'
    regex2 = r'^\s*oscillatory_shear_freq\s*=\s*[0-9.]+\s*;$'
    repl2  = '    oscillatory_shear_freq = ' + str(omega) + ';'
    # Read in the config file
    with open(configFilename, "r") as configFile:
        configLines = configFile.readlines()

    # Write it out, replacing lines as we go
    with open(configFilename, "w") as configFile:
        for line in configLines:
            configFile.write(re.sub(regex1, repl1,
                                    re.sub(regex2, repl2, line)))

def disableShear():
    regex1 = r'^\s*oscillatory_shear_amp\s*=\s*[0-9.]+\s*;$'
    repl1  = '    oscillatory_shear_amp = ' + '0' + ';'
    regex2 = r'^\s*oscillatory_shear_freq\s*=\s*[0-9.]+\s*;$'
    repl2  = '    oscillatory_shear_freq = ' + str(omega) + ';'
    # Read in the config file
    with open(configFilename, "r") as configFile:
        configLines = configFile.readlines()

    # Write it out, replacing lines as we go
    with open(configFilename, "w") as configFile:
        for line in configLines:
            configFile.write(re.sub(regex1, repl1,
                                    re.sub(regex2, repl2, line)))

def settfinal(tf):
    regex1 = '^\s*t_final\s*=\s*[0-9.]+.*'
    repl1 = '	t_final = ' + str(tf) + ';							// Simulation will stop if it reaches this time'
    # Read in the config file
    with open(configFilename, "r") as configFile:
        configLines = configFile.readlines()

    # Write it out, replacing lines as we go
    with open(configFilename, "w") as configFile:
        for line in configLines:
            configFile.write(re.sub(regex1, repl1, line))

def enableRestart():
    regex1 = r'^\s*restart\s*=\s*.*'
    repl1 = '\trestart = true;						// set to true to restart from a dump file. If no dump file is found, this option has no effect.'
    # Read in the config file
    with open(configFilename, "r") as configFile:
        configLines = configFile.readlines()

    # Write it out, replacing lines as we go
    with open(configFilename, "w") as configFile:
        for line in configLines:
            configFile.write(re.sub(regex1, repl1, line))

def disableRestart():
    regex1 = r'^\s*restart\s*=\s*.*'
    repl1 = '\trestart = false;						// set to true to restart from a dump file. If no dump file is found, this option has no effect.'
    # Read in the config file
    with open(configFilename, "r") as configFile:
        configLines = configFile.readlines()

    # Write it out, replacing lines as we go
    with open(configFilename, "w") as configFile:
        for line in configLines:
            configFile.write(re.sub(regex1, repl1, line))

def cleanSnoopy():
    FNULL = open(os.devnull, 'w')
    call(["rm *.dat time* ./data/*"],
         stdout=FNULL, stderr=FNULL, cwd=snoopyDir, shell=True)

def cleanScript():
    FNULL = open(os.devnull, 'w')
    call(["rm *.dmp *.dat"],
         stdout=FNULL, stderr=FNULL, cwd=scriptDir, shell=True)

def configureAndMake():
    try:
        FNULL = open(os.devnull, 'w')
        retcode = call(["./configure", "--enable-mpi", "--with-problem=OSABC"],
                       cwd=snoopyDir, stdout = FNULL, stderr = FNULL)
        if retcode != 0:
            print "configure returned a non-zero exit code!"
            sys.exit()
        retcode = call(["make"], cwd=snoopyDir, stdout = FNULL, stderr = FNULL)
        if retcode != 0:
            print "make returned a non-zero exit code!"
            sys.exit()
    except OSError as e:
        print "Configure and Make failed: ", e
        sys.exit()
	
def runSnoopy():
    try:
        FNULL = open(os.devnull, 'w')
        retcode = call(["mpirun", "-n", "4", "./snoopy"], stdout=FNULL,
                        stderr=FNULL, cwd=snoopyDir)
        if retcode != 0:
            print "Running snoopy returned a non-zero exit code!"
            sys.exit()
    except OSError as e:
        print "Running Snoopy failed: ", e
        sys.exit()

def saveDump():
    try:
        FNULL = open(os.devnull, 'w')
        retcode = call(["cp", "./dump.dmp", savedDumpFilename], cwd=snoopyDir,
                       stdout=FNULL, stderr=FNULL)
    except OSError as e:
        print "Saving the dump file failed: ", e
        sys.exit()

def restoreDump():
    try:
        FNULL = open(os.devnull, 'w')
        retcode = call(["cp", savedDumpFilename, "./dump.dmp"], cwd=snoopyDir,
                       stdout=FNULL, stderr=FNULL)
    except OSError as e:
        print "Saving the dump file failed: ", e
        sys.exit()

def exchangeDumps():
    try:
        FNULL = open(os.devnull, 'w')
        retcode = call(["cp", savedDumpFilename, tempDumpFilename], cwd=snoopyDir,
                       stdout=FNULL, stderr=FNULL)
        retcode = call(["cp", "./dump.dat", savedDumpFilename], cwd=snoopyDir,
                       stdout=FNULL, stderr=FNULL)
        retcode = call(["cp", tempDumpFilename, "./dump.dmp"], cwd=snoopyDir,
                       stdout=FNULL, stderr=FNULL)
        retcode = call(["rm", tempDumpFilename], cwd=snoopyDir,
                       stdout=FNULL, stderr=FNULL)
    except OSError as e:
        print "Saving the dump file failed: ", e
        sys.exit()

def appendRxy(shearRxyFilename):
    # If the file doesn't exist, write the header
    if not(os.path.isfile(shearRxyFilename)):
        outfile = open(shearRxyFilename, "w")
        outfile.write("#t\tRxy\tintSRxy\tintCRxy\tdisp_rate\n")
        with open(timevarFilename) as infile:
            inlines = infile.readlines()
        for line in inlines[1:]:
            sl = line.strip().split()
            outfile.write(sl[0] + "\t" + sl[9] + "\t" + sl[10] + "\t"
                          + sl[11] + "\t" + sl[12] + "\n")
        outfile.close()
    else:
        # Get the earliest time from the timevar file
        infile = open(timevarFilename, "r")
        line = infile.readline()
        infile.close()
        t0 = float(line.strip().split()[0])
        # Discard any of the outputfile with a time greater than this
        with open(shearRxyFilename, "r") as outfile:
            lines = outfile.readlines()
        with open(shearRxyFilename, "w") as outfile:
            outfile.write(lines[0])
            for line in lines[1:]:
                if float(line.strip().split()[0]) < t0:
                    outfile.write(line)
        # Get the values of the integrals so far (they are
        # lost on the restart)
        outfile = open(shearRxyFilename, "r")
        oldlines = outfile.readlines()
        outfile.close()
        intSR0 = float(oldlines[-1].strip().split()[2])
        intCR0 = float(oldlines[-1].strip().split()[3])
        outfile = open(shearRxyFilename, "a")
        with open(timevarFilename) as infile:
            inlines = infile.readlines()
        for line in inlines:
            sl = line.strip().split()
            intSR = str(float(sl[10]) + intSR0)
            intCR = str(float(sl[11]) + intCR0)
            outfile.write(sl[0] + "\t" + sl[9] + "\t" + intSR + "\t"
                          + intCR + "\t" + sl[12] + "\n")
        outfile.close()

def appendNoShearRxy():
    appendRxy(noshearRxyFilename)

def appendShearRxy():
    appendRxy(shearRxyFilename)

def constructDiffRxy():
    diffFile = open(diffRxyFilename, "w")
    diffFile.write("#t\tintSRxy\tintSCxy\n")
    with open(shearRxyFilename, "r") as shearFile:
        shearLines = shearFile.readlines()
    with open(noshearRxyFilename, "r") as noshearFile:
        noshearLines = noshearFile.readlines()
    tmax = min(float(shearLines[-2].strip().split()[0]),
               float(noshearLines[-2].strip().split()[0]))
    tmin = max(float(shearLines[2].strip().split()[0]),
               float(noshearLines[2].strip().split()[0]))
    t = tmin
    si = 1
    nsi = 1
    while t <= tmax:
        while float(shearLines[si].strip().split()[0]) < t:
           si += 1
           ta = float(shearLines[si-1].strip().split()[0])
           tb = float(shearLines[si].strip().split()[0])
           a = float(shearLines[si-1].strip().split()[2])
           b = float(shearLines[si].strip().split()[2])
           sintSRxy = a + (t - ta)*(b-a)/(tb-ta)
           a = float(shearLines[si-1].strip().split()[3])
           b = float(shearLines[si].strip().split()[3])
           sintCRxy = a + (t - ta)*(b-a)/(tb-ta)
        while float(noshearLines[nsi].strip().split()[0]) < t:
           nsi += 1
           ta = float(noshearLines[nsi-1].strip().split()[0])
           tb = float(noshearLines[nsi].strip().split()[0])
           a = float(noshearLines[nsi-1].strip().split()[2])
           b = float(noshearLines[nsi].strip().split()[2])
           nsintSRxy = a + (t - ta)*(b-a)/(tb-ta)
           a = float(noshearLines[nsi-1].strip().split()[3])
           b = float(noshearLines[nsi].strip().split()[3])
           nsintCRxy = a + (t - ta)*(b-a)/(tb-ta)
        diffFile.write(str(t) + "\t" + str(sintSRxy - nsintSRxy) + "\t" + str(sintCRxy - nsintCRxy) + "\n")
        t += delta_diff
    diffFile.close()
    

def runTot_turb():
    cleanScript()
    cleanSnoopy()
    disableRestart()
    if diff_from_shear: disableShear()
    else: enableShear()
    settfinal(t_turb)
    configureAndMake()
    runSnoopy()
    saveDump()
	
def runWithShear(t):
    cleanSnoopy()
    enableRestart()
    enableShear()
    restoreDump()
    settfinal(t+t_switch)
    configureAndMake()
    runSnoopy()
    if diff_from_shear: exchangeDumps()
    else: restoreDump()
    appendShearRxy()

def runWithoutShear(t):
    cleanSnoopy()
    enableRestart()
    disableShear()
    settfinal(t+t_switch)
    configureAndMake()
    runSnoopy()
    if not(diff_from_shear): saveDump()
    appendNoShearRxy()

# Not actually doing anything with gvars in this version of the code
#    with open(gvarsFilename, "r") as gvarsFile:
#        gvarsLines = gvarsFile.readlines()

#enableShear(configFilename, gvarsFilename)
#configureAndMake()
#settfinal(753)

def runAlternate():
   t = 0
   cleanScript()
   print "Starting run from time 0 to time " + str(t_turb)
   runTot_turb()
   t += t_turb
   while t < t_final:
       print "Starting sheared run from time " + str(t) + " to time " + str(t + t_switch)
       runWithShear(t)
       print "Starting unsheared run from time " + str(t) + " to time " + str(t + t_switch)
       runWithoutShear(t)
       t += t_switch
   constructDiffRxy()
   print "Done."

runAlternate()
