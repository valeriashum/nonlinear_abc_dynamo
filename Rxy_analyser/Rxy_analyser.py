#! /usr/bin/python

###########################################################
# Python script written by Harry Braviner in April 2013
#
# This script takes the 'timevar' output file from the
# 'SNOOPY' code (due to Geoffroy Lesur) and looks for a
# harmonic response in the x-y component of the Reynolds
# stress
#
# The script should be run with as follows:
# ./Rxy_analyser w a0 t1 t2 [inputFile]
# w - frequency of forcing mode
# a0 - amplitude of forcing mode
# (t1, t2) - interval over which Fourier integral is taken
# inputFile = 'timevar' file from snoopy
#
# You may freely modify this script as necessary
#
###########################################################

###########################################################
# Mathematical description of Rxy_analyser.py
# 
# This script extracts an amplitude and phase lag from the
# Reynolds stress by numerically calculating a pair of
# Fourier coefficients.
#
# Denote the Reynolds stress as R(t)
# Then define the integrals I_1 and I_2 as:
# I_1 is the integral of R(t)*sin(w*t) from t=t1 to t=t2
# I_2 is the integral of R(t)*cos(w*t) from t=t1 to t=t2
#
# If R(t) = Re( a_0 (i G_0 + G_1/w) exp(iwt) ),
# and t2-t1 = 2 pi n / w,
# where n is an integer, then we would have
# G_0 = - w I_1 / n pi
# G_1 = - w^2 I_2 / n pi
# FIXME - is the sign of G_1 incorrect?
#
# We numerically approximate I_1 and I_2 using the
# trapezoidal rule.
#
###########################################################

# Start of script parameters
# These may need to be changed to match your local
# environment

# Default location of the 'timevar' file
inputFileName = "./timevar"

# End of script parameters

from sys import argv, stderr, stdout
from math import pi, sin, cos, atan2

# Parse command-line arguments
if len(argv) < 5 or len(argv) > 6:
    stderr.write("Wrong number of command-line arguments. Exiting.\n")
    exit(1)
w = float(argv[1])
a0 = float(argv[2])
t1 = float(argv[3])
t2 = float(argv[4])
if  len(argv) == 6: inputFileName = argv[5]

try:
    inputFile = open(inputFileName, "r")
except IOError:
    stderr.write("Could not open %s for reading. Exiting.\n" % inputFileName)
    exit(1)

# Determine which columns of inputFile hold time and Rxy
line = inputFile.readline()
headerCount, tColumn, RxyColumn = 0, -1, -1
for header in line.split():
    if header == 't': tColumn = headerCount
    elif header == 'Rxy': RxyColumn = headerCount
    headerCount += 1

if tColumn == -1 or RxyColumn == -1:
    stderr.write("Failed to find either time (t) or Reynolds stress (Rxy) columns in %s. Exiting.\n" % inputFileName)
    exit(1)

# Read the time and Reynolds stress into memory
tList, RxyList = [], []
for line in inputFile:
    tList.append(float(line.split()[tColumn]))
    RxyList.append(float(line.split()[RxyColumn]))

# Truncate the lists to cover (just more than) an integer
# number of forcing periods. Then use linear interpolation
# to truncate it to exactly the required interval.
N = int(w*(t2-t1)/(2*pi))
if N < 1:
    stderr.write("The interval supplied does not cover at least once forcing period. Exiting.\n")
    exit(1)

intervalMask = [i for i,t in zip(range(0,len(tList)),tList) if t>= t1 and t<= t2]
if intervalMask[0] != 0:
    intervalMask.insert(0, intervalMask[0]-1)
if intervalMask[-1] != len(tList) - 1:
    intervalMask.append(intervalMask[-1] + 1)
tList = [tList[i] for i in intervalMask]
RxyList = [RxyList[i] for i in intervalMask]

RxyList[0] = RxyList[0] + (t1-tList[0])*(RxyList[1]-RxyList[0])/(tList[1]-tList[0])
tList[0] = t1
RxyList[-1] = RxyList[-2] + (t2-tList[-2])*(RxyList[-1]-RxyList[-2])/(tList[-1]-tList[-2])
tList[-1] = t2

# Numerically compute the integrals I_1 and I_2 using
# linear interpolation
I_1, I_2 = 0, 0
for i in range(0,len(tList)-1):
    I_1 += 0.5*(tList[i+1]-tList[i])*(RxyList[i]*sin(w*tList[i]) + RxyList[i+1]*sin(w*tList[i+1]))
    I_2 += 0.5*(tList[i+1]-tList[i])*(RxyList[i]*cos(w*tList[i]) + RxyList[i+1]*cos(w*tList[i+1]))

# Calculate and output G_0 and G_1
G_0 = -w*I_1/(N*pi*a0)
G_1 = -w*w*I_2/(N*pi*a0)

stdout.write("%f\t%f\n" % (G_0, G_1))

# Exit normally
exit(0)
