#! /usr/bin/python

"""
Code written by Harry Braviner 2014-05-14
Intended to take the R_xy_Fourier.dat files produced by SNOOPY
when the SHELL_R_XY option is defined, fit a linear regression line
for each |k| bin, and plot nu_e against |k|.
"""

import sys
import argparse

# Start of statistical function definitions
def mean(list):
    return sum(list) / len(list)

def variance(list):
    squaredList = [x*x for x in list]
    return sum(squaredList)/len(list) - (sum(list)/len(list))**2

def covariance(list1, list2):
    if len(list1) != len(list2):
        return 0
    else:
        return sum([x*y for x,y in zip(list1, list2)]) / len(list1) -  sum(list1)*sum(list2)/(len(list1)**2)
# End of statistical function definitions

parser = argparse.ArgumentParser(description = 'Produce spectrum of nu_E from R_xy_Fourier.dat file')
parser.add_argument('inputFilename', metavar = 'R_xy_Fourier.dat', type = str, help = 'input file, typically outputted by SNOOPY to R_xy_Fourier.dat')
parser.add_argument('--start', metavar='start time', type=float, help='time from which to start fitting the regression line')
parser.add_argument('--end', metavar='end time', type=float, help='time up to which we should fir the regression line')
parser.add_argument('--verbose', action='store_true', help = 'verbose mode')
parser.add_argument('--output', type=str, default ='shell_gradients.dat' , metavar='outputFilename', help='The file to which the gradients will be written')

args = parser.parse_args()
verbose = args.verbose
start = args.start; end = args.end

try:
	inputFile = open(args.inputFilename, 'r')
except IOError:
	print 'Unable to open ' + args.inputFilename + ' for reading. Exitting...'
	sys.exit(2)

# Count the number of k bins and record the nominal k value for each of them
line = inputFile.readline().strip().split()
if verbose and line[0] != '#t':
	print "Hmm, format at start of first line doesn't look quite right. Continuing anyway"
bin_k_values = []
for k in line[1:]:
	bin_k_values.append(float(k))

if len(bin_k_values)%2 == 1:
	print "Odd number of nominal k values, Exitting..."
	sys.exit(2)

for k1, k2 in zip(bin_k_values[::2], bin_k_values[1::2]):
	if k1 != k2:
		print "Bin for k = " + str(k1) + " does not have a matching value! Exitting..."
		sys.exit(2)

# Ok, so now we've established that half of them are duplicates we can throw them away...
bin_k_values = bin_k_values[::2]

# Table to hold the raw gradients from the linear fits
intCRxy_grads = []
intSRxy_grads = []

if verbose:
	print "Found " + str(len(bin_k_values)) + " k bins, nominal k values from " + str(bin_k_values[0]) + " to " + str(bin_k_values[-1])

# Loop over the k bins, finding regression lines for the effective viscosity
for col, k in zip(range(2, 2*len(bin_k_values) + 1, 2), bin_k_values):
	t = []; intCRxy = []
	inputFile.seek(0); inputFile.readline()	# Throw away the line with the k bins
	for line in inputFile:
		line = line.strip().split()
		t.append(float(line[0]))
		intCRxy.append(float(line[col]))
	# Cut out any points not within our time limits
	if start == None:
		start = 0
	if end ==None:
		end = t[-1]
	intCRxy = [j for i,j in zip(t, intCRxy) if (i >= start and i <= end)]
	t = [i for i in t if (i >= start and i <= end)]
	# Fit the regression line
	tVar = variance(t)
	tMean = mean(t)
	intCRxyMean = mean(intCRxy)
	intCRxyVar = variance(intCRxy)
	tintCRxyCov = covariance(t, intCRxy)
	w1 = tintCRxyCov / tVar
	w2 = intCRxyMean - tMean * tintCRxyCov / tVar
	if verbose:
		print "For (1-indexed) column " + str(col + 1) + ", bin k = " + str(k) + ", found regression line for intCRxy:"
		print str(w1) + " *x + " + str(w2)
	intCRxy_grads.append(w1)

# Loop over the k bins, finding regression lines for the effective elasticity
for col, k in zip(range(1, 2*len(bin_k_values), 2), bin_k_values):
	t = []; intSRxy = []
	inputFile.seek(0); inputFile.readline()	# Throw away the line with the k bins
	for line in inputFile:
		line = line.strip().split()
		t.append(float(line[0]))
		intSRxy.append(float(line[col]))
	# Cut out any points not within our time limits
	if start == None:
		start = 0
	if end ==None:
		end = t[-1]
	intSRxy = [j for i,j in zip(t, intSRxy) if (i >= start and i <= end)]
	t = [i for i in t if (i >= start and i <= end)]
	# Fit the regression line
	tVar = variance(t)
	tMean = mean(t)
	intSRxyMean = mean(intSRxy)
	intSRxyVar = variance(intSRxy)
	tintSRxyCov = covariance(t, intSRxy)
	w1 = tintSRxyCov / tVar
	w2 = intSRxyMean - tMean * tintSRxyCov / tVar
	if verbose:
		print "For (1-indexed) column " + str(col + 1) + ", bin k = " + str(k) + ", found regression line for intSRxy:"
		print str(w1) + " *x + " + str(w2)
	intSRxy_grads.append(w1)

# Tidy up (ie. close the input file)
inputFile.close()

# Write the gradients to a file
try:
	outputFile = open(args.output, "w")
except IOError:
	print 'Unable to open ' + args.output + ' for writing. Exitting...'
	sys.exit(2)

outputFile.write("#k\tintCRxy_grad\tintSRxy_grad\n")

for (k, C, S) in zip(bin_k_values, intCRxy_grads, intSRxy_grads):
	outputFile.write(str(k) + "\t" + str(C) + "\t" + str(S) + "\n")

outputFile.close()
