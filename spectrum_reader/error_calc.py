#! /usr/bin/python

"""
Program to quickly and easily calculate means and errors for
the spectra.
Take N input files of the form:
#x	y1	y2	...	ym
2.3	3.2	.3	...	0.4
...
3.9	3.1	7.	...	2.2
and produces an output file on columns of means and errors.
The 'error' is calculated as the standard deviation of the
mean, Var(mean yn) = Var(yn) / N, where Var(yn) is estimated
as the sample variance.
Different rows are handled independently from one another.
The formula above assumes that the numbers in different data
files are assumed to be independent and identically distributed.
"""

import sys
import argparse
from math import sqrt

# Start of statistical function definitions
def mean(list):
    return sum(list) / len(list)

def variance(list):
    squaredList = [x*x for x in list]
    return sum(squaredList)/len(list) - (sum(list)/len(list))**2

def covariance(list1, list2):	# Don't think I actually need it for tihs code, but leaving it here anyway
    if len(list1) != len(list2):
        return 0
    else:
        return sum([x*y for x,y in zip(list1, list2)]) / len(list1) -  sum(list1)*sum(list2)/(len(list1)**2)
# End of statistical function definitions

parser = argparse.ArgumentParser(description = 'Calculate means and errors (standard deviations of the mean) of data sets. This program is intended for use on the output files of Rxy_shell_spectrum.py, for runs with different random seeds but identical parameters.')
parser.add_argument('inputFilenames', metavar = 'input file', type = str, nargs = '+', help = 'Input files. Two is the minimum number for which it makes any sense to even run this progam!')
parser.add_argument('--output', metavar = 'output file', type = str, default = 'mean_error_output.dat', help = 'Output file. Defaults to mean_error_output.dat')
parser.add_argument('--verbose', action = 'store_true', help = 'Verbose mode.')

args = parser.parse_args()
verbose = args.verbose

# Open input files
inputFiles = []
for inputFilename in args.inputFilenames:
	try:
		inputFiles.append(open(inputFilename, 'r'))
	except IOError:
		print 'Unable to open ' + inputFilename + ' for reading. Exitting...'
		for inputFile in inputFiles:
			inputFile.close()
		sys.exit(2)

# Open output file
try:
	outputFile = open(args.output, 'w')
except IOError:
	print 'Unable to open ' + args.output +' for writing. Exitting...'
	for inputFile in inputFiles:
		inputFile.close()
	sys.exit(2)

#Read in the top line of each input file
headerLine = []
headerLine = inputFiles[0].readline().split()
for file in inputFiles[1:]:	# Throw away the header line of each file
	if headerLine != file.readline().split() and verbose:
		# Warn the user if the headers don't seem to match
		print 'Header line of ' + file.name + ' does not seem to match that of ' + inputFiles[0].name

# Create new header line: assume the first entry is dependent variable, with no mean or error to be calculated
if headerLine[0] == '#':
	headerLine[0] = headerLine[0] + headerLine[1]
	headerLine.__delitem__(1)

newHeaderLine = [headerLine[0]]

for item in headerLine[1:]:
	newHeaderLine.append(item + '_mean')
	newHeaderLine.append(item + '_error')

for item in newHeaderLine:
	outputFile.write(item + '\t')
outputFile.write('\n')

for line in inputFiles[0]:
	lines = [line.split()]
	t = float(lines[0][0])	# Independent variable
	for file in inputFiles[1:]:
		lines.append(file.readline().split())
		if verbose and float(lines[-1][0]) != t:
			print 'Warning, dependent variable values do not match!'
		if len(lines[-1]) != len(lines[0]):
			print 'Number of dependent variables for not match, exitting...'
			sys.exit(2)
	outputFile.write(str(t) + '\t')
	depVarTranspose = [[float(lines[j][i]) for j in range(0,len(lines))] for i in range(1,len(lines[0]))]
	for var in depVarTranspose:
		outputFile.write(str(mean(var)) + '\t')
		outputFile.write(str(sqrt(variance(var)/len(var))) + '\t')
	outputFile.write('\n')
	

# Clean up - close the output and input files
outputFile.close()

for inputFile in inputFiles:
	inputFile.close()
	sys.exit(2)
