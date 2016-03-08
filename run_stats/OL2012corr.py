#! /usr/bin/python

# This script is intended to estimate the auto-correlation coefficients
# gamma_k for the time-series of OL2012A and OL2012B
# This is under the assumption that they are stationary processes
# and that the outputs are equispaced in time (which is very nearly true)

from argparse import ArgumentParser

timevarFilename = '../timevar'
outfileName = './autocorrelation.dat'
tStartDefault = 50.0
conf95 = 1.959964   # The 97.5% percentile of the standard normal distribution

parser = ArgumentParser(description= "Compute and output the auto-correlation function for the OL2012A and OL2012B quantities in the timevar file in the parent directory")
parser.add_argument('--start', type=float, help='Time to begin taking statistics from')
parser.add_argument('--indepStep', type=int, help='Number of steps apart before we regard points as independent')

args = parser.parse_args()

if args.start is None:
    tStart = tStartDefault
else:
    tStart = args.start

indepStep = args.indepStep

timevarFile = open(timevarFilename, "r")
header = timevarFile.readline().split()

tIndex = header.index('t')
AIndex = header.index('OL2012A')

# Read in lists
tList = []
AList = []
for line in timevarFile:
    line = line.split()
    tList.append(float(line[tIndex]))
    AList.append(float(line[AIndex]))

timevarFile.close()

# Reducde lists to only times after tStart
AList = [a for a,t in zip(AList,tList) if t >= tStart]
tList = [t for t in tList if t >= tStart]
numPoints = len(tList)

# Compute the means
tMean = sum(tList) / numPoints
AMean = sum(AList) / numPoints

if indepStep is None:
    # Compute autocovariance function
    gammaAList = []
    for k in range (0, numPoints-1):
        gamma_k = 0
        for t in range(0, numPoints-k):
            gamma_k += (AList[t] - AMean)*(AList[t+k] - AMean)
        gamma_k /= (numPoints- k - 1)   # (numPoints-k) points summed over, and we want the unbiased estimator
        gammaAList.append(gamma_k)

    # Compute the autocorrelation function
    rhoAList = []
    for k in range(0, numPoints-1):
        rhoAList.append(gammaAList[k]/gammaAList[0])

    # Write the output to a file
    outfile = open(outfileName, "w")

    outfile.write('#k\tgamma_A_k\trho_A_k\n')

    for k in range(0, numPoints-1):
        outfile.write(str(k) + '\t' + str(gammaAList[k]) + '\t' + str(rhoAList[k]) + '\n')

    outfile.close()
else:
    # Strip out points to build a list that is almost statistically independent
    tSmallList = [tList[k] for k in range(numPoints) if k%indepStep == 0]
    ASmallList = [AList[k] for k in range(numPoints) if k%indepStep == 0]
    numSmallPoints = len(tSmallList)
    print 'Using {} points'.format(numSmallPoints)

    # Compute the sample variance of the ASmalList
    ASmallMean = sum(ASmallList)/numSmallPoints
    ASmallVar = sum([(A - ASmallMean)**2 for A in ASmallList])/(numSmallPoints - 1)
    # Estimate the standard deviation of the mean
    AMeanSD = (ASmallVar / numSmallPoints)**0.5
    # Compute the 95% confidence interval
    AConfMax = ASmallMean + AMeanSD * conf95
    AConfMin = ASmallMean - AMeanSD * conf95
    # Output this
    print "Mean of OL2012A is {}".format(ASmallMean)
    print "95% confidence interval for OL2012A is [{}, {}]".format(AConfMin, AConfMax)
