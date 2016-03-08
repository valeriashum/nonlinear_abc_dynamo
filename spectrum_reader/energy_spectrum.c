// Simplified version of my spectrum_reader program.
// This should take command line arguments and output a kinetic energy spectrum
// Correct usage is:
// ./energy_spectrum [spectrum_data_file] [start_time] [end_time] [num_average] [outfile_name]
// to average over [num_average] spectra roughly evenly spaced between [start_time] and [end_time]
// ./energy_spectrum [spectrum_data_file] [time] [outfile_name]
// to produce a spectrum for the ealiest time after [time] for which a spectrum is available

#include <stdio.h>
#include <stdlib.h>

//#define VERBOSE
#define MAX_LINE_LENGTH 1024

int main(int argc, char *argv[]){
	char *spectrumInFilename;
	FILE *spectrumInFile, *spectrumOutFile;
	char line[MAX_LINE_LENGTH];	// Holds lines read in from file

	unsigned short nbins;		// Number of wavenumber bins
	float *bink, *binMembers;	// Wavenumber value assigned to each bin and number of wavevectors in each bin
	double *binEnergy, *totalBinEnergy;			// Energy held in each bin in vx^2, vy^2 and vz^2, and the total energy
	unsigned int i, j, k, l;	// Dummy counter variables
	float time, min_time, max_time, next_time;				// Holds the time read from each line, and the minimum and maximum times specified on the command line, and the next time at which we should take a spectrum
	unsigned int numAverage;	// Number of spectra to average over

	if(argc != 6 && argc != 4){
		fprintf(stderr, "Incorrect number of arguments! Correct usage is:\n./energy_spectrum [spectrum_data_file] [start_time] [end_time] [num_average] [outfile_name]\nor:\n./energy_spectrum [spectrum_data_file] [time] [outfile_name]\n");
		return -1;
	} else {
		spectrumInFilename = argv[1];
	}

	spectrumInFile = fopen(spectrumInFilename, "r");
	if(spectrumInFile == NULL){
		fprintf(stderr, "Failed to open %s for reading. Quitting.\n", spectrumInFilename);
		return -1;
	}
	
	// Count how many k bins the file contains
	nbins=0;
	fgets(line, MAX_LINE_LENGTH, spectrumInFile);
	for(i=0; line[i]!='\0'; i++){
		if(line[i] == '\t') nbins++;
	}

#ifdef VERBOSE
	fprintf(stderr, "Found %hu wavenumber bins.\n", nbins);
#endif


	// Read what the nominal values of the kbins actually are
	bink = malloc(nbins*sizeof(float));
	j = 0;
	for(i=0; i<nbins; i++){
		sscanf(line + j*sizeof(char), "%f\t%n", &bink[i], &k);
		j +=k;
	}

#ifdef VERBOSE
	fprintf(stderr, "Found the following nominal k values for the bins:\n");
	for(i=0; i<nbins; i++)
		fprintf(stderr, "%f\n", bink[i]);
#endif

	// Read how many members are in each bin
	fgets(line, MAX_LINE_LENGTH, spectrumInFile);
	binMembers = malloc(nbins*sizeof(float));
	j = 0;
	for(i=0; i<nbins; i++){
		sscanf(line + j*sizeof(char), "%f\t%n", &binMembers[i], &k);
		j +=k;
	}

#ifdef VERBOSE
	fprintf(stderr, "Found the following number of members for the bins:\n");
	for(i=0; i<nbins; i++)
		fprintf(stderr, "%f\n", binMembers[i]);
#endif

	binEnergy = malloc(nbins*sizeof(double));
	totalBinEnergy = malloc(nbins*sizeof(double));

	if(argc == 4){
		// We've been asked to output the earliest spectrum at or after [time] given by argv[2]
		min_time = atof(argv[2]);
		do{
			if(fgets(line, MAX_LINE_LENGTH, spectrumInFile) == NULL){
				fprintf(stderr, "Reached enf of file before finding a spectrum from a time at or later than %f.\n", min_time);
				return -1;
			}
			sscanf(line, "%f\t%n", &time, &j);
		}while(time < min_time);

		// line should now hold the first spectrum for the requested time, which is the vx*vx spectrum
		spectrumOutFile = fopen(argv[3], "w");
		if(spectrumOutFile == NULL){
			fprintf(stderr, "Failed to open %s for writing. Quitting.\n", argv[3]);
			return -1;
		}
		// Write a header
		fprintf(spectrumOutFile, "#k\tE\n");

		for(i=0; i<nbins; i++){
			sscanf(line + j*sizeof(char), "%lf\t%n", &totalBinEnergy[i], &k);
			j += k;
		}

		// Now read the vy*vy spectrum and add it to the total
		fgets(line, MAX_LINE_LENGTH, spectrumInFile);
		sscanf(line, "%f\t%n", &time, &j);	// Need to be here to set j correctly
		for(i=0; i<nbins; i++){
			sscanf(line + j*sizeof(char), "%lf\t%n", &binEnergy[i], &k);
			j += k;
			totalBinEnergy[i] += binEnergy[i];
		}

		// Now read the vz*vz spectrum and add it to the total
		fgets(line, MAX_LINE_LENGTH, spectrumInFile);
		sscanf(line, "%f\t%n", &time, &j);	// Need to be here to set j correctly
		for(i=0; i<nbins; i++){
			sscanf(line + j*sizeof(char), "%lf\t%n", &binEnergy[i], &k);
			j += k;
			totalBinEnergy[i] += binEnergy[i];
			// Energy has a factor of 0.5
			totalBinEnergy[i] *= 0.5;
			// Write to the file
			fprintf(spectrumOutFile, "%f\t%08e\n", bink[i], totalBinEnergy[i]);
		}
		
		// Close the output file
		fclose(spectrumOutFile);

	}

	if(argc == 6){
		// We've been asked to output the earliest spectrum at or after [time] given by argv[2]
		min_time = atof(argv[2]);
		max_time = atof(argv[3]);
		numAverage = atoi(argv[4]);
		
		// Open the output file
		spectrumOutFile = fopen(argv[5], "w");
		if(spectrumOutFile == NULL){
			fprintf(stderr, "Failed to open %s for writing. Quitting.\n", argv[5]);
			return -1;
		}
		// Write a header
		fprintf(spectrumOutFile, "#k\tE\n");
		
		for(i=0; i<nbins; i++)
			totalBinEnergy[i] = 0.0;

		next_time = min_time;
		for(l=0; l< numAverage; l++){
			do{
				if(fgets(line, MAX_LINE_LENGTH, spectrumInFile) == NULL){
					fprintf(stderr, "Reached enf of file before finding a spectrum from a time at or later than %f.\n", next_time);
					return -1;
				}
				sscanf(line, "%f\t%n", &time, &j);
			}while(time < next_time);
	
			// line should now hold the first spectrum for the requested time, which is the vx*vx spectrum
			for(i=0; i<nbins; i++){
				sscanf(line + j*sizeof(char), "%lf\t%n", &binEnergy[i], &k);
				j += k;
				totalBinEnergy[i] += binEnergy[i];
			}
	
			// Now read the vy*vy spectrum and add it to the total
			fgets(line, MAX_LINE_LENGTH, spectrumInFile);
			sscanf(line, "%f\t%n", &time, &j);	// Need to be here to set j correctly
			for(i=0; i<nbins; i++){
				sscanf(line + j*sizeof(char), "%lf\t%n", &binEnergy[i], &k);
				j += k;
				totalBinEnergy[i] += binEnergy[i];
			}
	
			// Now read the vz*vz spectrum and add it to the total
			fgets(line, MAX_LINE_LENGTH, spectrumInFile);
			sscanf(line, "%f\t%n", &time, &j);	// Need to be here to set j correctly
			for(i=0; i<nbins; i++){
				sscanf(line + j*sizeof(char), "%lf\t%n", &binEnergy[i], &k);
				j += k;
				totalBinEnergy[i] += binEnergy[i];
			}
			
			next_time += (max_time - min_time) / numAverage;
		}
		
		for(i=0; i<nbins; i++){
			// Energy has a factor of 0.5
			totalBinEnergy[i] *= 0.5;
			// It's an average
			totalBinEnergy[i] /= numAverage;
			// Write to the file
			fprintf(spectrumOutFile, "%f\t%08e\n", bink[i], totalBinEnergy[i]);
		}

		// Close the output file
		fclose(spectrumOutFile);

	}
	

	return 0;
}
