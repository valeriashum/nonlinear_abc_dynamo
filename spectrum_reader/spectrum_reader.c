#include <stdio.h>
#include <stdlib.h>
#include "./tsvlib/tsvlib.h"

#define VERBOSE

int i;
float t_init, lastSpectrumTime;
long int numLines, numSpectra, numTimes;

char *specNames[] = {"vx*vx","vy*vy","vz*vz","bx*bx/tvx*tvx","by*by/tvy*tvy","bz*bz/tvz*tvz","th*th","vx*vy","bx*by/(2 tide spectra!)","?","?","?","0","0","0","0","0","0","0","0","0"};

void printMenu();	// Print the main menu
void findAndWriteSpectrum(FILE *spectrumInFile, long int numSpectra);
void writeMeanSpectrum(FILE *spectrumInFile, long int numSpectra, float lastSpectrumTime);

int main(int argc, char *argv[]){

	char *spectrumInFilename, menuChar;
	FILE *spectrumInFile;
	struct dataLine *spectrumLine = malloc(sizeof(struct dataLine));
	spectrumLine->charArray = NULL; spectrumLine->recordPos = NULL;
	struct floatRecords *spectrumFloats = malloc(sizeof(struct floatRecords));
	
	if(argc != 2){
		printf("Incorrect number of arguments! Correct usage is:\nspectrum_reader\t[spectrum datafile]\n");
		return 0;
	} else {
		spectrumInFilename = argv[1];
	}

	spectrumInFile = fopen(spectrumInFilename, "r");
	if(spectrumInFile == NULL){
		printf("Failed to open %s for reading. Quitting.\n",spectrumInFilename);
		return 0;
	}

	numLines = getNumLines(spectrumInFile);
	fseek(spectrumInFile, 0, SEEK_SET);
	for(i=1; i<numLines; i++){
		skipLine(spectrumInFile);
	}
	gobbleLine(spectrumInFile, spectrumLine);
	dataLineToFloats(spectrumLine, spectrumFloats);
	lastSpectrumTime = spectrumFloats->floatArray[0];
#ifdef VERBOSE
	printf("%s contains %ld lines.\n", spectrumInFilename, numLines);
	printf("The last spectrum in the file is from time %f\n", lastSpectrumTime);
#endif
	// Return to the beginning of the file, and chuck away the first two lines
	// (Line 1 is just a list of |k|, and line 2 is the number of modes in that bin.)
	fseek(spectrumInFile, 0, SEEK_SET);
	gobbleLine(spectrumInFile, spectrumLine);
	gobbleLine(spectrumInFile, spectrumLine);
	// Look at what the first time in the file is, the find out how many lines begin with this
	gobbleLine(spectrumInFile, spectrumLine);
	dataLineToFloats(spectrumLine, spectrumFloats);
	t_init = spectrumFloats->floatArray[0];
	for(i=0; t_init == spectrumFloats->floatArray[0]; i++){
		gobbleLine(spectrumInFile, spectrumLine);
		dataLineToFloats(spectrumLine, spectrumFloats);
	}
	numSpectra = i;
#ifdef VERBOSE
	printf("%s contains %ld different type of spectra (based on counting how many lines begin with the time %f).\n", spectrumInFilename, numSpectra, t_init);
	if(numSpectra == 21) printf("Good, since 21 spectra is what I expected to find!");
	printf("The next line begins with the time %f.\n",spectrumFloats->floatArray[0]);
#endif

	int menuQuit=0;
	while(menuQuit==0){
		printMenu();
		scanf("%c",&menuChar);
		if (menuChar=='1'){
			findAndWriteSpectrum(spectrumInFile, numSpectra);
		}else if(menuChar=='2')	writeMeanSpectrum(spectrumInFile, numSpectra, lastSpectrumTime);
		else if(menuChar=='3') menuQuit=1;
	}


	return 0;
}

void printMenu(){
	// FIXME - add code to give some information, eg the times the file holds
	fprintf(stdout,"%s\n%s\n%s\n%s\n%s",
	"What do you want to do?",
	"1) Output spectrum from a specified time.",
	"2) Output a time-averaged spectrum.",
	"3) Quit.",
	" > ");
}

void findAndWriteSpectrum(FILE *spectrumInFile, long int numSpectra){
	int i;
	float outputTime;
	char *outFileName = malloc(30*sizeof(char));
	FILE *outFile;
	int spectrumNumber;
	struct dataLine *spectrumLine = malloc(sizeof(struct dataLine));
	spectrumLine->charArray = NULL; spectrumLine->recordPos = NULL;
	struct floatRecords *spectrumFloats = malloc(sizeof(struct floatRecords));
	spectrumFloats->floatArray = NULL;
	struct floatRecords *kFloats = malloc(sizeof(struct floatRecords));
	kFloats->floatArray = NULL;
	struct floatRecords *binSizeFloats = malloc(sizeof(struct floatRecords));
	binSizeFloats->floatArray = NULL;

	fprintf(stdout, "%s\n%s",
	"Please specify the program time from which you want a spectrum:",
	" > ");
	scanf("%f", &outputTime);
	i=1;
	while(i){
		fprintf(stdout, "%s\n",
		"Please specify the spectrum number you want to output:");
		if(numSpectra == 21){
			fprintf(stdout,"%s\n","(If this file follows the standard SNOOPY output rules for spectra, they should be as follows: ");
			for(i=0; i<21; i++) fprintf(stdout,"%d: %s,  ", i+1, specNames[i]);
			fprintf(stdout, ")\n");
		}
		fprintf(stdout," > ");
		scanf("%d", &spectrumNumber);
		if(1 <= spectrumNumber && spectrumNumber <= numSpectra) i=0;
	}
#if 0
	// FIXME - generate the default filename
	defaultName[0] = 'S'; defaultName[1] = '\0';
	int keepGoing = 1;
	while(keepGoing){
		fprintf(stdout, "%s\n%s%s%s\n%s",
		"Please specify the filename to which the spectrum should be written",
		" (default is ", defaultName, ")",
		" > ");
		scanf("%s", outFileName);
		// Check that the filename is valid
		if(outFileName[0] == '\0'){
			outFile = fopen(defaultName, "w");
			if(outFile!=NULL) keepGoing = 0;
			else printf("Filename %s could not be opened for reading, please try again\n > ",defaultName);
		} else {
			outFile = fopen(outFileName, "w");
			if(outFile!=NULL) keepGoing = 0;
			else printf("Filename %s could not be opened for reading, please try again\n > ",outFileName);
		}
	}
#else
	outFileName = "theFileICouldntThinkOfANameFor.dat";
	outFile = fopen(outFileName, "w");
	printf("Outputting to %s\n", outFileName);
#endif

	// Move to the start of the file and then find the next time after the time specified
	fseek(spectrumInFile, 0, SEEK_SET);
	// Read in the first two lines
	gobbleLine(spectrumInFile, spectrumLine);
	dataLineToFloats(spectrumLine, kFloats);
	gobbleLine(spectrumInFile, spectrumLine);
	dataLineToFloats(spectrumLine, binSizeFloats);

	float currentSpectrumTime = -1.0;
	while(currentSpectrumTime < outputTime){
		// Update the time based on the first spectrum
		gobbleLine(spectrumInFile, spectrumLine);
		dataLineToFloats(spectrumLine, spectrumFloats);
		currentSpectrumTime = spectrumFloats->floatArray[0];
	}
	// Go to the correct spectrum
	for(i=1; i<spectrumNumber; i++){
		gobbleLine(spectrumInFile, spectrumLine);
	}
	dataLineToFloats(spectrumLine, spectrumFloats);

	fprintf(stdout, "%s%d%s%f%s%s\n",
	"Wrtiting spectrum number ", spectrumNumber,
	" at time ", currentSpectrumTime,
	" to the file ", outFileName);

	for(i=0; i<kFloats->numRecords; i++){
		fprintf(outFile,"%lf\t%e\n",
		kFloats->floatArray[i],
		spectrumFloats->floatArray[i+1]);
	}

	printf("\nSpectrum written.\n\n");
	fclose(outFile);

}

void writeMeanSpectrum(FILE* spectrumInFile, long int numSpectra, float lastSpectrumTime){
	int i, j;	// Dummy counter variables
	char *outFileName = malloc(30*sizeof(char));
	FILE *outFile;
	float startTime, endTime;
	float *meanSpectrum;
	int numAvOver;	// The number of spectra to average over
	int spectrumNumber; 	// Whether we want eg. vx*vx or bz*bz etc.
	struct dataLine *spectrumLine = malloc(sizeof(struct dataLine));
	spectrumLine->charArray = NULL; spectrumLine->recordPos = NULL;
	struct floatRecords *spectrumFloats = malloc(sizeof(struct floatRecords));
	spectrumFloats->floatArray = NULL;
	struct floatRecords *kFloats = malloc(sizeof(struct floatRecords));
	kFloats->floatArray = NULL;

	/* Start of user supplying information */
	printf("What time do you want the averaging to start from?\n > ");
	scanf("%f", &startTime);
	startTime = startTime > 0 ? startTime : 0;
	printf("What time do you want to average until? (The latest time in the file is %f)\n > ", lastSpectrumTime);
	scanf("%f", &endTime);
	endTime = (endTime > lastSpectrumTime ? lastSpectrumTime : endTime);
	printf("How many spectra do you want to average over? (Will try and get a uniform spread of times)\n > ");
	scanf("%d", &numAvOver);
	numAvOver = numAvOver > 1 ? numAvOver : 1;
	i=1;
	while(i){
		fprintf(stdout, "%s\n",
		"Please specify the spectrum number you want to output:");
		if(numSpectra == 21){
			fprintf(stdout,"%s\n","(If this file follows the standard SNOOPY output rules for spectra, they should be as follows: ");
			for(i=0; i<21; i++) fprintf(stdout,"%d: %s,  ", i+1, specNames[i]);
			fprintf(stdout, ")\n");
		}
		fprintf(stdout," > ");
		scanf("%d", &spectrumNumber);
		if(1 <= spectrumNumber && spectrumNumber <= numSpectra) i=0;
	}
	/* End of user supplying information */

	outFileName = "theFileICouldntThinkOfANameFor.dat";
	outFile = fopen(outFileName, "w");
	printf("Outputting to %s\n", outFileName);

	// Read in the k-value for each bin
	fseek(spectrumInFile, 0, SEEK_SET);
	gobbleLine(spectrumInFile, spectrumLine);
	dataLineToFloats(spectrumLine, kFloats);
	skipLine(spectrumInFile);	// Skip over the bin sizes

	// Allocate memory to hold the time-averaged spectrum
	meanSpectrum = malloc(sizeof(float)*kFloats->numRecords);
	for(j=0; j<kFloats->numRecords; j++) meanSpectrum[j] = 0.0;
	
	float currentSpectrumTime = -1.0;	// Make sure it's less than the first spectrum in the file
	float targetSpectrumTime = startTime;
	printf("Using spectra from program times");
	for(i=0; i<numAvOver; i++){	// Loop over the various times we want to average over
		while(currentSpectrumTime < targetSpectrumTime){	// Loop through until we find a time not less than the one we want
			// Update the time based on the first spectrum
			gobbleLine(spectrumInFile, spectrumLine);
			dataLineToFloats(spectrumLine, spectrumFloats);
			currentSpectrumTime = spectrumFloats->floatArray[0];
		}
		printf("\t%f", currentSpectrumTime);
		// Update the target time for the benefit of the next loop
		if(numAvOver>1) targetSpectrumTime += (endTime - startTime)/(numAvOver - 1);

		// Go to the correct spectrum
		for(j=1; j<spectrumNumber; j++){
			gobbleLine(spectrumInFile, spectrumLine);
		}
		dataLineToFloats(spectrumLine, spectrumFloats);

		// Add this spectrum to our average array
		for(j=0; j<kFloats->numRecords; j++) meanSpectrum[j] += spectrumFloats->floatArray[j+1];
	}
	printf("\n");
	// Turn the sum into the mean
	for(j=0; j<kFloats->numRecords; j++) meanSpectrum[j] /= numAvOver;

	// Write to the file
	for(j=0; j<kFloats->numRecords; j++){
		fprintf(outFile,"%lf\t%lf\n",
		kFloats->floatArray[j],
		meanSpectrum[j]);
	}

	// Clean-up
	fclose(outFile);
	free(meanSpectrum); free(spectrumLine); free(spectrumFloats); free(kFloats);
}
