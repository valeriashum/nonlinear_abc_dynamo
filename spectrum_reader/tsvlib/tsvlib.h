// My library of functions for dealing with tab separated value files
// FIXME - write something to do all the memory freeing for the structs automatically

#ifndef _TSVLIB_H_
#define _TSVLIB_H_

#include <stdio.h>
#include <stdlib.h>

struct dataLine{
	int numRecords;
	char *charArray;
	long int *recordPos;
};

struct floatRecords{
	int numRecords;
	float *floatArray;
};

#define TSV_COMMENT_CHAR '#'

long int getNumLines(FILE *stream);	// Get the number of (non-empty) lines in the file

int isComment(struct dataLine *lineIn);	// Determine whether or not the line is a comment FIXME - write this

/* The following functions return and int to indicate whether or not the function call was successful */
//	0	Success
//	1	EOF encountered
//	-1	stream pointer was null
int gobbleLine(FILE *stream, struct dataLine *lineOut);	// Processes a line of tsv data into lineOut, advancing the stream position.

int skipLine(FILE *stream);	// Advances the file position to the next line

int dataLineToFloats(struct dataLine *lineIn, struct floatRecords *floatsOut);	// Processes a dataLine object into an array of floats

#endif
