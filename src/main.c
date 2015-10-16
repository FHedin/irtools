/*
 * Copyright (c) 2013-2015, Florent Hédin, Pierre-André Cazade, Markus Meuwly, 
 * and the University of Basel, Switzerland. All rights reserved.
 * 
 * The 3-clause BSD license is applied to this software.
 * 
 * See LICENSE.txt
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <fftw3.h>

#include "globals.h"

#include "ir2d.h"

#ifdef USE_CUDA
#include "ir2dcuda.h"
#endif

// prints usage information
void help(FILE* STREAM, char* argv[])
{
#ifdef USE_CUDA
    fprintf(STREAM,"For a standard cpu use :\n");
    fprintf(STREAM,"Usage : %s -i inputFile -o outputFile\n\n",argv[0]);
    fprintf(STREAM,"For using a cuda gpu, add -cuda (or --cuda or -gpu or --gpu) :\n");
    fprintf(STREAM,"Usage : %s -i inputFile -o outputFile -cuda\n",argv[0]);
#else
    fprintf(STREAM,"Usage : %s -i inputFile -o outputFile\n",argv[0]);
#endif
}

int main(int argc, char* argv[])
{

    if (argc < 5)
    {
        help(stderr,argv);
        return EXIT_FAILURE;
    }
    char *inpName,*outName;

    bool useCuda = false;

    // arguments parsing
    for (int i=1; i<argc; i++)
    {
        // get name of input file
        if (!strcmp(argv[i],"-i"))
        {
            inpName = argv[++i];
        }
        // get name of output file
        else if (!strcmp(argv[i],"-o"))
        {
            outName = argv[++i];
        }
        // print help and proper exit
        else if ( !strcmp(argv[i],"-h") || !strcmp(argv[i],"--h")|| !strcmp(argv[i],"-help") || !strcmp(argv[i],"--help") )
        {
            help(stdout,argv);
            return EXIT_SUCCESS;
        }
#ifdef USE_CUDA
        // if CUDA enabled
        else if ( !strcmp(argv[i],"-cuda") || !strcmp(argv[i],"--cuda") || !strcmp(argv[i],"-gpu") || !strcmp(argv[i],"--gpu"))
        {
            useCuda = true;
        }
#endif
        // error if unknown command line option
        else
        {
            fprintf(stderr,"[Error] Argument '%s' is unknown.\n",argv[i]);
            help(stderr,argv);
            return EXIT_FAILURE;
        }
    }

    // open input and output file, with error check
    FILE* input=NULL;
    input = fopen(inpName,"rt");
    if(input==NULL)
    {
        fprintf(stderr,"Error while opening input file '%s' in text read mode. Check if file exists or if the path is not wrong ?\n",inpName);
        return EXIT_FAILURE;
    }

    FILE* output=NULL;
    output = fopen(outName,"wt");
    if(input==NULL)
    {
        fprintf(stderr,"Error while opening output file '%s' in text write mode. Check if file exists or if the path is not wrong ?\n",outName);
        return EXIT_FAILURE;
    }

    int status = EXIT_SUCCESS;
#ifdef USE_CUDA
    if (useCuda)
        status = ir2d_cuda(input,output);
    else
#endif
        status = ir2d(input,output);

    if(status != EXIT_SUCCESS)
    {

        fprintf(stderr,"Something wrong with calculation !!!! \n");
        return EXIT_FAILURE;

    }

    // if everything OK close files and exit properly
    fclose(input);
    fclose(output);

    return EXIT_SUCCESS;
}

