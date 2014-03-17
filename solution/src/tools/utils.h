/*
 * Do not remove.
 * MPI/OpenMP training courses
 * Adrien Cassagne, ASA - CINES, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef UTILS_H_
#define UTILS_H_

#include "stdlib.h"
#include "stdio.h"
#include "string.h"

#define TIME_DIFF(t1, t2) \
		((t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_usec - t1.tv_usec))

void printIterationTimeAndPerformance(unsigned long iIte, struct timeval t1, struct timeval  t2, struct timeval t3)
{
	if(Verbose && !RankMPI)
		printf("Processing step %ld takes %f ms (total time elapsed %f ms)\n", iIte,
		                                                                       TIME_DIFF(t2, t3)/1000.0,
		                                                                       TIME_DIFF(t1, t3)/1000.0);
}

void writeOuputFile(unsigned long iIte, const plan *p)
{
	if(WriteToFiles)
	{
		char outFileName[1024];
		sprintf(outFileName, "data/out/out.%ld.dat", iIte);

		for(int iRank = 0; iRank < SizeMPI; ++iRank)
		{
			if(iRank == RankMPI)
				writePlanIntoFile(p, outFileName);
#ifdef OPEN_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif
		}
	}
}

void argumentReader(int argc, char** argv)
{
	if(argc < 5)
	{
		printf("usage: %s -f fileName -i nIterations [-v] [-w]\n", argv[0]);
		printf("usage: %s -n nBody    -i nIterations [-v] [-w]\n", argv[0]);
		printf("       -f           root name of body file to read.\n");
		printf("       -n           approximate number of bodies for random generation.\n");
		printf("       -i           number of iterations to compute.\n");
		printf("       -v           activate verbose mode.\n");
		printf("       -w           write results into 'data/out/out.*.dat' files.\n");

		exit(0);
	}

	if(argc >= 5)
	{
		if(argv[1][0] == '-' && argv[1][1] == 'f')
			strcpy (FileName, argv[2]);
		else if (argv[1][0] == '-' && argv[1][1] == 'n')
			NBody = atoi(argv[2]);
		else
		{
			printf("Unexpected parameter \"%s\". Program exiting.\n", argv[1]);
			exit(0);
		}

		if(argv[3][0] == '-' && argv[3][1] == 'i')
			NIterations = atoi(argv[4]);
		else
		{
			printf("Unexpected parameter \"%s\". Program exiting.\n", argv[3]);
			exit(0);
		}

		if(argc >= 6)
			Verbose = (argv[5][0] == '-' && argv[5][1] == 'v') ? 1 : 0;
		else
			WriteToFiles = 0;

		if(argc >= 7)
			WriteToFiles = (argv[6][0] == '-' && argv[6][1] == 'w') ? 1 : 0;
		else
			WriteToFiles = 0;
	}
}

#endif /* UTILS_H_ */
