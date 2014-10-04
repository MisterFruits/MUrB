/*
 * Do not remove.
 * MPI/OpenMP training courses
 * Adrien Cassagne, ASA - CINES, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#include "math.h"
#include "stdio.h"
#include "unistd.h"
#include "stdlib.h"
#include "assert.h"
#include "sys/time.h"

#include "body.h"
#include "plan.h"

// global
char           FileName[2048];
unsigned long  NBody;
unsigned long  NIterations;
unsigned short WriteToFiles = 1;
unsigned short Verbose      = 0;

// global MPI
int SizeMPI = 1;
int RankMPI = 0;

#include "tools/utils.h"

int main(int argc, char** argv)
{
	struct timeval t1, t2, t3;

	// read arguments on command line
	argumentReader(argc, argv);

	// create plan in reading file or generate bodies randomly
	plan *p;
	if(FileName[0] == '\0')
	{
		p = createPlan(NBody);
		fillRandom(p);
	}
	else
	{
		p = createPlanWithFile(FileName);
		NBody = p->nBody;
	}

	// display simulation config
	if(!RankMPI)
	{
		printf("N body started !\n");
		if(FileName[0] != '\0')
			printf("  -> fileName    : %s.*.dat\n", FileName);
		printf("  -> nBody       : %ld\n", NBody);
		printf("  -> nIterations : %ld\n", NIterations);
		printf("  -> verbose     : %d\n", Verbose);
		printf("  -> writeToFiles: %d\n\n", WriteToFiles);
	}

	double dt;
	writeOuputFile(0, p);
	gettimeofday(&t1,NULL);

	if(!RankMPI)
		printf("Starting simulation...\n");
	for(unsigned long iIte = 1; iIte <= NIterations; ++iIte) {
		/*******************************/
		/*** Simulation computations ***/
		gettimeofday(&t2,NULL);

		computeAllLocalAcceleration(p);
		dt = findLocalDt(p);
		updateAllLocalPositionAndSpeed(p, dt);

		gettimeofday(&t3,NULL);
		/*** Simulation computations ***/
		/*******************************/

		printIterationTimeAndPerformance(iIte, t1, t2, t3);
		writeOuputFile(iIte, p);
	}
	if(!RankMPI)
		printf("Simulation end... exiting.\n");

	destroyPlan(p);

	return EXIT_SUCCESS;
}
