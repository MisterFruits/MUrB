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
#include "mpi/mpi.h"
#include "sys/time.h"

#include "body.h"
#include "plan.h"

// global
char           FileName[2048];
unsigned long  NBody;
unsigned long  NBodyPerNode;
unsigned long  NIterations;
unsigned short WriteToFiles = 1;
unsigned short Verbose      = 0;

// global MPI
int SizeMPI;
int RankMPI;

#include "tools/utils.h"

int main(int argc, char** argv)
{
	// MPI init
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &SizeMPI);
	MPI_Comm_rank(MPI_COMM_WORLD, &RankMPI);

	assert(SizeMPI == 1 || SizeMPI % 2 == 0);

	struct timeval t1, t2, t3;

	// read arguments on command line
	argumentReader(argc, argv);

	// create plan in reading file or generate bodies randomly
	plan *p;
	if(FileName[0] == '\0')
	{
		NBodyPerNode = NBody / SizeMPI;
		NBody = SizeMPI * NBodyPerNode;
		p = createPlan(NBodyPerNode);
		fillRandom(p);
	}
	else
	{
		p = createPlanWithFile(FileName);
		NBodyPerNode = p->nBody;
		NBody = p->nBody * SizeMPI;
	}

	// display simulation config
	if(!RankMPI)
	{
		printf("N body started !\n");
		if(FileName[0] != '\0')
			printf("  -> fileName      : %s.*.dat\n", FileName);
		printf("  -> sizeMPI       : %d\n", SizeMPI);
		printf("  -> rankMPI       : %d\n", RankMPI);
		printf("  -> nBody         : %ld\n", NBody);
		printf("  -> NBodyPerNode  : %ld\n", NBodyPerNode);
		printf("  -> nIterations   : %ld\n", NIterations);
		printf("  -> verbose       : %d\n", Verbose);
		printf("  -> writeToFiles  : %d\n", WriteToFiles);
		printf("  -> planSize (Ko) : %lf\n\n", p->nBody * (8 * sizeof(double)) / 1024.0);
	}

	// create MPI data structure
	MPI_Datatype bodyMPI;
	body b;
	int blockLengths[3] = {1, 1, 1};
	MPI_Aint displacements[3];
	displacements[0] = (char*) &b.mass - (char*) &b; // char* because sizeof(char) = 1
	displacements[1] = (char*) &b.posX - (char*) &b; // char* because sizeof(char) = 1
	displacements[2] = (char*) &b.posY - (char*) &b; // char* because sizeof(char) = 1
	MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
	MPI_Type_create_struct(3, blockLengths, displacements, types, &bodyMPI);

	// commit MPI data structure
	MPI_Type_commit(&bodyMPI);

	// who is next or previous MPI process ?
	int prevRankMPI = (RankMPI == 0) ? SizeMPI -1 : RankMPI -1;
	int nextRankMPI = (RankMPI +1) % SizeMPI;

	// allocate two MPI buffers which can contains all bodies of plan p
	body **bodyBufferMPI = (body**)malloc(2 * sizeof(body*));
	bodyBufferMPI[0]     = (body*) malloc(p->nBody * sizeof(body));
	bodyBufferMPI[1]     = (body*) malloc(p->nBody * sizeof(body));

	double localDt, dt;
	writeOuputFile(0, p);
	gettimeofday(&t1,NULL);

	if(!RankMPI)
		printf("Starting simulation...\n");
	for(unsigned long iIte = 1; iIte <= NIterations; ++iIte) {
		/*******************************/
		/*** Simulation computations ***/
		gettimeofday(&t2,NULL);

		// copy plan bodies into MPI buffer
		for(unsigned long iBody = 0; iBody < p->nBody; ++iBody)
			initBody(&(bodyBufferMPI[0][iBody]), p->lb[iBody].b->mass, p->lb[iBody].b->posX, p->lb[iBody].b->posY);

		// compute local bodies acceleration with local bodies
		computeAllLocalAcceleration(p);

		// compute local bodies acceleration with all bodies from others processes
		for(int iStep = 1; iStep < SizeMPI; ++iStep)
		{
			// double buffering
			if(RankMPI % 2)
			{
				MPI_Recv(bodyBufferMPI[iStep      %2], p->nBody, bodyMPI, prevRankMPI, 1234, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
				MPI_Send(bodyBufferMPI[(iStep +1) %2], p->nBody, bodyMPI, nextRankMPI, 1234, MPI_COMM_WORLD);
			}
			else
			{
				MPI_Send(bodyBufferMPI[(iStep +1) %2], p->nBody, bodyMPI, nextRankMPI, 1234, MPI_COMM_WORLD);
				MPI_Recv(bodyBufferMPI[iStep      %2], p->nBody, bodyMPI, prevRankMPI, 1234, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			}
			computeAllAcceleration(p, bodyBufferMPI[iStep %2]);
		}

		localDt = findLocalDt(p);
		MPI_Allreduce(&localDt, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

		updateAllLocalPositionAndSpeed(p, dt);

		gettimeofday(&t3,NULL);
		/*** Simulation computations ***/
		/*******************************/

		printIterationTimeAndPerformance(iIte, t1, t2, t3);
		writeOuputFile(iIte, p);
	}
	if(!RankMPI)
		printf("Simulation end... exiting.\n");

	MPI_Barrier(MPI_COMM_WORLD);

	// release allocations
	free(bodyBufferMPI[0]);
	free(bodyBufferMPI[1]);
	free(bodyBufferMPI);
	destroyPlan(p);

	// MPI deinit
	MPI_Type_free(&bodyMPI);
	MPI_Finalize();

	return EXIT_SUCCESS;
}
