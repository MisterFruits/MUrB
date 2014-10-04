/*
 * Do not remove.
 * MPI/OpenMP training courses
 * Adrien Cassagne, ASA - CINES, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#include "math.h"
#include "float.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "assert.h"

#include "plan.h"

extern double G;
extern int    RankMPI;

plan* createPlan(const unsigned long nBody)
{
	plan *p  = (plan*)malloc(1 * sizeof(plan));
	p->nBody = nBody;
	p->lb    = (localBody*)malloc(nBody * sizeof(localBody));
	return p;
}

plan* createPlanWithFile(const char *rootFileName)
{
	plan *p = NULL;

	char intFileName[1024];
	sprintf(intFileName, "%s.%d.dat", rootFileName, RankMPI);
	FILE *file = fopen(intFileName, "r");

	if(file == NULL)
		printf("Can't open \"%s\" file (reading).\n", intFileName);

	unsigned long nBody = 0;
	double mass, posX, posY, speedX, speedY;
	while(fscanf(file, "%lf %lf %lf %lf %lf\n", &mass, &posX, &posY, &speedX, &speedY) != EOF)
		nBody++;
	fseek(file, 0, SEEK_SET);

	if(nBody > 0)
	{
		p = createPlan(nBody);
		unsigned long iBody = 0;
		while(fscanf(file, "%lf %lf %lf %lf %lf\n", &mass, &posX, &posY, &speedX, &speedY) != EOF)
			initLocalBody(&(p->lb[iBody++]), mass, posX, posY, speedX, speedY);
	}
	else
		printf("Empty \"%s\" file (reading).\n", intFileName);

	fclose(file);

	return p;
}

void destroyPlan(plan *p)
{
	for(unsigned long iBody = 0; iBody < p->nBody; ++iBody)
		free(p->lb[iBody].b);

	free(p->lb);
	free(p);
}

// total flops ~= nBody^2 * 12
void computeAllLocalAcceleration(plan *p)
{
	// flops ~= nBody^2 * 12
	for(unsigned long iBody = 0; iBody < p->nBody; ++iBody)
		// flops ~= nBody * 12
		for(unsigned long jBody = 0; jBody < p->nBody; ++jBody)
			if(iBody != jBody)
				computeAcceleration(&(p->lb[iBody]), p->lb[jBody].b); // 12 flops
}

// total flops = nBody * 12
double findLocalDt(const plan *p)
{
	double dt = DBL_MAX;
	// flops = nBody * 12
	for(unsigned long iBody = 0; iBody < p->nBody; ++iBody)
	{
		const double newDt = computeDt(p->lb[iBody]); // 12 flops
		if(newDt < dt)
			dt = newDt;
	}

	return dt;
}

// total flops = nBody * 12
void updateAllLocalPositionAndSpeed(plan *p, double dt)
{
	// flops = nBody * 12
	for(unsigned long iBody = 0; iBody < p->nBody; ++iBody)
		updatePositionAndSpeed(&(p->lb[iBody]), dt); // 12 flops
}

void fillRandom(plan *p)
{
	srand(123 * RankMPI);
	for(unsigned long iBody = 0; iBody < p->nBody; ++iBody)
	{
		const double mass   = (rand() / (double) RAND_MAX) * 2000000;
		const double posX   = (rand() / (double) RAND_MAX) * 800;
		const double posY   = (rand() / (double) RAND_MAX) * 600;
		const double speedX = ((rand() - RAND_MAX/2) / (double) (RAND_MAX/2)) * 0.02;
		const double speedY = ((rand() - RAND_MAX/2) / (double) (RAND_MAX/2)) * 0.02;

		initLocalBody(&(p->lb[iBody]), mass, posX, posY, speedX, speedY);
	}
}

void writePlanIntoFile(const plan *p, const char *fileName)
{
	assert(p->nBody > 0);

	FILE *file = NULL;
	if(RankMPI == 0)
		file = fopen(fileName, "w");
	else
		file = fopen(fileName, "a");

	if(file == NULL)
	{
		printf("Can't open \"%s\" file (writing). Exiting...\n", fileName);
		exit(0);
	}

	for(unsigned long iBody = 0; iBody < p->nBody; ++iBody)
		fprintf(file, "%lf %lf %lf %lf %lf\n", p->lb[iBody].b->mass / G,
		                                       p->lb[iBody].b->posX,
		                                       p->lb[iBody].b->posY,
		                                       p->lb[iBody].speedX,
		                                       p->lb[iBody].speedY);

	fclose(file);
}
