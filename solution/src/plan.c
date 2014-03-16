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

#define MIN(a,b) (a<=b?a:b)

plan* createPlan(const unsigned long nBody)
{
	plan *p  = (plan*)malloc(1 * sizeof(plan));
	p->nBody = nBody;
	p->lb    = (localBody*)malloc(nBody * sizeof(localBody));
	return p;
}

plan* createPlanWithFile(const char *fileName)
{
	plan *p = NULL;
	FILE *file = fopen(fileName, "r");

	if(file == NULL)
		printf("Can't open \"%s\" file (reading).\n", fileName);

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
		printf("Empty \"%s\" file (reading).\n", fileName);

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
	unsigned long blockSize = 4096;
	unsigned long nBlock    = ceil(p->nBody / (blockSize * 1.0));

#pragma omp parallel for schedule(runtime)
	for(unsigned long iBlock = 0; iBlock < nBlock; ++iBlock)
	{
		const unsigned long realBlockSize = MIN(blockSize, p->nBody - iBlock * blockSize);

		for(unsigned long iBody = 0; iBody < p->nBody; ++iBody)
			for(unsigned long jBody = iBlock * blockSize; jBody < iBlock * blockSize + realBlockSize; ++jBody)
				if(iBody != jBody)
					computeAcceleration(&(p->lb[iBody]), p->lb[jBody].b);
	}
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
	srand(123);
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

	FILE *file = fopen(fileName, "w");

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
