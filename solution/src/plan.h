/*
 * Do not remove.
 * MPI/OpenMP training courses
 * Adrien Cassagne, ASA - CINES, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef PLAN_H_
#define PLAN_H_

#include "body.h"

typedef struct
{
	unsigned long nBody ; // bodies number
	localBody     *lb;    // local body array values
} plan;

/* Allocate and return a plan with bodies randomly created */
plan* createPlan(const unsigned long nBody);

/* Allocate and return a plan with bodies from a file */
plan* createPlanWithFile(const char *fileName);

/* Deallocate and destroy a plan (and all bodies in the plan) */
void destroyPlan(plan *p);

/* Compute acceleration for all bodies in local plan with local bodies */
void computeAllLocalAcceleration(plan *p);

/* Compute acceleration for all bodies in plan with others bodies from others processes */
void computeAllAcceleration(plan *p, body *b);

/* Return the lowest dt between all bodies in plan */
double findLocalDt(const plan *p);

/* Compute new position and speed for all bodies in plan */
void updateAllLocalPositionAndSpeed(plan *p, double dt);

/* Fill plan with random bodies */
void fillRandom(plan *p);

/* Write plan (and bodies) in file */
void writePlanIntoFile(const plan *p, const char *fileName);

#endif /* PLAN_H_ */
