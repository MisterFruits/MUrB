/*
 * Do not remove.
 * MPI/OpenMP training courses
 * Adrien CASSAGNE, ASA - CINES, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#include "math.h"
#include "stdio.h"
#include "float.h"
#include "assert.h"
#include "stdlib.h"

#include "body.h"

double G = 6.67384e-11;

void initLocalBody(localBody* lb, double mass, double posX, double posY, double speedX, double speedY)
{
	assert(lb  != NULL);
	assert(mass > 0);

	lb->speedX             = speedX;
	lb->speedY             = speedY;
	lb->accelerationX      = 0;
	lb->accelerationY      = 0;
	lb->closestNeighborLen = DBL_MAX;
	lb->b                  = (body*)malloc(1 * sizeof(body));
	lb->b->mass            = G * mass;
	lb->b->posX            = posX;
	lb->b->posY            = posY;
}

localBody* createLocalBody(double mass, double posX, double posY, double speedX, double speedY)
{
	localBody *lb = (localBody*)malloc(1 * sizeof(localBody));
	initLocalBody(lb, mass, posX, posY, speedX, speedY);
	return lb;
}

void destroyLocalBody(localBody *lb)
{
	assert(lb != NULL);
	free(lb->b);
	free(lb);
}

// total flops = 12
void computeAcceleration(localBody *lb, const body *b)
{

	const double vecX   = (b->posX - lb->b->posX); // 1 flop
	const double vecY   = (b->posY - lb->b->posY); // 1 flop
	const double vecLen = sqrt((vecX * vecX) + (vecY * vecY)); // 3 flops
	
	if(vecLen == 0) 
		printf("Collision at {%f, %f}\n", lb->b->posX, lb->b->posY);
	assert(vecLen != 0);
	
	const double acc  = b->mass / (vecLen * vecLen * vecLen); // 3 flops
	const double accX = acc * vecX; // 1 flop
	const double accY = acc * vecY; // 1 flop

	lb->accelerationX += accX; // 1 flop
	lb->accelerationY += accY; // 1 flop

	if(vecLen < lb->closestNeighborLen)
		lb->closestNeighborLen = vecLen;
}

// total flops = 12
double computeDt(const localBody lb)
{
	/* || lb.speed ||        */
	const double s = sqrt((lb.speedX * lb.speedX) + (lb.speedY * lb.speedY)); // 3 flops
	
	/* || lb.acceleration || */
	const double a = sqrt((lb.accelerationX * lb.accelerationX) + (lb.accelerationY * lb.accelerationY)); // 3 flops
	
	/* 
	 * compute dt
	 * solve:  (a/2)*dt^2 + s*dt + (-0.1)*ClosestNeighborLen = 0
	 * <=>     dt = [ (-s) +/-  sqrt( s^2 - 4 * (a/2) * (-0.1)*ClosestNeighborLen ) ] / [ 2 (a/2) ]
	 *
	 * dt should be positive (+/- becomes + because result of sqrt is positive)
	 * <=>     dt = [ -s + sqrt( s^2 + 0.2*ClosestNeighborLen*a) ] / a
	 */
	
	double dt = (sqrt(s * s + 0.2 * a * lb.closestNeighborLen) - s) / a; // 6 flops
	
	if(dt == 0)
		dt = DBL_EPSILON / a;

	return dt;
}

// total flops = 12
void updatePositionAndSpeed(localBody *lb, const double dt)
{
	double accXMultDt = lb->accelerationX * dt; // 1 flop
	double accYMultDt = lb->accelerationY * dt; // 1 flop

	lb->b->posX += (lb->speedX + accXMultDt / 2.0) * dt; // 4 flops
	lb->b->posY += (lb->speedY + accYMultDt / 2.0) * dt; // 4 flops

	lb->speedX += accXMultDt; // 1 flop
	lb->speedY += accYMultDt; // 1 flop

	lb->accelerationX      = 0;
	lb->accelerationY      = 0;
	lb->closestNeighborLen = DBL_MAX;
}
