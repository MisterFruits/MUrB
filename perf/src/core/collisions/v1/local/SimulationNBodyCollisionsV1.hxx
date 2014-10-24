/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#include <cmath>
#include <limits>
#include <string>
#include <cassert>
#include <fstream>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#else
#ifndef NO_OMP
#define NO_OMP
inline void omp_set_num_threads(int) {           }
inline int  omp_get_num_threads(   ) { return 1; }
inline int  omp_get_max_threads(   ) { return 1; }
inline int  omp_get_thread_num (   ) { return 0; }
#endif
#endif

#include "SimulationNBodyCollisionsV1.h"

template <typename T>
SimulationNBodyCollisionsV1<T>::SimulationNBodyCollisionsV1(const unsigned long nBodies)
	: SimulationNBodyCollisionsLocal<T>(nBodies)
{
	this->init();
}

template <typename T>
SimulationNBodyCollisionsV1<T>::SimulationNBodyCollisionsV1(const std::string inputFileName)
	: SimulationNBodyCollisionsLocal<T>(inputFileName)
{
	this->init();
}

template <typename T>
void SimulationNBodyCollisionsV1<T>::init()
{
	this->flopsPerIte = 20 * (this->bodies.getN() -1) * this->bodies.getN();
}

template <typename T>
SimulationNBodyCollisionsV1<T>::~SimulationNBodyCollisionsV1()
{
}

template <typename T>
void SimulationNBodyCollisionsV1<T>::initIteration()
{
	for(unsigned long iBody = 0; iBody < this->bodies.getN(); iBody++)
	{
		this->accelerations.x[iBody] = 0.0;
		this->accelerations.y[iBody] = 0.0;
		this->accelerations.z[iBody] = 0.0;

		this->closestNeighborDist[iBody] = std::numeric_limits<T>::infinity();

		this->collisions[iBody].clear();
		//if(this->collisions[iBody].capacity() < 4)
		//	this->collisions[iBody].resize(4);
	}
}

template <typename T>
void SimulationNBodyCollisionsV1<T>::computeLocalBodiesAcceleration()
{
#pragma omp parallel for schedule(runtime)
	for(unsigned long iBody = 0; iBody < this->bodies.getN(); iBody++)
		for(unsigned long jBody = 0; jBody < this->bodies.getN(); jBody++)
			if(iBody != jBody)
				//this->computeAccelerationBetweenTwoBodiesNaive(iBody, jBody);
				this->computeAccelerationBetweenTwoBodies(iBody, jBody);
}

// 25 flops
template <typename T>
void SimulationNBodyCollisionsV1<T>::computeAccelerationBetweenTwoBodiesNaive(const unsigned long &iBody, const unsigned long &jBody)
{
	assert(iBody != jBody);

	const T *masses     = this->bodies.getMasses();
	const T *radiuses   = this->bodies.getRadiuses();
	const T *positionsX = this->bodies.getPositionsX();
	const T *positionsY = this->bodies.getPositionsY();
	const T *positionsZ = this->bodies.getPositionsZ();

	const T diffPosX = positionsX[jBody] - positionsX[iBody]; // 1 flop
	const T diffPosY = positionsY[jBody] - positionsY[iBody]; // 1 flop
	const T diffPosZ = positionsZ[jBody] - positionsZ[iBody]; // 1 flop

	// compute distance between iBody and jBody: Dij
	const T dij = std::sqrt((diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ)); // 6 flops

	// we cannot divide by 0
	if(dij == 0)
	{
		std::cout << "Collision at {" << positionsX[jBody] << ", "
		                              << positionsY[jBody] << ", "
		                              << positionsZ[jBody] << "}" << std::endl;
		assert(dij != 0);
	}

	// detect collisions
	const T dist = dij - (radiuses[iBody] + radiuses[jBody]); // 2 flops
	if(dist <= 0)
		this->collisions[iBody].push_back(jBody);

	// compute the force value between iBody and jBody: || F || = G.mi.mj / Dij²
	const T force = (this->G * masses[iBody] * masses[jBody] / (dij * dij)); // 4 flops

	// compute the acceleration value: || a || = || F || / mi
	const T acc = force / masses[iBody]; // 1 flop

	// normalize and add acceleration value into acceleration vector: a += || a ||.u
	this->accelerations.x[iBody] += acc * (diffPosX / dij); // 3 flops
	this->accelerations.y[iBody] += acc * (diffPosY / dij); // 3 flops
	this->accelerations.z[iBody] += acc * (diffPosZ / dij); // 3 flops

	if(!this->dtConstant)
		if(dist < this->closestNeighborDist[iBody])
			this->closestNeighborDist[iBody] = dist;
}

// 20 flops
template <typename T>
void SimulationNBodyCollisionsV1<T>::computeAccelerationBetweenTwoBodies(const unsigned long &iBody, const unsigned long &jBody)
{
	assert(iBody != jBody);

	const T *masses     = this->bodies.getMasses();
	const T *radiuses   = this->bodies.getRadiuses();
	const T *positionsX = this->bodies.getPositionsX();
	const T *positionsY = this->bodies.getPositionsY();
	const T *positionsZ = this->bodies.getPositionsZ();

	const T diffPosX = positionsX[jBody] - positionsX[iBody]; // 1 flop
	const T diffPosY = positionsY[jBody] - positionsY[iBody]; // 1 flop
	const T diffPosZ = positionsZ[jBody] - positionsZ[iBody]; // 1 flop
	const T squareDist = (diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ); // 5 flops
	const T dij = std::sqrt(squareDist); // 1 flop
	assert(dij != 0);

	const T dist = dij - (radiuses[iBody] + radiuses[jBody]); // 2 flops
	if(dist <= 0)
		this->collisions[iBody].push_back(jBody);

	const T acc = this->G * masses[jBody] / (squareDist * dij); // 3 flops
	this->accelerations.x[iBody] += acc * diffPosX; // 2 flop
	this->accelerations.y[iBody] += acc * diffPosY; // 2 flop
	this->accelerations.z[iBody] += acc * diffPosZ; // 2 flop

	if(!this->dtConstant)
		if(dist < this->closestNeighborDist[iBody])
			this->closestNeighborDist[iBody] = dist;
}
