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

#include "SimulationNBodyCollisionV1.h"

template <typename T>
SimulationNBodyCollisionV1<T>::SimulationNBodyCollisionV1(const unsigned long nBodies)
	: SimulationNBodyCollisionLocal<T>(nBodies)
{
	this->init();
}

template <typename T>
SimulationNBodyCollisionV1<T>::SimulationNBodyCollisionV1(const std::string inputFileName)
	: SimulationNBodyCollisionLocal<T>(inputFileName)
{
	this->init();
}

template <typename T>
void SimulationNBodyCollisionV1<T>::init()
{
	this->flopsPerIte = 20 * (this->bodies.getN() -1) * this->bodies.getN();
}

template <typename T>
SimulationNBodyCollisionV1<T>::~SimulationNBodyCollisionV1()
{
}

template <typename T>
void SimulationNBodyCollisionV1<T>::initIteration()
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
void SimulationNBodyCollisionV1<T>::computeLocalBodiesAcceleration()
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
void SimulationNBodyCollisionV1<T>::computeAccelerationBetweenTwoBodiesNaive(const unsigned long &iBody, const unsigned long &jBody)
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
	const T dist = std::sqrt((diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ)); // 6 flops

	// we cannot divide by 0
	if(dist == 0)
	{
		std::cout << "Collision at {" << positionsX[jBody] << ", "
		                              << positionsY[jBody] << ", "
		                              << positionsZ[jBody] << "}" << std::endl;
		assert(dist != 0);
	}

	// detect collisions
	/*
	std::cout << "body n°" << iBody << ": "
	          << "dist = " << dist << ", "
	          << "dist - (radiuses[iBody] + radiuses[jBody]) = " << dist - (radiuses[iBody] + radiuses[jBody])
	          << std::endl;
	*/
	if(dist - (radiuses[iBody] + radiuses[jBody]) <= 0) // 2 flops
	{
		/*
		std::cout << "Collision detected at {" << positionsX[jBody] << ", "
		                                       << positionsY[jBody] << ", "
		                                       << positionsZ[jBody] << "} "
		          << "between body n°" << iBody << " and body n°" << jBody << "." << std::endl;
		*/

		this->collisions[iBody].push_back(jBody);
	}

	// compute the force value between iBody and jBody: || F || = G.mi.mj / Dij²
	const T force = (this->G * masses[iBody] * masses[jBody] / (dist * dist)); // 4 flops

	// compute the acceleration value: || a || = || F || / mi
	const T acc = force / masses[iBody]; // 1 flop

	/*
	std::cout << "body n°" << iBody << ": "
	          << "positionsX = " << positionsX[iBody] << ", "
	          << "positionsY = " << positionsY[iBody] << ", "
	          << "positionsZ = " << positionsZ[iBody] << ", "
	          << "dist = " << dist << ", "
	          << "acc = " << acc
 	          << std::endl;
	*/

	// normalize and add acceleration value into acceleration vector: a += || a ||.u
	this->accelerations.x[iBody] += acc * (diffPosX / dist); // 3 flops
	this->accelerations.y[iBody] += acc * (diffPosY / dist); // 3 flops
	this->accelerations.z[iBody] += acc * (diffPosZ / dist); // 3 flops

	if(!this->dtConstant)
		if(dist < this->closestNeighborDist[iBody])
			this->closestNeighborDist[iBody] = dist;
}

// 20 flops
template <typename T>
void SimulationNBodyCollisionV1<T>::computeAccelerationBetweenTwoBodies(const unsigned long &iBody, const unsigned long &jBody)
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
	const T dist = std::sqrt(squareDist); // 1 flop
	assert(dist != 0);

	if(dist - (radiuses[iBody] + radiuses[jBody]) <= 0) // 2 flops
		this->collisions[iBody].push_back(jBody);

	const T acc = this->G * masses[jBody] / (squareDist * dist); // 3 flops
	this->accelerations.x[iBody] += acc * diffPosX; // 2 flop
	this->accelerations.y[iBody] += acc * diffPosY; // 2 flop
	this->accelerations.z[iBody] += acc * diffPosZ; // 2 flop

	if(!this->dtConstant)
		if(dist < this->closestNeighborDist[iBody])
			this->closestNeighborDist[iBody] = dist;
}
