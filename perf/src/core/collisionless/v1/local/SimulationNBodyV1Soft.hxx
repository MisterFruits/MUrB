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

#include "SimulationNBodyV1Soft.h"

template <typename T>
SimulationNBodyV1Soft<T>::SimulationNBodyV1Soft(const unsigned long nBodies, T softening)
	: SimulationNBodyLocal<T>(nBodies), softeningSquared(softening * softening)
{
	this->init();
}

template <typename T>
SimulationNBodyV1Soft<T>::SimulationNBodyV1Soft(const std::string inputFileName, T softening)
	: SimulationNBodyLocal<T>(inputFileName), softeningSquared(softening * softening)
{
	this->init();
}

template <typename T>
void SimulationNBodyV1Soft<T>::init()
{
	this->flopsPerIte = 19 * (this->bodies->getN() -1) * this->bodies->getN();
}

template <typename T>
SimulationNBodyV1Soft<T>::~SimulationNBodyV1Soft()
{
}

template <typename T>
void SimulationNBodyV1Soft<T>::initIteration()
{
	for(unsigned long iBody = 0; iBody < this->bodies->getN(); iBody++)
	{
		this->accelerations.x[iBody] = 0.0;
		this->accelerations.y[iBody] = 0.0;
		this->accelerations.z[iBody] = 0.0;

		this->closestNeighborDist[iBody] = std::numeric_limits<T>::infinity();
	}
}

template <typename T>
void SimulationNBodyV1Soft<T>::computeLocalBodiesAcceleration()
{
#pragma omp parallel for schedule(runtime)
	for(unsigned long iBody = 0; iBody < this->bodies->getN(); iBody++)
		for(unsigned long jBody = 0; jBody < this->bodies->getN(); jBody++)
			//this->computeAccelerationBetweenTwoBodiesNaive(iBody, jBody);
			this->computeAccelerationBetweenTwoBodies(iBody, jBody);
}

// 19 flops
template <typename T>
void SimulationNBodyV1Soft<T>::computeAccelerationBetweenTwoBodies(const unsigned long &iBody, const unsigned long &jBody)
{
	assert(iBody != jBody);

	const T *masses     = this->bodies->getMasses();
	const T *positionsX = this->bodies->getPositionsX();
	const T *positionsY = this->bodies->getPositionsY();
	const T *positionsZ = this->bodies->getPositionsZ();

	const T rIJX = positionsX[jBody] - positionsX[iBody]; // 1 flop
	const T rIJY = positionsY[jBody] - positionsY[iBody]; // 1 flop
	const T rIJZ = positionsZ[jBody] - positionsZ[iBody]; // 1 flop
	const T squareDist = (rIJX * rIJX) + (rIJY * rIJY) + (rIJZ * rIJZ) + this->softeningSquared; // 6 flops
	const T dist = std::sqrt(squareDist); // 1 flop
	assert(dist != 0);

	const T acc = this->G * masses[jBody] / (squareDist * dist); // 3 flops

	this->accelerations.x[iBody] += acc * rIJX; // 2 flop
	this->accelerations.y[iBody] += acc * rIJY; // 2 flop
	this->accelerations.z[iBody] += acc * rIJZ; // 2 flop

	if(!this->dtConstant)
		if(dist < this->closestNeighborDist[iBody])
			this->closestNeighborDist[iBody] = dist;
}
