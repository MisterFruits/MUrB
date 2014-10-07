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
#ifndef _MYOPENMP
#define _MYOPENMP
inline void omp_set_num_threads(int) {           }
inline int  omp_get_num_threads(   ) { return 1; }
inline int  omp_get_max_threads(   ) { return 1; }
inline int  omp_get_thread_num (   ) { return 0; }
#endif
#endif

#include "SimulationNBodyV2.h"

template <typename T>
SimulationNBodyV2<T>::SimulationNBodyV2(const unsigned long nBodies)
	: SimulationNBody<T>(nBodies)
{
	this->allocateBuffers();
}

template <typename T>
SimulationNBodyV2<T>::SimulationNBodyV2(const std::string inputFileName)
	: SimulationNBody<T>(inputFileName)
{
	this->allocateBuffers();
}

template <typename T>
void SimulationNBodyV2<T>::allocateBuffers()
{
	this->accelerations.x = new T[this->bodies.getN() * omp_get_max_threads()];
	this->accelerations.y = new T[this->bodies.getN() * omp_get_max_threads()];
	this->accelerations.z = new T[this->bodies.getN() * omp_get_max_threads()];

	this->closestNeighborDist = new T[this->bodies.getN()];
}

template <typename T>
SimulationNBodyV2<T>::~SimulationNBodyV2()
{
}

template <typename T>
void SimulationNBodyV2<T>::initIteration()
{
	for(unsigned long iBody = 0; iBody < this->bodies.getN() * omp_get_max_threads(); iBody++)
	{
		this->accelerations.x[iBody] = 0.0;
		this->accelerations.y[iBody] = 0.0;
		this->accelerations.z[iBody] = 0.0;
	}

	for(unsigned long iBody = 0; iBody < this->bodies.getN(); iBody++)
		this->closestNeighborDist[iBody] = std::numeric_limits<T>::infinity();
}

template <typename T>
void SimulationNBodyV2<T>::computeBodiesAcceleration()
{
#pragma omp parallel for schedule(runtime)
	for(unsigned long iBody = 0; iBody < this->bodies.getN(); iBody++)
		for(unsigned long jBody = iBody +1; jBody < this->bodies.getN(); jBody++)
			//this->computeAccelerationBetweenTwoBodiesNaive(iBody, jBody);
			this->computeAccelerationBetweenTwoBodies(iBody, jBody);

	if(omp_get_max_threads() > 1)
		for(unsigned long iBody = 0; iBody < this->bodies.getN(); iBody++)
			for(unsigned short iThread = 1; iThread < omp_get_max_threads(); iThread++)
			{
				this->accelerations.x[iBody] += this->accelerations.x[iBody + iThread * this->bodies.getN()];
				this->accelerations.y[iBody] += this->accelerations.y[iBody + iThread * this->bodies.getN()];
				this->accelerations.z[iBody] += this->accelerations.z[iBody + iThread * this->bodies.getN()];
			}
}

// 32 flops
template <typename T>
void SimulationNBodyV2<T>::computeAccelerationBetweenTwoBodiesNaive(const unsigned long iBody, const unsigned long jBody)
{
	assert(iBody != jBody);

	const T *masses     = this->bodies.getMasses();
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

	// compute the force value between iBody and jBody: || F || = G.mi.mj / Dij²
	const T force = this->G * masses[iBody] * masses[jBody] / (dist * dist); // 4 flops

	// compute the acceleration value: || a || = || F || / mi
	T acc = force / masses[iBody]; // 1 flop

	// normalize and add acceleration value into acceleration vector: a += || a ||.u
	this->accelerations.x[iBody + omp_get_thread_num() * this->bodies.getN()] += acc * (diffPosX / dist); // 3 flops
	this->accelerations.y[iBody + omp_get_thread_num() * this->bodies.getN()] += acc * (diffPosY / dist); // 3 flops
	this->accelerations.z[iBody + omp_get_thread_num() * this->bodies.getN()] += acc * (diffPosZ / dist); // 3 flops

	// compute the acceleration value: || a || = || F || / mj
	acc = force / masses[jBody]; // 1 flop

	this->accelerations.x[jBody + omp_get_thread_num() * this->bodies.getN()] -= acc * (diffPosX / dist); // 3 flops
	this->accelerations.y[jBody + omp_get_thread_num() * this->bodies.getN()] -= acc * (diffPosY / dist); // 3 flops
	this->accelerations.z[jBody + omp_get_thread_num() * this->bodies.getN()] -= acc * (diffPosZ / dist); // 3 flops

	if(!this->dtConstant)
	{
		if(dist < this->closestNeighborDist[iBody])
			this->closestNeighborDist[iBody] = dist;

		if(dist < this->closestNeighborDist[jBody])
#pragma omp critical
			if(dist < this->closestNeighborDist[jBody])
				this->closestNeighborDist[jBody] = dist;
	}
}

// 25 flops
template <typename T>
void SimulationNBodyV2<T>::computeAccelerationBetweenTwoBodies(const unsigned long iBody, const unsigned long jBody)
{
	assert(iBody != jBody);

	const T *masses     = this->bodies.getMasses();
	const T *positionsX = this->bodies.getPositionsX();
	const T *positionsY = this->bodies.getPositionsY();
	const T *positionsZ = this->bodies.getPositionsZ();

	const T diffPosX  = positionsX[jBody] - positionsX[iBody]; // 1 flop
	const T diffPosY  = positionsY[jBody] - positionsY[iBody]; // 1 flop
	const T diffPosZ  = positionsZ[jBody] - positionsZ[iBody]; // 1 flop

	const T squareDist = (diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ); // 5 flops
	const T dist = std::sqrt(squareDist); // 1 flops
	assert(dist != 0);

	const T force = this->G / (squareDist * dist); // 2 flops

	T acc = force * masses[jBody]; // 1 flop
	this->accelerations.x[iBody + omp_get_thread_num() * this->bodies.getN()] += acc * diffPosX; // 2 flops
	this->accelerations.y[iBody + omp_get_thread_num() * this->bodies.getN()] += acc * diffPosY; // 2 flops
	this->accelerations.z[iBody + omp_get_thread_num() * this->bodies.getN()] += acc * diffPosZ; // 2 flops

	acc = force * masses[iBody]; // 1 flop
	this->accelerations.x[jBody + omp_get_thread_num() * this->bodies.getN()] -= acc * diffPosX; // 2 flops
	this->accelerations.y[jBody + omp_get_thread_num() * this->bodies.getN()] -= acc * diffPosY; // 2 flops
	this->accelerations.z[jBody + omp_get_thread_num() * this->bodies.getN()] -= acc * diffPosZ; // 2 flops

	if(!this->dtConstant)
	{
		if(dist < this->closestNeighborDist[iBody])
			this->closestNeighborDist[iBody] = dist;

		if(dist < this->closestNeighborDist[jBody])
#pragma omp critical
			if(dist < this->closestNeighborDist[jBody])
				this->closestNeighborDist[jBody] = dist;
	}
}
