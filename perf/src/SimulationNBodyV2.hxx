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
	: SimulationNBody<T>(nBodies), nMaxThreads(omp_get_max_threads())
{
	this->reAllocateBuffers();
}

template <typename T>
SimulationNBodyV2<T>::SimulationNBodyV2(const std::string inputFileName)
	: SimulationNBody<T>(inputFileName) , nMaxThreads(omp_get_max_threads())
{
	this->reAllocateBuffers();
}

template <typename T>
SimulationNBodyV2<T>::~SimulationNBodyV2()
{
	if(this->accelerations.x != nullptr) {
		delete[] this->accelerations.x;
		this->accelerations.x = nullptr;
	}
	if(this->accelerations.y != nullptr) {
		delete[] this->accelerations.y;
		this->accelerations.y = nullptr;
	}
	if(this->accelerations.z != nullptr) {
		delete[] this->accelerations.z;
		this->accelerations.z = nullptr;
	}
}

template <typename T>
void SimulationNBodyV2<T>::reAllocateBuffers()
{
	// TODO: this is not optimal to deallocate and to reallocate data
	if(this->nMaxThreads > 1)
		if(this->accelerations.x != nullptr)
			delete[] this->accelerations.x;
		if(this->accelerations.y != nullptr)
			delete[] this->accelerations.y;
		if(this->accelerations.z != nullptr)
			delete[] this->accelerations.z;

		this->accelerations.x = new T[this->bodies.getN() * this->nMaxThreads];
		this->accelerations.y = new T[this->bodies.getN() * this->nMaxThreads];
		this->accelerations.z = new T[this->bodies.getN() * this->nMaxThreads];
}

template <typename T>
void SimulationNBodyV2<T>::initIteration()
{
	for(unsigned long iBody = 0; iBody < this->bodies.getN() * this->nMaxThreads; iBody++)
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
#pragma omp parallel
{
	const unsigned tid = omp_get_thread_num();

#pragma omp for schedule(runtime)
	for(unsigned long iBody = 0; iBody < this->bodies.getN(); iBody++)
		for(unsigned long jBody = iBody +1; jBody < this->bodies.getN(); jBody++)
			//this->computeAccelerationBetweenTwoBodiesNaive(iBody, jBody, tid);
			this->computeAccelerationBetweenTwoBodies(iBody, jBody, tid);
}

	if(this->nMaxThreads > 1)
	{
		for(unsigned long iBody = 0; iBody < this->bodies.getN(); iBody++)
			for(unsigned iThread = 1; iThread < this->nMaxThreads; iThread++)
			{
				this->accelerations.x[iBody] += this->accelerations.x[iBody + iThread * this->bodies.getN()];
				this->accelerations.y[iBody] += this->accelerations.y[iBody + iThread * this->bodies.getN()];
				this->accelerations.z[iBody] += this->accelerations.z[iBody + iThread * this->bodies.getN()];
			}
	}
}

// 32 flops
template <typename T>
void SimulationNBodyV2<T>::computeAccelerationBetweenTwoBodiesNaive(const unsigned long &iBody,
                                                                    const unsigned long &jBody,
                                                                    const unsigned      &tid)
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
	this->accelerations.x[iBody + tid * this->bodies.getN()] += acc * (diffPosX / dist); // 3 flops
	this->accelerations.y[iBody + tid * this->bodies.getN()] += acc * (diffPosY / dist); // 3 flops
	this->accelerations.z[iBody + tid * this->bodies.getN()] += acc * (diffPosZ / dist); // 3 flops

	// compute the acceleration value: || a || = || F || / mj
	acc = force / masses[jBody]; // 1 flop

	this->accelerations.x[jBody + tid * this->bodies.getN()] -= acc * (diffPosX / dist); // 3 flops
	this->accelerations.y[jBody + tid * this->bodies.getN()] -= acc * (diffPosY / dist); // 3 flops
	this->accelerations.z[jBody + tid * this->bodies.getN()] -= acc * (diffPosZ / dist); // 3 flops

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
void SimulationNBodyV2<T>::computeAccelerationBetweenTwoBodies(const unsigned long &iBody,
                                                               const unsigned long &jBody,
                                                               const unsigned      &tid)
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

	const unsigned long idIBody = iBody + tid* this->bodies.getN();
	const unsigned long idJBody = jBody + tid* this->bodies.getN();


	this->accelerations.x[idIBody] += acc * diffPosX; // 2 flops
	this->accelerations.y[idIBody] += acc * diffPosY; // 2 flops
	this->accelerations.z[idIBody] += acc * diffPosZ; // 2 flops

	acc = force * masses[iBody]; // 1 flop
	this->accelerations.x[idJBody] -= acc * diffPosX; // 2 flops
	this->accelerations.y[idJBody] -= acc * diffPosY; // 2 flops
	this->accelerations.z[idJBody] -= acc * diffPosZ; // 2 flops

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
