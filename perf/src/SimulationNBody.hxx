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
inline void omp_set_num_threads(int) {           }
inline int  omp_get_num_threads(   ) { return 1; }
inline int  omp_get_max_threads(   ) { return 1; }
inline int  omp_get_thread_num (   ) { return 0; }
#endif

#include "SimulationNBody.h"

template <typename T>
SimulationNBody<T>::SimulationNBody(const unsigned long nBodies)
	: bodies(nBodies),
	  closestNeighborDist(NULL),
	  dt                 (std::numeric_limits<T>::infinity()),
	  dtConstant         (false)
{
	assert(nBodies > 0);
	allocateBuffers();
}

template <typename T>
SimulationNBody<T>::SimulationNBody(const std::string inputFileName)
	: bodies             (inputFileName),
	  closestNeighborDist(NULL),
	  dt                 (std::numeric_limits<T>::infinity()),
	  dtConstant         (false)
{
	allocateBuffers();
}

template <typename T>
void SimulationNBody<T>::allocateBuffers()
{
	this->accelerations.x = new T[this->bodies.n * omp_get_max_threads()];
	this->accelerations.y = new T[this->bodies.n * omp_get_max_threads()];
	this->accelerations.z = new T[this->bodies.n * omp_get_max_threads()];

	this->closestNeighborDist = new T[this->bodies.n];
}

template <typename T>
SimulationNBody<T>::~SimulationNBody() {
	if(this->accelerations.x)
		delete[] this->accelerations.x;
	if(this->accelerations.y)
		delete[] this->accelerations.y;
	if(this->accelerations.z)
		delete[] this->accelerations.z;

	if(this->closestNeighborDist)
		delete[] this->closestNeighborDist;
}

template <typename T>
void SimulationNBody<T>::setDtConstant(T dtVal)
{
	this->dtConstant = true;
	this->dt = dtVal;
}

template <typename T>
void SimulationNBody<T>::setDtVariable()
{
	this->dtConstant = false;
	this->dt = std::numeric_limits<T>::infinity();
}

template <typename T>
T SimulationNBody<T>:: getDt()
{
	return this->dt;
}

template <typename T>
void SimulationNBody<T>::computeIteration()
{
	for(unsigned long iBody = 0; iBody < this->bodies.n * omp_get_max_threads(); iBody++)
	{
		this->accelerations.x[iBody] = 0.0;
		this->accelerations.y[iBody] = 0.0;
		this->accelerations.z[iBody] = 0.0;
	}

	for(unsigned long iBody = 0; iBody < this->bodies.n; iBody++)
		this->closestNeighborDist[iBody] = std::numeric_limits<T>::infinity();

	this->computeBodiesAccelerationV2();
	if(!this->dtConstant)
		this->findTimeStep();
	this->bodies.updatePositionsAndVelocities(this->accelerations, this->dt);
}

template <typename T>
void SimulationNBody<T>::computeBodiesAcceleration()
{
#pragma omp parallel for schedule(runtime)
	for(unsigned long iBody = 0; iBody < this->bodies.n; iBody++)
		for(unsigned long jBody = 0; jBody < this->bodies.n; jBody++)
			if(iBody != jBody)
				//this->computeAccelerationBetweenTwoBodiesNaive(iBody, jBody);
				this->computeAccelerationBetweenTwoBodies(iBody, jBody);
}

/* 
	AI  = (23 * blockSize * nBodies * nBlocks)  / ((4 * blockSize + 7 * nBodies) * nBlocks) <=>
	AI  = (23 * blockSize * nBodies)            /  (4 * blockSize + 7 * nBodies)            <=>
	AI  = (23 * blockSize * nBlock * blockSize) /  (4 * blockSize + 7 * nBlock * blockSize) <=>
	AI  = (23 * nBlock * blockSize²)            / ((4 + 7 * nBlock) * blockSize)            <=>
	AI  = (23 * nBlock * blockSize)             /  (4 + 7 * nBlock)                         <=>
	AI ~= (23 * nBlock * blockSize)             /      (7 * nBlock)                         <=>
	AI ~= (23 * blockSize)                      /       7
	-------------------------------------------------------------------------------------------
	OI  = AI                                    /      sizeof(T)                            <=>
	OI  = (23 * blockSize)                      / (7 * sizeof(T))
*/
template <typename T>
void SimulationNBody<T>::computeBodiesAccelerationCB()
{
	unsigned long blockSize = 512;
	// flops  = 23 * blockSize * nBodies      * nBlocks
	// memops = (4 * blockSize + 7 * nBodies) * nBlocks
	for(unsigned long jOff = 0; jOff < this->bodies.n; jOff += blockSize)
	{
		blockSize = std::min(blockSize, this->bodies.n - jOff);
		// flops  = 23 * blockSize * nBodies
		// memops =  4 * blockSize + 7 * nBodies 
#pragma omp parallel for schedule(runtime)
		for(unsigned long iBody = 0; iBody < this->bodies.n; iBody++)
			// flops  = 23 * blockSize
			// memops =  4 * blockSize + 7 
			for(unsigned long jBody = jOff; jBody < jOff + blockSize; jBody++)
				if(iBody != jBody)
					this->computeAccelerationBetweenTwoBodies(iBody, jBody);
					//this->computeAccelerationBetweenTwoBodiesNaive(iBody, jBody);
	}
}

// 23 flops
template <typename T>
void SimulationNBody<T>::computeAccelerationBetweenTwoBodiesNaive(const unsigned long iBody, const unsigned long jBody)
{
	assert(iBody != jBody);

	const T diffPosX = this->bodies.positions.x[jBody] - this->bodies.positions.x[iBody]; // 1 flop
	const T diffPosY = this->bodies.positions.y[jBody] - this->bodies.positions.y[iBody]; // 1 flop
	const T diffPosZ = this->bodies.positions.z[jBody] - this->bodies.positions.z[iBody]; // 1 flop

	// compute distance between iBody and jBody: Dij
	const T dist = std::sqrt((diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ)); // 6 flops

	// we cannot divide by 0
	if(dist == 0)
	{
		std::cout << "Collision at {" << this->bodies.positions.x[jBody] << ", "
		                              << this->bodies.positions.y[jBody] << ", "
		                              << this->bodies.positions.z[jBody] << "}" << std::endl;
		assert(dist != 0);
	}

	// compute the force value between iBody and jBody: || F || = G.mi.mj / Dij²
	const T force = G * this->bodies.masses[iBody] * this->bodies.masses[jBody] / (dist * dist); // 4 flops

	// compute the acceleration value: || a || = || F || / mi
	const T acc = force / this->bodies.masses[iBody]; // 1 flop

	// normalize and add acceleration value into acceleration vector: a += || a ||.u
	this->accelerations.x[iBody] += acc * (diffPosX / dist); // 3 flops
	this->accelerations.y[iBody] += acc * (diffPosY / dist); // 3 flops
	this->accelerations.z[iBody] += acc * (diffPosZ / dist); // 3 flops

	if(!this->dtConstant)
		if(dist < this->closestNeighborDist[iBody])
			this->closestNeighborDist[iBody] = dist;
}

// 18 flops
template <typename T>
void SimulationNBody<T>::computeAccelerationBetweenTwoBodies(const unsigned long iBody, const unsigned long jBody)
{
	assert(iBody != jBody);

	const T diffPosX = this->bodies.positions.x[jBody] - this->bodies.positions.x[iBody]; // 1 flop
	const T diffPosY = this->bodies.positions.y[jBody] - this->bodies.positions.y[iBody]; // 1 flop
	const T diffPosZ = this->bodies.positions.z[jBody] - this->bodies.positions.z[iBody]; // 1 flop
	const T squareDist = (diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ); // 5 flops
	const T dist = std::sqrt(squareDist); // 1 flop
	assert(dist != 0);

	const T acc = G * this->bodies.masses[jBody] / (squareDist * dist); // 3 flops
	this->accelerations.x[iBody] += acc * diffPosX; // 2 flop
	this->accelerations.y[iBody] += acc * diffPosY; // 2 flop
	this->accelerations.z[iBody] += acc * diffPosZ; // 2 flop

	if(!this->dtConstant)
		if(dist < this->closestNeighborDist[iBody])
			this->closestNeighborDist[iBody] = dist;
}

template <typename T>
void SimulationNBody<T>::computeBodiesAccelerationV2()
{
#pragma omp parallel for schedule(runtime)
	for(unsigned long iBody = 0; iBody < this->bodies.n; iBody++)
		for(unsigned long jBody = iBody +1; jBody < this->bodies.n; jBody++)
			//this->computeAccelerationBetweenTwoBodiesNaiveV2(iBody, jBody);
			this->computeAccelerationBetweenTwoBodiesV2(iBody, jBody);

	if(omp_get_max_threads() > 1)
		for(unsigned long iBody = 0; iBody < this->bodies.n; iBody++)
			for(unsigned short iThread = 1; iThread < omp_get_max_threads(); iThread++)
			{
				this->accelerations.x[iBody] += this->accelerations.x[iBody + iThread * this->bodies.n];
				this->accelerations.y[iBody] += this->accelerations.y[iBody + iThread * this->bodies.n];
				this->accelerations.z[iBody] += this->accelerations.z[iBody + iThread * this->bodies.n];
			}
}

// 33 flops
template <typename T>
void SimulationNBody<T>::computeAccelerationBetweenTwoBodiesNaiveV2(const unsigned long iBody, const unsigned long jBody)
{
	assert(iBody != jBody);

	const T diffPosX = this->bodies.positions.x[jBody] - this->bodies.positions.x[iBody]; // 1 flop
	const T diffPosY = this->bodies.positions.y[jBody] - this->bodies.positions.y[iBody]; // 1 flop
	const T diffPosZ = this->bodies.positions.z[jBody] - this->bodies.positions.z[iBody]; // 1 flop

	// compute distance between iBody and jBody: Dij
	const T dist = std::sqrt((diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ)); // 6 flops

	// we cannot divide by 0
	if(dist == 0)
	{
		std::cout << "Collision at {" << this->bodies.positions.x[jBody] << ", "
		                              << this->bodies.positions.y[jBody] << ", "
		                              << this->bodies.positions.z[jBody] << "}" << std::endl;
		assert(dist != 0);
	}

	// compute the force value between iBody and jBody: || F || = G.mi.mj / Dij²
	const T force = G * this->bodies.masses[iBody] * this->bodies.masses[jBody] / (dist * dist); // 4 flops

	// compute the acceleration value: || a || = || F || / mi
	T acc = force / this->masses[iBody]; // 1 flop

	// normalize and add acceleration value into acceleration vector: a += || a ||.u
	this->accelerations.x[iBody + omp_get_thread_num() * this->bodies.n] += acc * (diffPosX / dist); // 3 flops
	this->accelerations.y[iBody + omp_get_thread_num() * this->bodies.n] += acc * (diffPosY / dist); // 3 flops
	this->accelerations.z[iBody + omp_get_thread_num() * this->bodies.n] += acc * (diffPosZ / dist); // 3 flops

	// compute the acceleration value: || a || = || F || / mj
	acc = force / this->masses[jBody]; // 1 flop

	this->accelerations.x[jBody + omp_get_thread_num() * this->bodies.n] -= acc * (diffPosX / dist); // 3 flops
	this->accelerations.y[jBody + omp_get_thread_num() * this->bodies.n] -= acc * (diffPosY / dist); // 3 flops
	this->accelerations.z[jBody + omp_get_thread_num() * this->bodies.n] -= acc * (diffPosZ / dist); // 3 flops

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
void SimulationNBody<T>::computeAccelerationBetweenTwoBodiesV2(const unsigned long iBody, const unsigned long jBody)
{
	assert(iBody != jBody);

	const T diffPosX  = this->bodies.positions.x[jBody] - this->bodies.positions.x[iBody]; // 1 flop
	const T diffPosY  = this->bodies.positions.y[jBody] - this->bodies.positions.y[iBody]; // 1 flop
	const T diffPosZ  = this->bodies.positions.z[jBody] - this->bodies.positions.z[iBody]; // 1 flop

	const T squareDist = (diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ); // 5 flops
	const T dist = std::sqrt(squareDist); // 1 flops
	assert(dist != 0);

	const T force = G / (squareDist * dist); // 2 flops

	T acc = force * this->bodies.masses[jBody]; // 1 flop
	this->accelerations.x[iBody + omp_get_thread_num() * this->bodies.n] += acc * diffPosX; // 2 flops
	this->accelerations.y[iBody + omp_get_thread_num() * this->bodies.n] += acc * diffPosY; // 2 flops
	this->accelerations.z[iBody + omp_get_thread_num() * this->bodies.n] += acc * diffPosZ; // 2 flops

	acc = force * this->bodies.masses[iBody]; // 1 flop
	this->accelerations.x[jBody + omp_get_thread_num() * this->bodies.n] -= acc * diffPosX; // 2 flops
	this->accelerations.y[jBody + omp_get_thread_num() * this->bodies.n] -= acc * diffPosY; // 2 flops
	this->accelerations.z[jBody + omp_get_thread_num() * this->bodies.n] -= acc * diffPosZ; // 2 flops

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

template <typename T>
void SimulationNBody<T>::findTimeStep()
{
	if(!this->dtConstant)
	{
		this->dt = std::numeric_limits<T>::infinity();
		for(unsigned long iBody = 0; iBody < this->bodies.n; iBody++)
		{
			const T newDt = computeTimeStep(iBody);

			if(newDt < this->dt)
				this->dt = newDt;
		}
	}
}

template <typename T>
T SimulationNBody<T>::computeTimeStep(const unsigned long iBody)
{
	// || velocity[iBody] ||
	const T v = std::sqrt((this->bodies.velocities.x[iBody] * this->bodies.velocities.x[iBody]) +
	                      (this->bodies.velocities.y[iBody] * this->bodies.velocities.y[iBody]) +
	                      (this->bodies.velocities.z[iBody] * this->bodies.velocities.z[iBody]));

	// || acceleration[iBody] ||
	const T a = std::sqrt((this->accelerations.x[iBody] * this->accelerations.x[iBody]) +
	                      (this->accelerations.y[iBody] * this->accelerations.y[iBody]) +
	                      (this->accelerations.z[iBody] * this->accelerations.z[iBody]));

	/*
	 * compute dt
	 * solve:  (a/2)*dt^2 + v*dt + (-0.1)*ClosestNeighborDist = 0
	 * <=>     dt = [ (-v) +/-  sqrt( v^2 - 4 * (a/2) * (-0.1)*ClosestNeighborDist ) ] / [ 2 (a/2) ]
	 *
	 * dt should be positive (+/- becomes + because result of sqrt is positive)
	 * <=>     dt = [ -v + sqrt( v^2 + 0.2*ClosestNeighborDist*a) ] / a
	 */
	T dt = (std::sqrt(v * v + 0.2 * a * this->closestNeighborDist[iBody]) - v) / a;

	if(dt == 0)
		dt = std::numeric_limits<T>::epsilon() / a;

	return dt;
}
