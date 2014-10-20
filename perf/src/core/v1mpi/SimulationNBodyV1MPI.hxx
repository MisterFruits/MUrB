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

#include "SimulationNBodyV1MPI.h"

template <typename T>
SimulationNBodyV1MPI<T>::SimulationNBodyV1MPI(const unsigned long nBodies)
	: SimulationNBodyMPI<T>(nBodies)
{
	this->init();
}

template <typename T>
SimulationNBodyV1MPI<T>::SimulationNBodyV1MPI(const std::string inputFileName)
	: SimulationNBodyMPI<T>(inputFileName)
{
	this->init();
}

template <typename T>
void SimulationNBodyV1MPI<T>::init()
{
	this->flopsPerIte = 18 * ((this->bodies.getN() * this->MPISize) -1) * (this->bodies.getN() * this->MPISize);
}

template <typename T>
SimulationNBodyV1MPI<T>::~SimulationNBodyV1MPI()
{
}

template <typename T>
void SimulationNBodyV1MPI<T>::initIteration()
{
	for(unsigned long iBody = 0; iBody < this->bodies.getN(); iBody++)
	{
		this->accelerations.x[iBody] = 0.0;
		this->accelerations.y[iBody] = 0.0;
		this->accelerations.z[iBody] = 0.0;

		this->closestNeighborDist[iBody] = std::numeric_limits<T>::infinity();
	}
}

template <typename T>
void SimulationNBodyV1MPI<T>::computeLocalBodiesAcceleration()
{
	const T *masses     = this->bodies.getMasses();
	const T *positionsX = this->bodies.getPositionsX();
	const T *positionsY = this->bodies.getPositionsY();
	const T *positionsZ = this->bodies.getPositionsZ();

#pragma omp parallel for schedule(runtime)
	for(unsigned long iBody = 0; iBody < this->bodies.getN(); iBody++)
		for(unsigned long jBody = 0; jBody < this->bodies.getN(); jBody++)
			if(iBody != jBody)
				this->computeAccelerationBetweenTwoBodies(positionsX               [iBody],
				                                          positionsY               [iBody],
				                                          positionsZ               [iBody],
				                                          this->accelerations.x    [iBody],
				                                          this->accelerations.y    [iBody],
				                                          this->accelerations.z    [iBody],
				                                          this->closestNeighborDist[iBody],
				                                          masses                   [jBody],
				                                          positionsX               [jBody],
				                                          positionsY               [jBody],
				                                          positionsZ               [jBody]);
}

template <typename T>
void SimulationNBodyV1MPI<T>::computeNeighborBodiesAcceleration()
{
	//const T *masses     = this->bodies.getMasses();
	const T *positionsX = this->bodies.getPositionsX();
	const T *positionsY = this->bodies.getPositionsY();
	const T *positionsZ = this->bodies.getPositionsZ();

	const T *neighMasses     = this->neighborBodies->getMasses();
	const T *neighPositionsX = this->neighborBodies->getPositionsX();
	const T *neighPositionsY = this->neighborBodies->getPositionsY();
	const T *neighPositionsZ = this->neighborBodies->getPositionsZ();

#pragma omp parallel for schedule(runtime)
	for(unsigned long iBody = 0; iBody < this->bodies.getN(); iBody++)
		for(unsigned long jBody = 0; jBody < this->bodies.getN(); jBody++)
			this->computeAccelerationBetweenTwoBodies(positionsX               [iBody],
			                                          positionsY               [iBody],
			                                          positionsZ               [iBody],
			                                          this->accelerations.x    [iBody],
			                                          this->accelerations.y    [iBody],
			                                          this->accelerations.z    [iBody],
			                                          this->closestNeighborDist[iBody],
			                                          neighMasses              [jBody],
			                                          neighPositionsX          [jBody],
			                                          neighPositionsY          [jBody],
			                                          neighPositionsZ          [jBody]);
}

// 23 flops
template <typename T>
void SimulationNBodyV1MPI<T>::computeAccelerationBetweenTwoBodiesNaive(const T &iMasses,
                                                                       const T &iPosX, const T &iPosY, const T &iPosZ,
                                                                             T &iAccsX,      T &iAccsY,      T &iAccsZ,
                                                                             T &iClosNeiDist,
                                                                       const T &jMasses,
                                                                       const T &jPosX, const T &jPosY, const T &jPosZ)
{
	assert(iBody != jBody);

	const T diffPosX = jPosX - iPosX; // 1 flop
	const T diffPosY = jPosY - iPosY; // 1 flop
	const T diffPosZ = jPosZ - iPosZ; // 1 flop

	// compute distance between iBody and jBody: Dij
	const T dist = std::sqrt((diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ)); // 6 flops

	// we cannot divide by 0
	if(dist == 0)
	{
		std::cout << "Collision at {" << jPosX << ", "
		                              << jPosY << ", "
		                              << jPosZ << "}" << std::endl;
		assert(dist != 0);
	}

	// compute the force value between iBody and jBody: || F || = G.mi.mj / Dij²
	const T force = (this->G * iMasses * jMasses / (dist * dist)); // 4 flops

	// compute the acceleration value: || a || = || F || / mi
	const T acc = force / iMasses; // 1 flop

	// normalize and add acceleration value into acceleration vector: a += || a ||.u
	iAccsX += acc * (diffPosX / dist); // 3 flops
	iAccsY += acc * (diffPosY / dist); // 3 flops
	iAccsZ += acc * (diffPosZ / dist); // 3 flops

	if(!this->dtConstant)
		if(dist < iClosNeiDist)
			iClosNeiDist = dist;
}

// 18 flops
template <typename T>
void SimulationNBodyV1MPI<T>::computeAccelerationBetweenTwoBodies(const T &iPosX, const T &iPosY, const T &iPosZ,
                                                                        T &iAccsX,      T &iAccsY,      T &iAccsZ,
                                                                        T &iClosNeiDist,
                                                                  const T &jMasses,
                                                                  const T &jPosX, const T &jPosY, const T &jPosZ)
{
	assert(iBody != jBody);

	const T diffPosX = jPosX - iPosX; // 1 flop
	const T diffPosY = jPosY - iPosY; // 1 flop
	const T diffPosZ = jPosZ - iPosZ; // 1 flop
	const T squareDist = (diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ); // 5 flops
	const T dist = std::sqrt(squareDist); // 1 flop
	assert(dist != 0);

	const T acc = this->G * jMasses / (squareDist * dist); // 3 flops
	iAccsX += acc * diffPosX; // 2 flop
	iAccsY += acc * diffPosY; // 2 flop
	iAccsZ += acc * diffPosZ; // 2 flop

	if(!this->dtConstant)
		if(dist < iClosNeiDist)
			iClosNeiDist = dist;
}
