/*!
 * \file    SimulationNBodyV2.hxx
 * \brief   Implementation of SimulationNBodyLocal (n²/2 computations).
 * \author  A. Cassagne
 * \date    2014
 *
 * \section LICENSE
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode).
 */
#include <cmath>
#include <limits>
#include <string>
#include <cassert>
#include <fstream>
#include <iostream>
#include <mipp.h>

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

#include "SimulationNBodyV2.h"

template <typename T>
SimulationNBodyV2<T>::SimulationNBodyV2(const unsigned long nBodies)
	: SimulationNBodyLocal<T>(nBodies)
{
	this->reAllocateBuffers();
}

template <typename T>
SimulationNBodyV2<T>::SimulationNBodyV2(const std::string inputFileName)
	: SimulationNBodyLocal<T>(inputFileName)
{
	this->reAllocateBuffers();
}

template <typename T>
SimulationNBodyV2<T>::~SimulationNBodyV2()
{
#ifdef __ARM_NEON__
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
#else
	if(this->accelerations.x != nullptr) {
		_mm_free(this->accelerations.x);
		this->accelerations.x = nullptr;
	}
	if(this->accelerations.y != nullptr) {
		_mm_free(this->accelerations.y);
		this->accelerations.y = nullptr;
	}
	if(this->accelerations.z != nullptr) {
		_mm_free(this->accelerations.z);
		this->accelerations.z = nullptr;
	}
#endif
}

template <typename T>
void SimulationNBodyV2<T>::reAllocateBuffers()
{
	if(this->nMaxThreads > 1)
	{
		// TODO: this is not optimal to deallocate and to reallocate data
#ifdef __ARM_NEON__
		if(this->accelerations.x != nullptr)
			delete[] this->accelerations.x;
		if(this->accelerations.y != nullptr)
			delete[] this->accelerations.y;
		if(this->accelerations.z != nullptr)
			delete[] this->accelerations.z;

		this->accelerations.x = new T[(this->bodies->getN() + this->bodies->getPadding()) * this->nMaxThreads];
		this->accelerations.y = new T[(this->bodies->getN() + this->bodies->getPadding()) * this->nMaxThreads];
		this->accelerations.z = new T[(this->bodies->getN() + this->bodies->getPadding()) * this->nMaxThreads];
#else
		if(this->accelerations.x != nullptr)
			_mm_free(this->accelerations.x);
		if(this->accelerations.y != nullptr)
			_mm_free(this->accelerations.y);
		if(this->accelerations.z != nullptr)
			_mm_free(this->accelerations.z);

		this->accelerations.x = (T*)_mm_malloc((this->bodies->getN() + this->bodies->getPadding()) *
		                                        this->nMaxThreads * sizeof(T), mipp::RequiredAlignment);
		this->accelerations.y = (T*)_mm_malloc((this->bodies->getN() + this->bodies->getPadding()) *
		                                        this->nMaxThreads * sizeof(T), mipp::RequiredAlignment);
		this->accelerations.z = (T*)_mm_malloc((this->bodies->getN() + this->bodies->getPadding()) *
		                                        this->nMaxThreads * sizeof(T), mipp::RequiredAlignment);
#endif
		this->allocatedBytes += (this->bodies->getN() + this->bodies->getPadding()) *
		                        sizeof(T) * (this->nMaxThreads - 1) * 3;
	}

	this->flopsPerIte = 25.f * ((float)this->bodies->getN() * 0.5f) * (float)this->bodies->getN();
}

template <typename T>
void SimulationNBodyV2<T>::initIteration()
{
	for(unsigned long iBody = 0; iBody < this->bodies->getN() * this->nMaxThreads; iBody++)
	{
		this->accelerations.x[iBody] = 0.0;
		this->accelerations.y[iBody] = 0.0;
		this->accelerations.z[iBody] = 0.0;
	}

	for(unsigned long iBody = 0; iBody < this->bodies->getN(); iBody++)
		this->closestNeighborDist[iBody] = std::numeric_limits<T>::infinity();
}

template <typename T>
void SimulationNBodyV2<T>::computeLocalBodiesAcceleration()
{
	const T *masses     = this->bodies->getMasses();
	const T *positionsX = this->bodies->getPositionsX();
	const T *positionsY = this->bodies->getPositionsY();
	const T *positionsZ = this->bodies->getPositionsZ();

#pragma omp parallel
{
	const unsigned int  tid     = omp_get_thread_num();
	const unsigned long tStride = tid * this->bodies->getN();

#pragma omp for schedule(runtime)
	for(unsigned long iBody = 0; iBody < this->bodies->getN(); iBody++)
		for(unsigned long jBody = iBody +1; jBody < this->bodies->getN(); jBody++)
			SimulationNBodyV2<T>::computeAccelerationBetweenTwoBodies(this->G,
			                                                          masses                   [iBody          ],
			                                                          positionsX               [iBody          ],
			                                                          positionsY               [iBody          ],
			                                                          positionsZ               [iBody          ],
			                                                          this->accelerations.x    [iBody + tStride],
			                                                          this->accelerations.y    [iBody + tStride],
			                                                          this->accelerations.z    [iBody + tStride],
			                                                          this->closestNeighborDist[iBody          ],
			                                                          masses                   [jBody          ],
			                                                          positionsX               [jBody          ],
			                                                          positionsY               [jBody          ],
			                                                          positionsZ               [jBody          ],
			                                                          this->accelerations.x    [jBody + tStride],
			                                                          this->accelerations.y    [jBody + tStride],
			                                                          this->accelerations.z    [jBody + tStride],
			                                                          this->closestNeighborDist[jBody          ]);
}

	if(this->nMaxThreads > 1)
		for(unsigned long iBody = 0; iBody < this->bodies->getN(); iBody++)
			for(unsigned iThread = 1; iThread < this->nMaxThreads; iThread++)
			{
				this->accelerations.x[iBody] += this->accelerations.x[iBody + iThread * this->bodies->getN()];
				this->accelerations.y[iBody] += this->accelerations.y[iBody + iThread * this->bodies->getN()];
				this->accelerations.z[iBody] += this->accelerations.z[iBody + iThread * this->bodies->getN()];
			}
}

// 25 flops
template <typename T>
void SimulationNBodyV2<T>::computeAccelerationBetweenTwoBodies(const T &G,
                                                               const T &mi,
                                                               const T &qiX, const T &qiY, const T &qiZ,
                                                                     T &aiX,       T &aiY,       T &aiZ,
                                                                     T &closNeighi,
                                                               const T &mj,
                                                               const T &qjX, const T &qjY, const T &qjZ,
                                                                     T &ajX,       T &ajY,       T &ajZ,
                                                                     T &closNeighj)
{
	const T rijX = qjX - qiX; // 1 flop
	const T rijY = qjY - qiY; // 1 flop
	const T rijZ = qjZ - qiZ; // 1 flop

	// compute the distance || rij ||² between body i and body j
	const T rijSquared = (rijX * rijX) + (rijY * rijY) + (rijZ * rijZ); // 5 flops

	// compute the distance || rij ||
	const T rij = std::sqrt(rijSquared); // 1 flop

	// we cannot divide by 0
	assert(rij != 0);

	// compute the acceleration value between body i and body j: || ai || = G.mj / || rij ||^3
	const T aTmp = G / (rijSquared * rij); // 2 flops
	const T ai = aTmp * mj; // 1 flops

	// add the acceleration value into the acceleration vector: ai += || ai ||.rij
	aiX += ai * rijX; // 2 flops
	aiY += ai * rijY; // 2 flops
	aiZ += ai * rijZ; // 2 flops

	// compute the acceleration value between body j and body i: || aj || = G.mi / || rij ||^3
	const T aj = aTmp * mi; // 1 flops

	// add the acceleration value into the acceleration vector: aj -= || aj ||.rij
	ajX -= aj * rijX; // 2 flops
	ajY -= aj * rijY; // 2 flops
	ajZ -= aj * rijZ; // 2 flops

	closNeighi = std::min(closNeighi, rij);
	if(rij < closNeighj)
#pragma omp critical
		closNeighj = std::min(closNeighj, rij);
}
