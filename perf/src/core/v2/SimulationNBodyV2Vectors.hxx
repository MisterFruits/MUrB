/*
 * Do not remove.
 * Gabriel Hautreux, CINES, gabrielhautreux@gmail.com
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

#include "../../utils/myIntrinsicsPlusPlus.h"

#include "SimulationNBodyV2Vectors.h"

template <typename T>
SimulationNBodyV2Vectors<T>::SimulationNBodyV2Vectors(const unsigned long nBodies)
	: SimulationNBody<T>(nBodies)
{
	this->reAllocateBuffers();
}

template <typename T>
SimulationNBodyV2Vectors<T>::SimulationNBodyV2Vectors(const std::string inputFileName)
	: SimulationNBody<T>(inputFileName)
{
	this->reAllocateBuffers();
}

template <typename T>
SimulationNBodyV2Vectors<T>::~SimulationNBodyV2Vectors()
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
void SimulationNBodyV2Vectors<T>::reAllocateBuffers()
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

		this->accelerations.x = new T[(this->bodies.getN() + this->bodies.getPadding()) * this->nMaxThreads];
		this->accelerations.y = new T[(this->bodies.getN() + this->bodies.getPadding()) * this->nMaxThreads];
		this->accelerations.z = new T[(this->bodies.getN() + this->bodies.getPadding()) * this->nMaxThreads];
#else
		if(this->accelerations.x != nullptr)
			_mm_free(this->accelerations.x);
		if(this->accelerations.y != nullptr)
			_mm_free(this->accelerations.y);
		if(this->accelerations.z != nullptr)
			_mm_free(this->accelerations.z);

		this->accelerations.x = (T*)_mm_malloc((this->bodies.getN() + this->bodies.getPadding()) *
		                                        this->nMaxThreads * sizeof(T), mipp::RequiredAlignement);
		this->accelerations.y = (T*)_mm_malloc((this->bodies.getN() + this->bodies.getPadding()) *
		                                        this->nMaxThreads * sizeof(T), mipp::RequiredAlignement);
		this->accelerations.z = (T*)_mm_malloc((this->bodies.getN() + this->bodies.getPadding()) *
		                                        this->nMaxThreads * sizeof(T), mipp::RequiredAlignement);
#endif
		this->allocatedBytes += (this->bodies.getN() + this->bodies.getPadding()) *
		                        sizeof(T) * (this->nMaxThreads - 1) * 3;
	}

	this->flopsPerIte = 25 * (this->bodies.getN() * 0.5) * this->bodies.getN();
}

template <typename T>
void SimulationNBodyV2Vectors<T>::initIteration()
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
void SimulationNBodyV2Vectors<T>::computeBodiesAcceleration()
{
	const T *masses = this->getBodies().getMasses();

	const T *positionsX = this->getBodies().getPositionsX();
	const T *positionsY = this->getBodies().getPositionsY();
	const T *positionsZ = this->getBodies().getPositionsZ();

#pragma omp parallel
{
	const unsigned long thStride = omp_get_thread_num() * (this->bodies.getN() + this->bodies.getPadding());

#pragma omp for schedule(runtime)
	for(unsigned long iVec = 0; iVec < this->bodies.getNVecs(); iVec++)
	{
		const unsigned long iVecOff = iVec * mipp::vectorSize<T>();
  
		// computation of the vector number iVec with itself
		for(unsigned short iVecPos = 0; iVecPos < mipp::vectorSize<T>(); iVecPos++)
		{
			const unsigned long iBody = iVecPos + iVecOff;
			for(unsigned short jVecPos = iVecPos +1; jVecPos < mipp::vectorSize<T>(); jVecPos++)
			{
				const unsigned long jBody = jVecPos + iVecOff;
					this->computeAccelerationBetweenTwoBodies(positionsX               [iBody           ],
					                                          positionsY               [iBody           ],
					                                          positionsZ               [iBody           ],
					                                          this->accelerations.x    [iBody + thStride],
					                                          this->accelerations.y    [iBody + thStride],
					                                          this->accelerations.z    [iBody + thStride],
					                                          this->closestNeighborDist[iBody           ],
					                                          masses                   [iBody           ],
					                                          positionsX               [jBody           ],
					                                          positionsY               [jBody           ],
					                                          positionsZ               [jBody           ],
					                                          this->accelerations.x    [jBody + thStride],
					                                          this->accelerations.y    [jBody + thStride],
					                                          this->accelerations.z    [jBody + thStride],
					                                          this->closestNeighborDist[jBody           ],
					                                          masses                   [jBody           ]);
			}
		}

		// computation of the vector number iVec with the following other vectors
		for(unsigned long jVec = iVec +1; jVec < this->bodies.getNVecs(); jVec++){
					const unsigned long jo = jVec * mipp::vectorSize<T>();
			for(unsigned short iVecPos = 0; iVecPos < mipp::vectorSize<T>(); iVecPos++)
			{
				const unsigned long iBody = iVecPos + iVecOff;
				for(unsigned short jVecPos = 0; jVecPos < mipp::vectorSize<T>(); jVecPos++)
				{
					const unsigned long jBody = jVecPos + jVec * mipp::vectorSize<T>();
					this->computeAccelerationBetweenTwoBodies(positionsX               [iBody           ],
					                                          positionsY               [iBody           ],
					                                          positionsZ               [iBody           ],
					                                          this->accelerations.x    [iBody + thStride],
					                                          this->accelerations.y    [iBody + thStride],
					                                          this->accelerations.z    [iBody + thStride],
					                                          this->closestNeighborDist[iBody           ],
					                                          masses                   [iBody           ],
					                                          positionsX               [jBody           ],
					                                          positionsY               [jBody           ],
					                                          positionsZ               [jBody           ],
					                                          this->accelerations.x    [jBody + thStride],
					                                          this->accelerations.y    [jBody + thStride],
					                                          this->accelerations.z    [jBody + thStride],
					                                          this->closestNeighborDist[jBody           ],
					                                          masses                   [jBody           ]);
				}
			}
	}
	}
}

	if(this->nMaxThreads > 1)
		for(unsigned long iBody = 0; iBody < this->bodies.getN(); iBody++)
			for(unsigned iThread = 1; iThread < this->nMaxThreads; iThread++)
			{
				this->accelerations.x[iBody] +=
				            this->accelerations.x[iBody + iThread * (this->bodies.getN() + this->bodies.getPadding())];
				this->accelerations.y[iBody] +=
				            this->accelerations.y[iBody + iThread * (this->bodies.getN() + this->bodies.getPadding())];
				this->accelerations.z[iBody] +=
				            this->accelerations.z[iBody + iThread * (this->bodies.getN() + this->bodies.getPadding())];
			}
}

// 25 flops
template <typename T>
void SimulationNBodyV2Vectors<T>::computeAccelerationBetweenTwoBodies(const T &iPosX, const T &iPosY, const T &iPosZ,
                                                                            T &iAccsX,      T &iAccsY,      T &iAccsZ,
                                                                            T &iClosNeiDist,
                                                                      const T &iMasses,
                                                                      const T &jPosX, const T &jPosY, const T &jPosZ,
                                                                            T &jAccsX,      T &jAccsY,      T &jAccsZ,
                                                                            T &jClosNeiDist,
                                                                      const T &jMasses)
{
	const T diffPosX = jPosX - iPosX; // 1 flop
	const T diffPosY = jPosY - iPosY; // 1 flop
	const T diffPosZ = jPosZ - iPosZ; // 1 flop
	const T squareDist = (diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ); // 5 flops
	const T dist = std::sqrt(squareDist); // 1 flop

	const T force = this->G / (squareDist * dist); // 2 flops

	T acc = force * jMasses; // 1 flop

	iAccsX += acc * diffPosX; // 2 flop
	iAccsY += acc * diffPosY; // 2 flop
	iAccsZ += acc * diffPosZ; // 2 flop

	acc = force * iMasses; // 1 flop

	jAccsX -= acc * diffPosX; // 2 flop
	jAccsY -= acc * diffPosY; // 2 flop
	jAccsZ -= acc * diffPosZ; // 2 flop
	
	if(!this->dtConstant)
	{
		iClosNeiDist = std::min(iClosNeiDist, dist);
		if(dist < jClosNeiDist)
#pragma omp critical
			jClosNeiDist = std::min(jClosNeiDist, dist);
	}
}
