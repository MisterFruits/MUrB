/*!
 * \file    SimulationNBodyV2Vectors.hxx
 * \brief   Implementation of SimulationNBodyLocal with vector size stride loops (n²/2 computations).
 * \author  G. Hautreux
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

#include "../../../../../common/utils/mipp.h"

#include "SimulationNBodyV2Vectors.h"

template <typename T>
SimulationNBodyV2Vectors<T>::SimulationNBodyV2Vectors(const unsigned long nBodies)
	: SimulationNBodyV2<T>(nBodies)
{
	this->reAllocateBuffers();
}

template <typename T>
SimulationNBodyV2Vectors<T>::SimulationNBodyV2Vectors(const std::string inputFileName)
	: SimulationNBodyV2<T>(inputFileName)
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
		                                        this->nMaxThreads * sizeof(T), mipp::RequiredAlignement);
		this->accelerations.y = (T*)_mm_malloc((this->bodies->getN() + this->bodies->getPadding()) *
		                                        this->nMaxThreads * sizeof(T), mipp::RequiredAlignement);
		this->accelerations.z = (T*)_mm_malloc((this->bodies->getN() + this->bodies->getPadding()) *
		                                        this->nMaxThreads * sizeof(T), mipp::RequiredAlignement);
#endif
		this->allocatedBytes += (this->bodies->getN() + this->bodies->getPadding()) *
		                        sizeof(T) * (this->nMaxThreads - 1) * 3;
	}

	this->flopsPerIte = 25 * (this->bodies->getN() * 0.5) * this->bodies->getN();
}

template <typename T>
void SimulationNBodyV2Vectors<T>::initIteration()
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
void SimulationNBodyV2Vectors<T>::computeLocalBodiesAcceleration()
{
	const T *masses = this->bodies->getMasses();

	const T *positionsX = this->bodies->getPositionsX();
	const T *positionsY = this->bodies->getPositionsY();
	const T *positionsZ = this->bodies->getPositionsZ();

#pragma omp parallel
{
	const unsigned int  tid     = omp_get_thread_num();
	const unsigned long tStride = tid * (this->bodies->getN() + this->bodies->getPadding());

#pragma omp for schedule(runtime)
	for(unsigned long iVec = 0; iVec < this->bodies->getNVecs(); iVec++)
	{
		const unsigned long iVecOff = iVec * mipp::vectorSize<T>();
  
		// computation of the vector number iVec with itself
		for(unsigned short iVecPos = 0; iVecPos < mipp::vectorSize<T>(); iVecPos++)
		{
			const unsigned long iBody = iVecPos + iVecOff;
			for(unsigned short jVecPos = iVecPos +1; jVecPos < mipp::vectorSize<T>(); jVecPos++)
			{
				const unsigned long jBody = jVecPos + iVecOff;
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
		}

		// computation of the vector number iVec with the following other vectors
		for(unsigned long jVec = iVec +1; jVec < this->bodies->getNVecs(); jVec++)
			for(unsigned short iVecPos = 0; iVecPos < mipp::vectorSize<T>(); iVecPos++)
			{
				const unsigned long iBody = iVecPos + iVecOff;
				for(unsigned short jVecPos = 0; jVecPos < mipp::vectorSize<T>(); jVecPos++)
				{
					const unsigned long jBody = jVecPos + jVec * mipp::vectorSize<T>();
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
			}
	}
}

	if(this->nMaxThreads > 1)
		for(unsigned long iBody = 0; iBody < this->bodies->getN(); iBody++)
			for(unsigned iThread = 1; iThread < this->nMaxThreads; iThread++)
			{
				this->accelerations.x[iBody] +=
				          this->accelerations.x[iBody + iThread * (this->bodies->getN() + this->bodies->getPadding())];
				this->accelerations.y[iBody] +=
				          this->accelerations.y[iBody + iThread * (this->bodies->getN() + this->bodies->getPadding())];
				this->accelerations.z[iBody] +=
				          this->accelerations.z[iBody + iThread * (this->bodies->getN() + this->bodies->getPadding())];
			}
}
