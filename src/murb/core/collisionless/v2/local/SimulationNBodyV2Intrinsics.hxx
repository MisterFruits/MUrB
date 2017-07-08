/*!
 * \file    SimulationNBodyV2Intrinsics.hxx
 * \brief   Implementation of SimulationNBodyLocal with intrinsic function calls (n²/2 computations).
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

#include "SimulationNBodyV2Intrinsics.h"

template <typename T>
SimulationNBodyV2Intrinsics<T>::SimulationNBodyV2Intrinsics(const unsigned long nBodies)
	: SimulationNBodyV2<T>(nBodies)
{
	this->reAllocateBuffers();
}

template <typename T>
SimulationNBodyV2Intrinsics<T>::SimulationNBodyV2Intrinsics(const std::string inputFileName)
	: SimulationNBodyV2<T>(inputFileName)
{
	this->reAllocateBuffers();
}

template <typename T>
SimulationNBodyV2Intrinsics<T>::~SimulationNBodyV2Intrinsics()
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
void SimulationNBodyV2Intrinsics<T>::_reAllocateBuffers()
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
}


template <typename T>
void SimulationNBodyV2Intrinsics<T>::reAllocateBuffers()
{
	this->_reAllocateBuffers();
	this->flopsPerIte = 26.f * ((float)this->bodies->getN() * 0.5f) * (float)this->bodies->getN();
}

template <>
void SimulationNBodyV2Intrinsics<float>::reAllocateBuffers()
{
	this->_reAllocateBuffers();
	this->flopsPerIte = 27.f * ((float)this->bodies->getN() * 0.5f) * (float)this->bodies->getN();
}

template <typename T>
void SimulationNBodyV2Intrinsics<T>::_initIteration()
{
	for(unsigned long iBody = 0; iBody < this->bodies->getN() * this->nMaxThreads; iBody++)
	{
		this->accelerations.x[iBody] = 0.0;
		this->accelerations.y[iBody] = 0.0;
		this->accelerations.z[iBody] = 0.0;
	}
}

template <typename T>
void SimulationNBodyV2Intrinsics<T>::initIteration()
{
	this->_initIteration();

	for(unsigned long iBody = 0; iBody < this->bodies->getN(); iBody++)
		this->closestNeighborDist[iBody] = std::numeric_limits<T>::infinity();
}

template <>
void SimulationNBodyV2Intrinsics<float>::initIteration()
{
	this->_initIteration();

	for(unsigned long iBody = 0; iBody < this->bodies->getN(); iBody++)
		this->closestNeighborDist[iBody] = 0;
}

template <typename T>
void SimulationNBodyV2Intrinsics<T>::_computeLocalBodiesAcceleration()
{
	const T *masses = this->bodies->getMasses();

	const T *positionsX = this->bodies->getPositionsX();
	const T *positionsY = this->bodies->getPositionsY();
	const T *positionsZ = this->bodies->getPositionsZ();

	const mipp::reg rG = mipp::set1<T>(this->G);

#pragma omp parallel firstprivate(rG)
{
	const unsigned int  tid     = omp_get_thread_num();
	const unsigned long tStride = tid * (this->bodies->getN() + this->bodies->getPadding());

#pragma omp for schedule(runtime)
	for(unsigned long iVec = 0; iVec < this->bodies->getNVecs(); iVec++)
	{
		// computation of the vector number iVec with itself
		for(unsigned short iVecPos = 0; iVecPos < mipp::N<T>(); iVecPos++)
		{
			const unsigned long iBody = iVecPos + iVec * mipp::N<T>();
			for(unsigned short jVecPos = iVecPos +1; jVecPos < mipp::N<T>(); jVecPos++)
			{
				const unsigned long jBody = jVecPos + iVec * mipp::N<T>();
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
		
		// load vectors
		const mipp::reg rmi  = mipp::load<T>(masses     + iVec * mipp::N<T>());
		const mipp::reg rqiX = mipp::load<T>(positionsX + iVec * mipp::N<T>());
		const mipp::reg rqiY = mipp::load<T>(positionsY + iVec * mipp::N<T>());
		const mipp::reg rqiZ = mipp::load<T>(positionsZ + iVec * mipp::N<T>());

		mipp::reg raiX = mipp::load<T>(this->accelerations.x + iVec * mipp::N<T>() + tStride);
		mipp::reg raiY = mipp::load<T>(this->accelerations.y + iVec * mipp::N<T>() + tStride);
		mipp::reg raiZ = mipp::load<T>(this->accelerations.z + iVec * mipp::N<T>() + tStride);

		mipp::reg rclosNeighi;
		if(!this->dtConstant)
			rclosNeighi = mipp::load<T>(this->closestNeighborDist + iVec * mipp::N<T>());
		
		// computation of the vector number iVec with the following other vectors
		for(unsigned long jVec = iVec +1; jVec < this->bodies->getNVecs(); jVec++)
		{
			mipp::reg rmj  = mipp::load<T>(masses     + jVec * mipp::N<T>());
			mipp::reg rqjX = mipp::load<T>(positionsX + jVec * mipp::N<T>());
			mipp::reg rqjY = mipp::load<T>(positionsY + jVec * mipp::N<T>());
			mipp::reg rqjZ = mipp::load<T>(positionsZ + jVec * mipp::N<T>());

			mipp::reg rajX = mipp::load<T>(this->accelerations.x + jVec * mipp::N<T>() + tStride);
			mipp::reg rajY = mipp::load<T>(this->accelerations.y + jVec * mipp::N<T>() + tStride);
			mipp::reg rajZ = mipp::load<T>(this->accelerations.z + jVec * mipp::N<T>() + tStride);

			mipp::reg rclosNeighj;
			if(!this->dtConstant)
				rclosNeighj = mipp::load<T>(this->closestNeighborDist + jVec * mipp::N<T>());
	
			for(unsigned short iRot = 0; iRot < mipp::N<T>(); iRot++)
			{
				SimulationNBodyV2Intrinsics<T>::computeAccelerationBetweenTwoBodies(rG,
				                                                                    rmi,
				                                                                    rqiX, rqiY, rqiZ,
				                                                                    raiX, raiY, raiZ,
				                                                                    rclosNeighi,
				                                                                    rmj,
				                                                                    rqjX, rqjY, rqjZ,
				                                                                    rajX, rajY, rajZ,
				                                                                    rclosNeighj);

				rmj  = mipp::rrot<T>(rmj);
				rqjX = mipp::rrot<T>(rqjX); rqjY = mipp::rrot<T>(rqjY); rqjZ = mipp::rrot<T>(rqjZ);
				rajX = mipp::rrot<T>(rajX); rajY = mipp::rrot<T>(rajY); rajZ = mipp::rrot<T>(rajZ);
				rclosNeighj = mipp::rrot<T>(rclosNeighj);
			}
			
			mipp::store<T>(this->accelerations.x + jVec * mipp::N<T>() + tStride, rajX);
			mipp::store<T>(this->accelerations.y + jVec * mipp::N<T>() + tStride, rajY);
			mipp::store<T>(this->accelerations.z + jVec * mipp::N<T>() + tStride, rajZ);
		
			if(!this->dtConstant)
#pragma omp critical //TODO: this critical section is bad for performances when we use variable time step
				mipp::store<T>(this->closestNeighborDist + jVec * mipp::N<T>(), rclosNeighj);
		}

		mipp::store<T>(this->accelerations.x + iVec * mipp::N<T>() + tStride, raiX);
		mipp::store<T>(this->accelerations.y + iVec * mipp::N<T>() + tStride, raiY);
		mipp::store<T>(this->accelerations.z + iVec * mipp::N<T>() + tStride, raiZ);

		if(!this->dtConstant)
			mipp::store<T>(this->closestNeighborDist + iVec * mipp::N<T>(), rclosNeighi);
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


template <typename T>
void SimulationNBodyV2Intrinsics<T>::computeLocalBodiesAcceleration()
{
	this->_computeLocalBodiesAcceleration();
}

template <>
void SimulationNBodyV2Intrinsics<float>::computeLocalBodiesAcceleration()
{
	this->_computeLocalBodiesAcceleration();

	for(unsigned long iBody = 0; iBody < this->bodies->getN(); iBody++)
		this->closestNeighborDist[iBody] = 1.0 / this->closestNeighborDist[iBody];
}

// 26 flops
template <typename T>
void SimulationNBodyV2Intrinsics<T>::computeAccelerationBetweenTwoBodies(const mipp::reg &rG,
                                                                         const mipp::reg &rmi,
                                                                         const mipp::reg &rqiX,
                                                                         const mipp::reg &rqiY,
                                                                         const mipp::reg &rqiZ,
                                                                               mipp::reg &raiX,
                                                                               mipp::reg &raiY,
                                                                               mipp::reg &raiZ,
                                                                               mipp::reg &rclosNeighi,
                                                                         const mipp::reg &rmj,
                                                                         const mipp::reg &rqjX,
                                                                         const mipp::reg &rqjY,
                                                                         const mipp::reg &rqjZ,
                                                                               mipp::reg &rajX,
                                                                               mipp::reg &rajY,
                                                                               mipp::reg &rajZ,
                                                                               mipp::reg &rclosNeighj)
{
	mipp::reg rrijX = mipp::sub<T>(rqjX, rqiX);
	mipp::reg rrijY = mipp::sub<T>(rqjY, rqiY);
	mipp::reg rrijZ = mipp::sub<T>(rqjZ, rqiZ);

	mipp::reg rrijSquared = mipp::set1<T>((T)0);	
	rrijSquared = mipp::fmadd<T>(rrijX, rrijX, rrijSquared); // 2 flops
	rrijSquared = mipp::fmadd<T>(rrijY, rrijY, rrijSquared); // 2 flops
	rrijSquared = mipp::fmadd<T>(rrijZ, rrijZ, rrijSquared); // 2 flops

	mipp::reg rrij = mipp::sqrt<T>(rrijSquared); // 1 flop

	mipp::reg raTmp = mipp::div<T>(rG, mipp::mul<T>(rrij, rrijSquared)); // 2 flops

	mipp::reg rai = mipp::mul<T>(raTmp, rmj); // 1 flop

	raiX = mipp::fmadd<T>(rai, rrijX, raiX); // 2 flops
	raiY = mipp::fmadd<T>(rai, rrijY, raiY); // 2 flops
	raiZ = mipp::fmadd<T>(rai, rrijZ, raiZ); // 2 flops

	mipp::reg raj = mipp::mul<T>(raTmp, rmi); // 1 flop

	rajX = mipp::fnmadd<T>(raj, rrijX, rajX); // 2 flops
	rajY = mipp::fnmadd<T>(raj, rrijY, rajY); // 2 flops
	rajZ = mipp::fnmadd<T>(raj, rrijZ, rajZ); // 2 flops

	rclosNeighi = mipp::min<T>(rrij, rclosNeighi);
	rclosNeighj = mipp::min<T>(rrij, rclosNeighj); // TODO: this second minimum is not threads safe !!!
}

// 27 flops
template <>
void SimulationNBodyV2Intrinsics<float>::computeAccelerationBetweenTwoBodies(const mipp::reg &rG,
                                                                             const mipp::reg &rmi,
                                                                             const mipp::reg &rqiX,
                                                                             const mipp::reg &rqiY,
                                                                             const mipp::reg &rqiZ,
                                                                                   mipp::reg &raiX,
                                                                                   mipp::reg &raiY,
                                                                                   mipp::reg &raiZ,
                                                                                   mipp::reg &rclosNeighi,
                                                                             const mipp::reg &rmj,
                                                                             const mipp::reg &rqjX,
                                                                             const mipp::reg &rqjY,
                                                                             const mipp::reg &rqjZ,
                                                                                   mipp::reg &rajX,
                                                                                   mipp::reg &rajY,
                                                                                   mipp::reg &rajZ,
                                                                                   mipp::reg &rclosNeighj)
{
	mipp::reg rrijX = mipp::sub<float>(rqjX, rqiX); // 1 flop
	mipp::reg rrijY = mipp::sub<float>(rqjY, rqiY); // 1 flop
	mipp::reg rrijZ = mipp::sub<float>(rqjZ, rqiZ); // 1 flop

	mipp::reg rrijSquared = mipp::set1<float>(0.f);	
	rrijSquared = mipp::fmadd<float>(rrijX, rrijX, rrijSquared); // 2 flops
	rrijSquared = mipp::fmadd<float>(rrijY, rrijY, rrijSquared); // 2 flops
	rrijSquared = mipp::fmadd<float>(rrijZ, rrijZ, rrijSquared); // 2 flops

	mipp::reg rrijInv = mipp::rsqrt<float>(rrijSquared); // 1 flop

	mipp::reg raTmp = mipp::mul<float>(rG, mipp::mul<float>(mipp::mul<float>(rrijInv, rrijInv), rrijInv)); // 3 flops

	mipp::reg rai = mipp::mul<float>(raTmp, rmj); // 1 flop

	raiX = mipp::fmadd<float>(rai, rrijX, raiX); // 2 flops
	raiY = mipp::fmadd<float>(rai, rrijY, raiY); // 2 flops
	raiZ = mipp::fmadd<float>(rai, rrijZ, raiZ); // 2 flops

	mipp::reg raj = mipp::mul<float>(raTmp, rmi); // 1 flop

	rajX = mipp::fnmadd<float>(raj, rrijX, rajX); // 2 flops
	rajY = mipp::fnmadd<float>(raj, rrijY, rajY); // 2 flops
	rajZ = mipp::fnmadd<float>(raj, rrijZ, rajZ); // 2 flops
	
	rclosNeighi = mipp::max<float>(rrijInv, rclosNeighi);
	rclosNeighj = mipp::max<float>(rrijInv, rclosNeighj); // TODO: this second minimum is not threads safe !!!
}
