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
#ifndef NO_OMP
#define NO_OMP
inline void omp_set_num_threads(int) {           }
inline int  omp_get_num_threads(   ) { return 1; }
inline int  omp_get_max_threads(   ) { return 1; }
inline int  omp_get_thread_num (   ) { return 0; }
#endif
#endif

#include "../../utils/myIntrinsicsPlusPlus.h"

#include "SimulationNBodyV2FineTuned.h"

template <typename T>
SimulationNBodyV2FineTuned<T>::SimulationNBodyV2FineTuned(const unsigned long nBodies)
	: SimulationNBodyLocal<T>(nBodies)
{
	this->reAllocateBuffers();
}

template <typename T>
SimulationNBodyV2FineTuned<T>::SimulationNBodyV2FineTuned(const std::string inputFileName)
	: SimulationNBodyLocal<T>(inputFileName)
{
	this->reAllocateBuffers();
}

template <typename T>
SimulationNBodyV2FineTuned<T>::~SimulationNBodyV2FineTuned()
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
void SimulationNBodyV2FineTuned<T>::_reAllocateBuffers()
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
}


template <typename T>
void SimulationNBodyV2FineTuned<T>::reAllocateBuffers()
{
	this->_reAllocateBuffers();
	this->flopsPerIte = 26 * (this->bodies.getN() * 0.5) * this->bodies.getN();
}

template <>
void SimulationNBodyV2FineTuned<float>::reAllocateBuffers()
{
	this->_reAllocateBuffers();
	this->flopsPerIte = 27 * (this->bodies.getN() * 0.5) * this->bodies.getN();
}

template <typename T>
void SimulationNBodyV2FineTuned<T>::_initIteration()
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
void SimulationNBodyV2FineTuned<T>::initIteration()
{
	this->_initIteration();

	for(unsigned long iBody = 0; iBody < this->bodies.getN(); iBody++)
		this->closestNeighborDist[iBody] = std::numeric_limits<T>::infinity();
}

template <>
void SimulationNBodyV2FineTuned<float>::initIteration()
{
	this->_initIteration();

	for(unsigned long iBody = 0; iBody < this->bodies.getN(); iBody++)
		this->closestNeighborDist[iBody] = 0;
}

template <typename T>
void SimulationNBodyV2FineTuned<T>::_computeLocalBodiesAcceleration()
{
	const T *masses = this->getBodies().getMasses();

	const T *positionsX = this->getBodies().getPositionsX();
	const T *positionsY = this->getBodies().getPositionsY();
	const T *positionsZ = this->getBodies().getPositionsZ();

  const mipp::vec rG = mipp::set1<T>(this->G);

#pragma omp parallel firstprivate(rG)
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
				this->computeAccelerationBetweenTwoBodiesSelf(positionsX               [iBody           ],
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
		
		// load vectors
		const mipp::vec rIMass = mipp::load<T>(masses     + iVecOff);
		const mipp::vec rIPosX = mipp::load<T>(positionsX + iVecOff);
		const mipp::vec rIPosY = mipp::load<T>(positionsY + iVecOff);
		const mipp::vec rIPosZ = mipp::load<T>(positionsZ + iVecOff);

		mipp::vec rIAccX = mipp::load<T>(this->accelerations.x + iVecOff + thStride);
		mipp::vec rIAccY = mipp::load<T>(this->accelerations.y + iVecOff + thStride);
		mipp::vec rIAccZ = mipp::load<T>(this->accelerations.z + iVecOff + thStride);
		mipp::vec rIClosNeiDist;
		if(!this->dtConstant)
			rIClosNeiDist = mipp::load<T>(this->closestNeighborDist + iVecOff);
		
		// computation of the vector number iVec with the following other vectors
		for(unsigned long jVec = iVec +1; jVec < this->bodies.getNVecs(); jVec++)
		{
			const unsigned long jVecOff = jVec * mipp::vectorSize<T>();
			mipp::vec rJMass = mipp::load<T>(masses     + jVecOff);
			mipp::vec rJPosX = mipp::load<T>(positionsX + jVecOff); 
			mipp::vec rJPosY = mipp::load<T>(positionsY + jVecOff);
			mipp::vec rJPosZ = mipp::load<T>(positionsZ + jVecOff);

			mipp::vec rJAccX = mipp::load<T>(this->accelerations.x + jVecOff + thStride);
			mipp::vec rJAccY = mipp::load<T>(this->accelerations.y + jVecOff + thStride);
			mipp::vec rJAccZ = mipp::load<T>(this->accelerations.z + jVecOff + thStride);
			mipp::vec rJClosNeiDist;
			if(!this->dtConstant)
				rJClosNeiDist = mipp::load<T>(this->closestNeighborDist + jVecOff);
	
			for(unsigned short iRot = 0; iRot < mipp::vectorSize<T>(); iRot++)
			{
				this->computeAccelerationBetweenTwoBodies(rG,
				                                          rIPosX, rIPosY, rIPosZ,
				                                          rIAccX, rIAccY, rIAccZ,
				                                          rIClosNeiDist,
				                                          rIMass,
				                                          rJPosX, rJPosY, rJPosZ,
				                                          rJAccX, rJAccY, rJAccZ,
				                                          rJClosNeiDist,
				                                          rJMass, this->dtConstant);

				rJMass = mipp::rot<T>(rJMass);
				rJPosX = mipp::rot<T>(rJPosX); rJPosY = mipp::rot<T>(rJPosY); rJPosZ = mipp::rot<T>(rJPosZ);
				rJAccX = mipp::rot<T>(rJAccX); rJAccY = mipp::rot<T>(rJAccY); rJAccZ = mipp::rot<T>(rJAccZ);
				rJClosNeiDist = mipp::rot<T>(rJClosNeiDist);
			}
			
			mipp::store<T>(this->accelerations.x + jVecOff + thStride, rJAccX);
			mipp::store<T>(this->accelerations.y + jVecOff + thStride, rJAccY);
			mipp::store<T>(this->accelerations.z + jVecOff + thStride, rJAccZ);
		
			if(!this->dtConstant)
			{
#pragma omp critical
				mipp::store<T>(this->closestNeighborDist + jVecOff, rJClosNeiDist);
			}
		}

		mipp::store<T>(this->accelerations.x + iVecOff + thStride, rIAccX);
		mipp::store<T>(this->accelerations.y + iVecOff + thStride, rIAccY);
		mipp::store<T>(this->accelerations.z + iVecOff + thStride, rIAccZ);

		if(!this->dtConstant)
			mipp::store<T>(this->closestNeighborDist + iVecOff, rIClosNeiDist);
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


template <typename T>
void SimulationNBodyV2FineTuned<T>::computeLocalBodiesAcceleration()
{
	this->_computeLocalBodiesAcceleration();
}

template <>
void SimulationNBodyV2FineTuned<float>::computeLocalBodiesAcceleration()
{
	this->_computeLocalBodiesAcceleration();

	for(unsigned long iBody = 0; iBody < this->bodies.getN(); iBody++)
		this->closestNeighborDist[iBody] = 1.0 / this->closestNeighborDist[iBody];
}

// 25 flops
template <typename T>
void SimulationNBodyV2FineTuned<T>::computeAccelerationBetweenTwoBodiesSelf(const T &iPosX, const T &iPosY, const T &iPosZ,
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

// 26 flops
template <typename T>
void SimulationNBodyV2FineTuned<T>::computeAccelerationBetweenTwoBodies(const mipp::vec &rG, 
                                                                         const mipp::vec &rIPosX,
                                                                         const mipp::vec &rIPosY,
                                                                         const mipp::vec &rIPosZ,
                                                                               mipp::vec &rIAccX,
                                                                               mipp::vec &rIAccY,
                                                                               mipp::vec &rIAccZ,
                                                                               mipp::vec &rIClosNeiDist,
                                                                         const mipp::vec &rIMass,
                                                                         const mipp::vec &rJPosX,
                                                                         const mipp::vec &rJPosY,
                                                                         const mipp::vec &rJPosZ,
                                                                               mipp::vec &rJAccX,
                                                                               mipp::vec &rJAccY,
                                                                               mipp::vec &rJAccZ,
                                                                               mipp::vec &rJClosNeiDist,
                                                                         const mipp::vec &rJMass,
                                                                         const bool dtConstant)
{
	//const T diffPosX = jPosX - iPosX; // 1 flop
	mipp::vec rDiffPosX = mipp::sub<T>(rJPosX, rIPosX);
	//const T diffPosY = jPosY - iPosY; // 1 flop
	mipp::vec rDiffPosY = mipp::sub<T>(rJPosY, rIPosY);
	//const T diffPosZ = jPosZ - iPosZ; // 1 flop
	mipp::vec rDiffPosZ = mipp::sub<T>(rJPosZ, rIPosZ);

	//const T squareDist += (diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ); // 6 flops
	mipp::vec rSquareDist = mipp::set1<T>(0);	
	rSquareDist = mipp::fmadd<T>(rDiffPosX, rDiffPosX, rSquareDist); // 2 flops
	rSquareDist = mipp::fmadd<T>(rDiffPosY, rDiffPosY, rSquareDist); // 2 flops
	rSquareDist = mipp::fmadd<T>(rDiffPosZ, rDiffPosZ, rSquareDist); // 2 flops

	//const T dist = std::sqrt(squareDist); // 1 flop
	mipp::vec rDist = mipp::sqrt<T>(rSquareDist);

	//const T force = this->G / (squareDist * dist); // 2 flops
	mipp::vec rForce = mipp::div<T>(rG, mipp::mul<T>(rDist, rSquareDist));

	//T acc = force * jMasses; // 1 flop
	mipp::vec rAcc = mipp::mul<T>(rForce, rJMass);

	//iAccsX += acc * diffPosX; // 2 flops
	rIAccX = mipp::fmadd<T>(rAcc, rDiffPosX, rIAccX);
	//iAccsY += acc * diffPosY; // 2 flops
	rIAccY = mipp::fmadd<T>(rAcc, rDiffPosY, rIAccY);
	//iAccsZ += acc * diffPosZ; // 2 flops
	rIAccZ = mipp::fmadd<T>(rAcc, rDiffPosZ, rIAccZ);

	//acc = force * iMasses; // 1 flop
	rAcc = mipp::mul<T>(rForce, rIMass);

	//jAccsX -= acc * diffPosX; // 2 flop
	rJAccX = mipp::fnmadd<T>(rAcc, rDiffPosX, rJAccX);
	//jAccsY -= acc * diffPosY; // 2 flop
	rJAccY = mipp::fnmadd<T>(rAcc, rDiffPosY, rJAccY);
	//jAccsZ -= acc * diffPosZ; // 2 flop
	rJAccZ = mipp::fnmadd<T>(rAcc, rDiffPosZ, rJAccZ);

	//  min(iClosNeiDist, dist);
	if(!dtConstant)
		rIClosNeiDist = mipp::min<T>(rDist, rIClosNeiDist);
}

// 27 flops
template <>
void SimulationNBodyV2FineTuned<float>::computeAccelerationBetweenTwoBodies(const mipp::vec &rG, 
                                                                             const mipp::vec &rIPosX,
                                                                             const mipp::vec &rIPosY,
                                                                             const mipp::vec &rIPosZ,
                                                                                   mipp::vec &rIAccX,
                                                                                   mipp::vec &rIAccY,
                                                                                   mipp::vec &rIAccZ,
                                                                                   mipp::vec &rIClosNeiDist,
                                                                             const mipp::vec &rIMass,
                                                                             const mipp::vec &rJPosX,
                                                                             const mipp::vec &rJPosY,
                                                                             const mipp::vec &rJPosZ,
                                                                                   mipp::vec &rJAccX,
                                                                                   mipp::vec &rJAccY,
                                                                                   mipp::vec &rJAccZ,
                                                                                   mipp::vec &rJClosNeiDist,
                                                                             const mipp::vec &rJMass,
                                                                             const bool dtConstant)
{
	//const T diffPosX = jPosX - iPosX; // 1 flop
	mipp::vec rDiffPosX = mipp::sub<float>(rJPosX, rIPosX);
	//const T diffPosY = jPosY - iPosY; // 1 flop
	mipp::vec rDiffPosY = mipp::sub<float>(rJPosY, rIPosY);
	//const T diffPosZ = jPosZ - iPosZ; // 1 flop
	mipp::vec rDiffPosZ = mipp::sub<float>(rJPosZ, rIPosZ);

	//const T squareDist += (diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ); // 6 flops
	mipp::vec rSquareDist = mipp::set1<float>(0);	
	rSquareDist = mipp::fmadd<float>(rDiffPosX, rDiffPosX, rSquareDist); 
	rSquareDist = mipp::fmadd<float>(rDiffPosY, rDiffPosY, rSquareDist);
	rSquareDist = mipp::fmadd<float>(rDiffPosZ, rDiffPosZ, rSquareDist);

	//const T invDist =  1.0 / std::sqrt(squareDist); // 1 flop
	mipp::vec rInvDist = mipp::rsqrt<float>(rSquareDist);

	//const T force = this->G / (squareDist * dist); // 3 flops
	mipp::vec rForce = mipp::mul<float>(rG, mipp::mul<float>(mipp::mul<float>(rInvDist, rInvDist), rInvDist));

	//T acc = force * jMasses; // 1 flop
	mipp::vec rAcc = mipp::mul<float>(rForce, rJMass);

	//iAccsX += acc * diffPosX; // 2 flops
	rIAccX = mipp::fmadd<float>(rAcc, rDiffPosX, rIAccX);
	//iAccsY += acc * diffPosY; // 2 flops
	rIAccY = mipp::fmadd<float>(rAcc, rDiffPosY, rIAccY);
	//iAccsZ += acc * diffPosZ; // 2 flops
	rIAccZ = mipp::fmadd<float>(rAcc, rDiffPosZ, rIAccZ);

	//acc = force * iMasses; // 1 flop
	rAcc = mipp::mul<float>(rForce, rIMass);

	//jAccsX -= acc * diffPosX; // 2 flop
	rJAccX = mipp::fnmadd<float>(rAcc, rDiffPosX, rJAccX);
	//jAccsY -= acc * diffPosY; // 2 flop
	rJAccY = mipp::fnmadd<float>(rAcc, rDiffPosY, rJAccY);
	//jAccsZ -= acc * diffPosZ; // 2 flop
	rJAccZ = mipp::fnmadd<float>(rAcc, rDiffPosZ, rJAccZ);
	
	//max(iClosNeiDist, invDist);
	if(!dtConstant)
		rIClosNeiDist = mipp::max<float>(rInvDist, rIClosNeiDist);
}
