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

#include "SimulationNBodyMPIV1Intrinsics.h"

template <typename T>
SimulationNBodyMPIV1Intrinsics<T>::SimulationNBodyMPIV1Intrinsics(const unsigned long nBodies)
	: SimulationNBodyMPI<T>(nBodies)
{
	this->init();
}

template <typename T>
SimulationNBodyMPIV1Intrinsics<T>::SimulationNBodyMPIV1Intrinsics(const std::string inputFileName)
	: SimulationNBodyMPI<T>(inputFileName)
{
	this->init();
}

template <typename T>
void SimulationNBodyMPIV1Intrinsics<T>::init()
{
	this->flopsPerIte = 19 * ((this->bodies->getN() * this->MPISize) -1) * (this->bodies->getN() * this->MPISize);
}

template <>
void SimulationNBodyMPIV1Intrinsics<float>::init()
{
	this->flopsPerIte = 20 * ((this->bodies->getN() * this->MPISize) -1) * (this->bodies->getN() * this->MPISize);
}

template <typename T>
SimulationNBodyMPIV1Intrinsics<T>::~SimulationNBodyMPIV1Intrinsics()
{
}

template <typename T>
void SimulationNBodyMPIV1Intrinsics<T>::initIteration()
{
	for(unsigned long iBody = 0; iBody < this->bodies->getN(); iBody++)
	{
		this->accelerations.x[iBody] = 0.0;
		this->accelerations.y[iBody] = 0.0;
		this->accelerations.z[iBody] = 0.0;

		this->closestNeighborDist[iBody] = std::numeric_limits<T>::infinity();
	}
}

template <>
void SimulationNBodyMPIV1Intrinsics<float>::initIteration()
{
	for(unsigned long iBody = 0; iBody < this->bodies->getN(); iBody++)
	{
		this->accelerations.x[iBody] = 0.0;
		this->accelerations.y[iBody] = 0.0;
		this->accelerations.z[iBody] = 0.0;

		this->closestNeighborDist[iBody] = 0.0;
	}
}

template <typename T>
void SimulationNBodyMPIV1Intrinsics<T>::_computeLocalBodiesAcceleration()
{
	// TODO: be careful with the V1Intrinsics version: with fake bodies added at the end of the last vector, the
	//       dynamic time step is broken.
	//       It is necessary to launch the simulation with a number of bodies multiple of mipp::vectorSize<T>()!
	assert(this->dtConstant || (this->bodies->getN() % mipp::vectorSize<T>() == 0));

	const T *masses = this->bodies->getMasses();

	const T *positionsX = this->bodies->getPositionsX();
	const T *positionsY = this->bodies->getPositionsY();
	const T *positionsZ = this->bodies->getPositionsZ();

	const mipp::vec rG = mipp::set1<T>(this->G);

#pragma omp parallel for schedule(runtime) firstprivate(rG)
	for(unsigned long iVec = 0; iVec < this->bodies->getNVecs(); iVec++)
	{
		// load vectors
		const mipp::vec rIPosX = mipp::load<T>(positionsX + iVec * mipp::vectorSize<T>());
		const mipp::vec rIPosY = mipp::load<T>(positionsY + iVec * mipp::vectorSize<T>());
		const mipp::vec rIPosZ = mipp::load<T>(positionsZ + iVec * mipp::vectorSize<T>());

		mipp::vec rIAccX = mipp::load<T>(this->accelerations.x + iVec * mipp::vectorSize<T>());
		mipp::vec rIAccY = mipp::load<T>(this->accelerations.y + iVec * mipp::vectorSize<T>());
		mipp::vec rIAccZ = mipp::load<T>(this->accelerations.z + iVec * mipp::vectorSize<T>());

		mipp::vec rIClosNeiDist = mipp::set1<T>(0.0);
		if(!this->dtConstant)
			rIClosNeiDist = mipp::load<T>(this->closestNeighborDist + iVec * mipp::vectorSize<T>());

		for(unsigned long jVec = 0; jVec < this->bodies->getNVecs(); jVec++)
		{
			// load vectors
			mipp::vec rJMass = mipp::load<T>(masses     + jVec * mipp::vectorSize<T>());
			mipp::vec rJPosX = mipp::load<T>(positionsX + jVec * mipp::vectorSize<T>());
			mipp::vec rJPosY = mipp::load<T>(positionsY + jVec * mipp::vectorSize<T>());
			mipp::vec rJPosZ = mipp::load<T>(positionsZ + jVec * mipp::vectorSize<T>());

			if(iVec != jVec)
			{
				for(unsigned short iRot = 0; iRot < mipp::vectorSize<T>(); iRot++)
				{
					this->computeAccelerationBetweenTwoBodies(rG,
					                                          rIPosX, rIPosY, rIPosZ,
					                                          rIAccX, rIAccY, rIAccZ,
					                                          rIClosNeiDist,
					                                          rJMass,
					                                          rJPosX, rJPosY, rJPosZ);

					// we make one useless rotate in the last iteration...
					rJMass = mipp::rot<T>(rJMass);
					rJPosX = mipp::rot<T>(rJPosX); rJPosY = mipp::rot<T>(rJPosY); rJPosZ = mipp::rot<T>(rJPosZ);
				}
			}
			else
			{
				for(unsigned short iRot = 1; iRot < mipp::vectorSize<T>(); iRot++)
				{
					rJMass = mipp::rot<T>(rJMass);
					rJPosX = mipp::rot<T>(rJPosX); rJPosY = mipp::rot<T>(rJPosY); rJPosZ = mipp::rot<T>(rJPosZ);

					this->computeAccelerationBetweenTwoBodies(rG,
					                                          rIPosX, rIPosY, rIPosZ,
					                                          rIAccX, rIAccY, rIAccZ,
					                                          rIClosNeiDist,
					                                          rJMass,
					                                          rJPosX, rJPosY, rJPosZ);
				}
			}
		}

		// store vectors
		mipp::store<T>(this->accelerations.x + iVec * mipp::vectorSize<T>(), rIAccX);
		mipp::store<T>(this->accelerations.y + iVec * mipp::vectorSize<T>(), rIAccY);
		mipp::store<T>(this->accelerations.z + iVec * mipp::vectorSize<T>(), rIAccZ);
		if(!this->dtConstant)
			mipp::store<T>(this->closestNeighborDist + iVec * mipp::vectorSize<T>(), rIClosNeiDist);
	}
}

template <typename T>
void SimulationNBodyMPIV1Intrinsics<T>::computeLocalBodiesAcceleration()
{
	this->_computeLocalBodiesAcceleration();
}

template <>
void SimulationNBodyMPIV1Intrinsics<float>::computeLocalBodiesAcceleration()
{
	this->_computeLocalBodiesAcceleration();

	for(unsigned long iBody = 0; iBody < this->bodies->getN(); iBody++)
		this->closestNeighborDist[iBody] = 1.0 / this->closestNeighborDist[iBody];
}

template <typename T>
void SimulationNBodyMPIV1Intrinsics<T>::_computeNeighborBodiesAcceleration()
{
	// TODO: be careful with the V1Intrinsics version: with fake bodies added at the end of the last vector, the
	//       dynamic time step is broken.
	//       It is necessary to launch the simulation with a number of bodies multiple of mipp::vectorSize<T>()!
	assert(this->dtConstant || (this->bodies->getN() % mipp::vectorSize<T>() == 0));

	const T *positionsX = this->bodies->getPositionsX();
	const T *positionsY = this->bodies->getPositionsY();
	const T *positionsZ = this->bodies->getPositionsZ();

	const T *neighMasses     = this->neighborBodies->getMasses();
	const T *neighPositionsX = this->neighborBodies->getPositionsX();
	const T *neighPositionsY = this->neighborBodies->getPositionsY();
	const T *neighPositionsZ = this->neighborBodies->getPositionsZ();

	const mipp::vec rG = mipp::set1<T>(this->G);

#pragma omp parallel for schedule(runtime) firstprivate(rG)
	for(unsigned long iVec = 0; iVec < this->bodies->getNVecs(); iVec++)
	{
		// load vectors
		const mipp::vec rIPosX = mipp::load<T>(positionsX + iVec * mipp::vectorSize<T>());
		const mipp::vec rIPosY = mipp::load<T>(positionsY + iVec * mipp::vectorSize<T>());
		const mipp::vec rIPosZ = mipp::load<T>(positionsZ + iVec * mipp::vectorSize<T>());

		mipp::vec rIAccX = mipp::load<T>(this->accelerations.x + iVec * mipp::vectorSize<T>());
		mipp::vec rIAccY = mipp::load<T>(this->accelerations.y + iVec * mipp::vectorSize<T>());
		mipp::vec rIAccZ = mipp::load<T>(this->accelerations.z + iVec * mipp::vectorSize<T>());

		mipp::vec rIClosNeiDist = mipp::set1<T>(0.0);
		if(!this->dtConstant)
			rIClosNeiDist = mipp::load<T>(this->closestNeighborDist + iVec * mipp::vectorSize<T>());

		for(unsigned long jVec = 0; jVec < this->bodies->getNVecs(); jVec++)
		{
			// load vectors
			mipp::vec rJMass = mipp::load<T>(neighMasses     + jVec * mipp::vectorSize<T>());
			mipp::vec rJPosX = mipp::load<T>(neighPositionsX + jVec * mipp::vectorSize<T>());
			mipp::vec rJPosY = mipp::load<T>(neighPositionsY + jVec * mipp::vectorSize<T>());
			mipp::vec rJPosZ = mipp::load<T>(neighPositionsZ + jVec * mipp::vectorSize<T>());

			for(unsigned short iRot = 0; iRot < mipp::vectorSize<T>(); iRot++)
			{
				this->computeAccelerationBetweenTwoBodies(rG,
				                                          rIPosX, rIPosY, rIPosZ,
				                                          rIAccX, rIAccY, rIAccZ,
				                                          rIClosNeiDist,
				                                          rJMass,
				                                          rJPosX, rJPosY, rJPosZ);

				// we make one useless rotate in the last iteration...
				rJMass = mipp::rot<T>(rJMass);
				rJPosX = mipp::rot<T>(rJPosX); rJPosY = mipp::rot<T>(rJPosY); rJPosZ = mipp::rot<T>(rJPosZ);
			}
		}

		// store vectors
		mipp::store<T>(this->accelerations.x + iVec * mipp::vectorSize<T>(), rIAccX);
		mipp::store<T>(this->accelerations.y + iVec * mipp::vectorSize<T>(), rIAccY);
		mipp::store<T>(this->accelerations.z + iVec * mipp::vectorSize<T>(), rIAccZ);
		if(!this->dtConstant)
			mipp::store<T>(this->closestNeighborDist + iVec * mipp::vectorSize<T>(), rIClosNeiDist);
	}
}

template <typename T>
void SimulationNBodyMPIV1Intrinsics<T>::computeNeighborBodiesAcceleration()
{
	this->_computeNeighborBodiesAcceleration();
}

template <>
void SimulationNBodyMPIV1Intrinsics<float>::computeNeighborBodiesAcceleration()
{
	this->_computeNeighborBodiesAcceleration();

	for(unsigned long iBody = 0; iBody < this->bodies->getN(); iBody++)
		this->closestNeighborDist[iBody] = 1.0 / this->closestNeighborDist[iBody];
}

// 19 flops
template <typename T>
void SimulationNBodyMPIV1Intrinsics<T>::computeAccelerationBetweenTwoBodies(const mipp::vec &rG,
                                                                            const mipp::vec &rIPosX,
                                                                            const mipp::vec &rIPosY,
                                                                            const mipp::vec &rIPosZ,
                                                                                  mipp::vec &rIAccX,
                                                                                  mipp::vec &rIAccY,
                                                                                  mipp::vec &rIAccZ,
                                                                                  mipp::vec &rIClosNeiDist,
                                                                            const mipp::vec &rJMass,
                                                                            const mipp::vec &rJPosX,
                                                                            const mipp::vec &rJPosY,
                                                                            const mipp::vec &rJPosZ)
{
	//const T diffPosX = jPosX - iPosX;
	mipp::vec rDiffPosX = mipp::sub<T>(rJPosX, rIPosX); // 1 flop
	//const T diffPosY = jPosY - iPosY;
	mipp::vec rDiffPosY = mipp::sub<T>(rJPosY, rIPosY); // 1 flop
	//const T diffPosZ = jPosZ - iPosZ;
	mipp::vec rDiffPosZ = mipp::sub<T>(rJPosZ, rIPosZ); // 1 flop

	//const T squareDist = (diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ);
	mipp::vec rSquareDist = mipp::set1<T>(0);
	rSquareDist = mipp::fmadd<T>(rDiffPosX, rDiffPosX, rSquareDist); // 2 flops
	rSquareDist = mipp::fmadd<T>(rDiffPosY, rDiffPosY, rSquareDist); // 2 flops
	rSquareDist = mipp::fmadd<T>(rDiffPosZ, rDiffPosZ, rSquareDist); // 2 flops

	//const T dist = std::sqrt(squareDist);
	mipp::vec rDist = mipp::sqrt<T>(rSquareDist); // 1 flop

	//const T acc = G * jMasses / (squareDist * dist);
	mipp::vec rAcc = mipp::div<T>(mipp::mul<T>(rG, rJMass), mipp::mul<T>(rDist, rSquareDist)); // 3 flops

	//iAccsX += acc * diffPosX;
	rIAccX = mipp::fmadd<T>(rAcc, rDiffPosX, rIAccX); // 2 flops
	//iAccsY += acc * diffPosY;
	rIAccY = mipp::fmadd<T>(rAcc, rDiffPosY, rIAccY); // 2 flops
	//iAccsZ += acc * diffPosZ;
	rIAccZ = mipp::fmadd<T>(rAcc, rDiffPosZ, rIAccZ); // 2 flops

	//if(!this->dtConstant)
	//	min(iClosNeiDist, dist);
	rIClosNeiDist = mipp::min<T>(rDist, rIClosNeiDist);
}

// 20 flops
template <>
void SimulationNBodyMPIV1Intrinsics<float>::computeAccelerationBetweenTwoBodies(const mipp::vec &rG,
                                                                                const mipp::vec &rIPosX,
                                                                                const mipp::vec &rIPosY,
                                                                                const mipp::vec &rIPosZ,
                                                                                      mipp::vec &rIAccX,
                                                                                      mipp::vec &rIAccY,
                                                                                      mipp::vec &rIAccZ,
                                                                                      mipp::vec &rIClosNeiDist,
                                                                                const mipp::vec &rJMass,
                                                                                const mipp::vec &rJPosX,
                                                                                const mipp::vec &rJPosY,
                                                                                const mipp::vec &rJPosZ)
{
	//const T diffPosX = jPosX - iPosX;
	mipp::vec rDiffPosX = mipp::sub<float>(rJPosX, rIPosX); // 1 flop
	//const T diffPosY = jPosY - iPosY;
	mipp::vec rDiffPosY = mipp::sub<float>(rJPosY, rIPosY); // 1 flop
	//const T diffPosZ = jPosZ - iPosZ;
	mipp::vec rDiffPosZ = mipp::sub<float>(rJPosZ, rIPosZ); // 1 flop

	//const T squareDist = (diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ);
	mipp::vec rSquareDist = mipp::set1<float>(0);
	rSquareDist = mipp::fmadd<float>(rDiffPosX, rDiffPosX, rSquareDist); // 2 flops
	rSquareDist = mipp::fmadd<float>(rDiffPosY, rDiffPosY, rSquareDist); // 2 flops
	rSquareDist = mipp::fmadd<float>(rDiffPosZ, rDiffPosZ, rSquareDist); // 2 flops

	//const T invDist = 1.0 / std::sqrt(squareDist);
	mipp::vec rInvDist = mipp::rsqrt<float>(rSquareDist); // 1 flop

	// const T acc = G * jMasses / (dist * dist * dist) <=>
	// const T acc = G * jMasses * (invDist * invDist * invDist)
	mipp::vec rAcc = mipp::mul<float>(mipp::mul<float>(rG, rJMass),
	                                  mipp::mul<float>(mipp::mul<float>(rInvDist, rInvDist), rInvDist)); // 4 flops

	//iAccsX += acc * diffPosX
	rIAccX = mipp::fmadd<float>(rAcc, rDiffPosX, rIAccX); // 2 flops
	//iAccsY += acc * diffPosY
	rIAccY = mipp::fmadd<float>(rAcc, rDiffPosY, rIAccY); // 2 flops
	//iAccsZ += acc * diffPosZ
	rIAccZ = mipp::fmadd<float>(rAcc, rDiffPosZ, rIAccZ); // 2 flops

	//if(!this->dtConstant)
	//	min(iClosNeiDist, dist) <=> max(iClosNeiDist, rInvDist)
	rIClosNeiDist = mipp::max<float>(rInvDist, rIClosNeiDist);
}
