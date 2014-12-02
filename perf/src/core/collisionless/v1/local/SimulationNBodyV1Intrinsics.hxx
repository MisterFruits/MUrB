/*!
 * \file    SimulationNBodyV1Intrinsics.hxx
 * \brief   Implementation of SimulationNBodyLocal with intrinsic function calls (n² computations).
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

#include "SimulationNBodyV1Intrinsics.h"

template <typename T>
SimulationNBodyV1Intrinsics<T>::SimulationNBodyV1Intrinsics(const unsigned long nBodies)
	: SimulationNBodyLocal<T>(nBodies)
{
	this->init();
}

template <typename T>
SimulationNBodyV1Intrinsics<T>::SimulationNBodyV1Intrinsics(const std::string inputFileName)
	: SimulationNBodyLocal<T>(inputFileName)
{
	this->init();
}

template <typename T>
void SimulationNBodyV1Intrinsics<T>::init()
{
	this->flopsPerIte = 19 * (this->bodies->getN() -1) * this->bodies->getN();
}

template <>
void SimulationNBodyV1Intrinsics<float>::init()
{
	this->flopsPerIte = 20 * (this->bodies->getN() -1) * this->bodies->getN();
}

template <typename T>
SimulationNBodyV1Intrinsics<T>::~SimulationNBodyV1Intrinsics()
{
}

template <typename T>
void SimulationNBodyV1Intrinsics<T>::initIteration()
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
void SimulationNBodyV1Intrinsics<float>::initIteration()
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
void SimulationNBodyV1Intrinsics<T>::_computeLocalBodiesAcceleration()
{
	// TODO: be careful with the V1Intrinsics version: with fake bodies added at the end of the last vector, the
	//       dynamic time step is broken.
	//       It is necessary to launch the simulation with a number of bodies multiple of mipp::vectorSize<T>()!
	assert(this->dtConstant || (this->bodies->getN() % mipp::vectorSize<T>() == 0));

	const T *masses = this->getBodies()->getMasses();

	const T *positionsX = this->getBodies()->getPositionsX();
	const T *positionsY = this->getBodies()->getPositionsY();
	const T *positionsZ = this->getBodies()->getPositionsZ();

	const mipp::vec rG = mipp::set1<T>(this->G);

#pragma omp parallel for schedule(runtime) firstprivate(rG)
	for(unsigned long iVec = 0; iVec < this->bodies->getNVecs(); iVec++)
	{
		// load vectors
		const mipp::vec rqiX = mipp::load<T>(positionsX + iVec * mipp::vectorSize<T>());
		const mipp::vec rqiY = mipp::load<T>(positionsY + iVec * mipp::vectorSize<T>());
		const mipp::vec rqiZ = mipp::load<T>(positionsZ + iVec * mipp::vectorSize<T>());

		mipp::vec raiX = mipp::load<T>(this->accelerations.x + iVec * mipp::vectorSize<T>());
		mipp::vec raiY = mipp::load<T>(this->accelerations.y + iVec * mipp::vectorSize<T>());
		mipp::vec raiZ = mipp::load<T>(this->accelerations.z + iVec * mipp::vectorSize<T>());

		mipp::vec rclosNeighi = mipp::set1<T>(0.0);
		if(!this->dtConstant)
			rclosNeighi = mipp::load<T>(this->closestNeighborDist + iVec * mipp::vectorSize<T>());

		for(unsigned long jVec = 0; jVec < this->bodies->getNVecs(); jVec++)
		{
			// load vectors
			mipp::vec rmj  = mipp::load<T>(masses     + jVec * mipp::vectorSize<T>());
			mipp::vec rqjX = mipp::load<T>(positionsX + jVec * mipp::vectorSize<T>());
			mipp::vec rqjY = mipp::load<T>(positionsY + jVec * mipp::vectorSize<T>());
			mipp::vec rqjZ = mipp::load<T>(positionsZ + jVec * mipp::vectorSize<T>());

			if(iVec != jVec)
			{
				for(unsigned short iRot = 0; iRot < mipp::vectorSize<T>(); iRot++)
				{
					SimulationNBodyV1Intrinsics<T>::computeAccelerationBetweenTwoBodies(rG,
					                                                                    rqiX, rqiY, rqiZ,
					                                                                    raiX, raiY, raiZ,
					                                                                    rclosNeighi,
					                                                                    rmj,
					                                                                    rqjX, rqjY, rqjZ);

					// we make one useless rotate in the last iteration...
					rmj  = mipp::rot<T>(rmj);
					rqjX = mipp::rot<T>(rqjX); rqjY = mipp::rot<T>(rqjY); rqjZ = mipp::rot<T>(rqjZ);
				}
			}
			else
			{
				for(unsigned short iRot = 1; iRot < mipp::vectorSize<T>(); iRot++)
				{
					rmj  = mipp::rot<T>(rmj);
					rqjX = mipp::rot<T>(rqjX); rqjY = mipp::rot<T>(rqjY); rqjZ = mipp::rot<T>(rqjZ);

					SimulationNBodyV1Intrinsics<T>::computeAccelerationBetweenTwoBodies(rG,
					                                                                    rqiX, rqiY, rqiZ,
					                                                                    raiX, raiY, raiZ,
					                                                                    rclosNeighi,
					                                                                    rmj,
					                                                                    rqjX, rqjY, rqjZ);
				}
			}
		}

		// store vectors
		mipp::store<T>(this->accelerations.x + iVec * mipp::vectorSize<T>(), raiX);
		mipp::store<T>(this->accelerations.y + iVec * mipp::vectorSize<T>(), raiY);
		mipp::store<T>(this->accelerations.z + iVec * mipp::vectorSize<T>(), raiZ);
		if(!this->dtConstant)
			mipp::store<T>(this->closestNeighborDist + iVec * mipp::vectorSize<T>(), rclosNeighi);
	}
}

template <typename T>
void SimulationNBodyV1Intrinsics<T>::computeLocalBodiesAcceleration()
{
	this->_computeLocalBodiesAcceleration();
}

template <>
void SimulationNBodyV1Intrinsics<float>::computeLocalBodiesAcceleration()
{
	this->_computeLocalBodiesAcceleration();

	for(unsigned long iBody = 0; iBody < this->bodies->getN(); iBody++)
		this->closestNeighborDist[iBody] = 1.0 / this->closestNeighborDist[iBody];
}

// 19 flops
template <typename T>
void SimulationNBodyV1Intrinsics<T>::computeAccelerationBetweenTwoBodies(const mipp::vec &rG,
                                                                         const mipp::vec &rqiX,
                                                                         const mipp::vec &rqiY,
                                                                         const mipp::vec &rqiZ,
                                                                               mipp::vec &raiX,
                                                                               mipp::vec &raiY,
                                                                               mipp::vec &raiZ,
                                                                               mipp::vec &rclosNeighi,
                                                                         const mipp::vec &rmj,
                                                                         const mipp::vec &rqjX,
                                                                         const mipp::vec &rqjY,
                                                                         const mipp::vec &rqjZ)
{
	mipp::vec rrijX = mipp::sub<T>(rqjX, rqiX); // 1 flop
	mipp::vec rrijY = mipp::sub<T>(rqjY, rqiY); // 1 flop
	mipp::vec rrijZ = mipp::sub<T>(rqjZ, rqiZ); // 1 flop

	mipp::vec rrijSquared = mipp::set1<T>(0);
	rrijSquared = mipp::fmadd<T>(rrijX, rrijX, rrijSquared); // 2 flops
	rrijSquared = mipp::fmadd<T>(rrijY, rrijY, rrijSquared); // 2 flops
	rrijSquared = mipp::fmadd<T>(rrijZ, rrijZ, rrijSquared); // 2 flops

	mipp::vec rrij = mipp::sqrt<T>(rrijSquared); // 1 flop

	mipp::vec rai = mipp::div<T>(mipp::mul<T>(rG, rmj), mipp::mul<T>(rrij, rrijSquared)); // 3 flops

	raiX = mipp::fmadd<T>(rai, rrijX, raiX); // 2 flops
	raiY = mipp::fmadd<T>(rai, rrijY, raiY); // 2 flops
	raiZ = mipp::fmadd<T>(rai, rrijZ, raiZ); // 2 flops

	rclosNeighi = mipp::min<T>(rclosNeighi, rrij);
}

// 20 flops
template <>
void SimulationNBodyV1Intrinsics<float>::computeAccelerationBetweenTwoBodies(const mipp::vec &rG,
                                                                             const mipp::vec &rqiX,
                                                                             const mipp::vec &rqiY,
                                                                             const mipp::vec &rqiZ,
                                                                                   mipp::vec &raiX,
                                                                                   mipp::vec &raiY,
                                                                                   mipp::vec &raiZ,
                                                                                   mipp::vec &rclosNeighi,
                                                                             const mipp::vec &rmj,
                                                                             const mipp::vec &rqjX,
                                                                             const mipp::vec &rqjY,
                                                                             const mipp::vec &rqjZ)
{
	mipp::vec rrijX = mipp::sub<float>(rqjX, rqiX); // 1 flop
	mipp::vec rrijY = mipp::sub<float>(rqjY, rqiY); // 1 flop
	mipp::vec rrijZ = mipp::sub<float>(rqjZ, rqiZ); // 1 flop

	mipp::vec rrijSquared = mipp::set1<float>(0);
	rrijSquared = mipp::fmadd<float>(rrijX, rrijX, rrijSquared); // 2 flops
	rrijSquared = mipp::fmadd<float>(rrijY, rrijY, rrijSquared); // 2 flops
	rrijSquared = mipp::fmadd<float>(rrijZ, rrijZ, rrijSquared); // 2 flops

	mipp::vec rrijInv = mipp::rsqrt<float>(rrijSquared); // 1 flop

	// || ai || = G * mj / (rij   * rij   * rij  ) <=>
	// || ai || = G * mj * (1/rij * 1/rij * 1/rij)
	mipp::vec rai = mipp::mul<float>(mipp::mul<float>(rG, rmj),
	                                 mipp::mul<float>(mipp::mul<float>(rrijInv, rrijInv), rrijInv)); // 4 flops

	raiX = mipp::fmadd<float>(rai, rrijX, raiX); // 2 flops
	raiY = mipp::fmadd<float>(rai, rrijY, raiY); // 2 flops
	raiZ = mipp::fmadd<float>(rai, rrijZ, raiZ); // 2 flops

	rclosNeighi = mipp::max<float>(rclosNeighi, rrijInv);
}
