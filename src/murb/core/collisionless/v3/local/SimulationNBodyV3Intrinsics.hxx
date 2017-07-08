/*!
 * \file    SimulationNBodyV3Intrinsics.hxx
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

#include "SimulationNBodyV3Intrinsics.h"

template <typename T>
SimulationNBodyV3Intrinsics<T>::SimulationNBodyV3Intrinsics(const unsigned long nBodies, T softening)
	: SimulationNBodyV3<T>(nBodies, softening)
{
	this->init();
}

template <typename T>
SimulationNBodyV3Intrinsics<T>::SimulationNBodyV3Intrinsics(const std::string inputFileName, T softening)
	: SimulationNBodyV3<T>(inputFileName, softening)
{
	this->init();
}

template <typename T>
void SimulationNBodyV3Intrinsics<T>::init()
{
	this->flopsPerIte = 19.f * ((float)this->bodies->getN() -1.f) * (float)this->bodies->getN();
}

template <>
void SimulationNBodyV3Intrinsics<float>::init()
{
	this->flopsPerIte = 20.f * ((float)this->bodies->getN() -1.f) * (float)this->bodies->getN();
}

template <typename T>
SimulationNBodyV3Intrinsics<T>::~SimulationNBodyV3Intrinsics()
{
}

template <typename T>
void SimulationNBodyV3Intrinsics<T>::initIteration()
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
void SimulationNBodyV3Intrinsics<float>::initIteration()
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
void SimulationNBodyV3Intrinsics<T>::_computeLocalBodiesAcceleration()
{
	// TODO: be careful with the V1Intrinsics version: with fake bodies added at the end of the last vector, the
	//       dynamic time step is broken.
	//       It is necessary to launch the simulation with a number of bodies multiple of mipp::N<T>()!
	assert(this->dtConstant || (this->bodies->getN() % mipp::N<T>() == 0));

	const T *masses = this->getBodies()->getMasses();

	const T *positionsX = this->getBodies()->getPositionsX();
	const T *positionsY = this->getBodies()->getPositionsY();
	const T *positionsZ = this->getBodies()->getPositionsZ();

	const mipp::reg rG           = mipp::set1<T>(this->G);
	const mipp::reg rSoftSquared = mipp::set1<T>(this->softeningSquared);

#pragma omp parallel for schedule(runtime) firstprivate(rG)
	for(unsigned long iVec = 0; iVec < this->bodies->getNVecs(); iVec++)
	{
		// load vectors
		const mipp::reg rqiX = mipp::load<T>(positionsX + iVec * mipp::N<T>());
		const mipp::reg rqiY = mipp::load<T>(positionsY + iVec * mipp::N<T>());
		const mipp::reg rqiZ = mipp::load<T>(positionsZ + iVec * mipp::N<T>());

		mipp::reg raiX = mipp::load<T>(this->accelerations.x + iVec * mipp::N<T>());
		mipp::reg raiY = mipp::load<T>(this->accelerations.y + iVec * mipp::N<T>());
		mipp::reg raiZ = mipp::load<T>(this->accelerations.z + iVec * mipp::N<T>());

		mipp::reg rclosNeighi = mipp::set1<T>(0.0);
		if(!this->dtConstant)
			rclosNeighi = mipp::load<T>(this->closestNeighborDist + iVec * mipp::N<T>());

		for(unsigned long jVec = 0; jVec < this->bodies->getNVecs(); jVec++)
		{
			// load vectors
			mipp::reg rmj  = mipp::load<T>(masses     + jVec * mipp::N<T>());
			mipp::reg rqjX = mipp::load<T>(positionsX + jVec * mipp::N<T>());
			mipp::reg rqjY = mipp::load<T>(positionsY + jVec * mipp::N<T>());
			mipp::reg rqjZ = mipp::load<T>(positionsZ + jVec * mipp::N<T>());

			for(unsigned short iRot = 0; iRot < mipp::N<T>(); iRot++)
			{
				SimulationNBodyV3Intrinsics<T>::computeAccelerationBetweenTwoBodies(rG, rSoftSquared,
				                                                                    rqiX, rqiY, rqiZ,
				                                                                    raiX, raiY, raiZ,
				                                                                    rclosNeighi,
				                                                                    rmj,
				                                                                    rqjX, rqjY, rqjZ);

				// we make one useless rotate in the last iteration...
				rmj  = mipp::rrot<T>(rmj);
				rqjX = mipp::rrot<T>(rqjX); rqjY = mipp::rrot<T>(rqjY); rqjZ = mipp::rrot<T>(rqjZ);
			}
		}

		// store vectors
		mipp::store<T>(this->accelerations.x + iVec * mipp::N<T>(), raiX);
		mipp::store<T>(this->accelerations.y + iVec * mipp::N<T>(), raiY);
		mipp::store<T>(this->accelerations.z + iVec * mipp::N<T>(), raiZ);
		if(!this->dtConstant)
			mipp::store<T>(this->closestNeighborDist + iVec * mipp::N<T>(), rclosNeighi);
	}
}

template <typename T>
void SimulationNBodyV3Intrinsics<T>::computeLocalBodiesAcceleration()
{
	this->_computeLocalBodiesAcceleration();
}

template <>
void SimulationNBodyV3Intrinsics<float>::computeLocalBodiesAcceleration()
{
	this->_computeLocalBodiesAcceleration();

	for(unsigned long iBody = 0; iBody < this->bodies->getN(); iBody++)
		this->closestNeighborDist[iBody] = 1.0 / this->closestNeighborDist[iBody];
}

// 19 flops
template <typename T>
void SimulationNBodyV3Intrinsics<T>::computeAccelerationBetweenTwoBodies(const mipp::reg &rG,
                                                                         const mipp::reg &rSoftSquared,
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
                                                                         const mipp::reg &rqjZ)
{
	mipp::reg rrijX = mipp::sub<T>(rqjX, rqiX); // 1 flop
	mipp::reg rrijY = mipp::sub<T>(rqjY, rqiY); // 1 flop
	mipp::reg rrijZ = mipp::sub<T>(rqjZ, rqiZ); // 1 flop

	mipp::reg rrijSquared = rSoftSquared;
	rrijSquared = mipp::fmadd<T>(rrijX, rrijX, rrijSquared); // 2 flops
	rrijSquared = mipp::fmadd<T>(rrijY, rrijY, rrijSquared); // 2 flops
	rrijSquared = mipp::fmadd<T>(rrijZ, rrijZ, rrijSquared); // 2 flops

	mipp::reg rrij = mipp::sqrt<T>(rrijSquared); // 1 flop

	mipp::reg rai = mipp::div<T>(mipp::mul<T>(rG, rmj), mipp::mul<T>(rrij, rrijSquared)); // 3 flops

	raiX = mipp::fmadd<T>(rai, rrijX, raiX); // 2 flops
	raiY = mipp::fmadd<T>(rai, rrijY, raiY); // 2 flops
	raiZ = mipp::fmadd<T>(rai, rrijZ, raiZ); // 2 flops

	rclosNeighi = mipp::min<T>(rclosNeighi, rrij);
}

// 20 flops
template <>
void SimulationNBodyV3Intrinsics<float>::computeAccelerationBetweenTwoBodies(const mipp::reg &rG,
                                                                             const mipp::reg &rSoftSquared,
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
                                                                             const mipp::reg &rqjZ)
{
	mipp::reg rrijX = mipp::sub<float>(rqjX, rqiX); // 1 flop
	mipp::reg rrijY = mipp::sub<float>(rqjY, rqiY); // 1 flop
	mipp::reg rrijZ = mipp::sub<float>(rqjZ, rqiZ); // 1 flop

	mipp::reg rrijSquared = rSoftSquared;
	rrijSquared = mipp::fmadd<float>(rrijX, rrijX, rrijSquared); // 2 flops
	rrijSquared = mipp::fmadd<float>(rrijY, rrijY, rrijSquared); // 2 flops
	rrijSquared = mipp::fmadd<float>(rrijZ, rrijZ, rrijSquared); // 2 flops

	mipp::reg rrijInv = mipp::rsqrt<float>(rrijSquared); // 1 flop

	// || ai || = G * mj / (rij   * rij   * rij  ) <=>
	// || ai || = G * mj * (1/rij * 1/rij * 1/rij)
	mipp::reg rai = mipp::mul<float>(mipp::mul<float>(rG, rmj),
	                                 mipp::mul<float>(mipp::mul<float>(rrijInv, rrijInv), rrijInv)); // 4 flops

	raiX = mipp::fmadd<float>(rai, rrijX, raiX); // 2 flops
	raiY = mipp::fmadd<float>(rai, rrijY, raiY); // 2 flops
	raiZ = mipp::fmadd<float>(rai, rrijZ, raiZ); // 2 flops

	rclosNeighi = mipp::max<float>(rclosNeighi, rrijInv);
}
