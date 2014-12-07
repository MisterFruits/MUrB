/*!
 * \file    SimulationNBodyV3IntrinsicsBH.hxx
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

#include "SimulationNBodyV3IntrinsicsBH.h"

template <typename T>
SimulationNBodyV3IntrinsicsBH<T>::SimulationNBodyV3IntrinsicsBH(const unsigned long nBodies, T softening)
	: SimulationNBodyV3Intrinsics<T>(nBodies, softening)
{
	this->init();
}

template <typename T>
SimulationNBodyV3IntrinsicsBH<T>::SimulationNBodyV3IntrinsicsBH(const std::string inputFileName, T softening)
	: SimulationNBodyV3Intrinsics<T>(inputFileName, softening)
{
	this->init();
}

template <typename T>
void SimulationNBodyV3IntrinsicsBH<T>::init()
{
	this->flopsPerIte = 19 * (this->bodies->getN() -1) * this->bodies->getN();
}

template <>
void SimulationNBodyV3IntrinsicsBH<float>::init()
{
	this->flopsPerIte = 20 * (this->bodies->getN() -1) * this->bodies->getN();
}

template <typename T>
SimulationNBodyV3IntrinsicsBH<T>::~SimulationNBodyV3IntrinsicsBH()
{
}

template <typename T>
void SimulationNBodyV3IntrinsicsBH<T>::initIteration()
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
void SimulationNBodyV3IntrinsicsBH<float>::initIteration()
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
void SimulationNBodyV3IntrinsicsBH<T>::_computeLocalBodiesAcceleration()
{
	// TODO: be careful with the V1Intrinsics version: with fake bodies added at the end of the last vector, the
	//       dynamic time step is broken.
	//       It is necessary to launch the simulation with a number of bodies multiple of mipp::vectorSize<T>()!
	assert(this->dtConstant || (this->bodies->getN() % mipp::vectorSize<T>() == 0));

	const T *masses = this->getBodies()->getMasses();

	const T *positionsX = this->getBodies()->getPositionsX();
	const T *positionsY = this->getBodies()->getPositionsY();
	const T *positionsZ = this->getBodies()->getPositionsZ();

	const mipp::vec rG           = mipp::set1<T>(this->G);
	const mipp::vec rSoftSquared = mipp::set1<T>(this->softeningSquared);

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

			for(unsigned short iRot = 0; iRot < mipp::vectorSize<T>(); iRot++)
			{
				SimulationNBodyV3Intrinsics<T>::computeAccelerationBetweenTwoBodies(rG, rSoftSquared,
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

		// store vectors
		mipp::store<T>(this->accelerations.x + iVec * mipp::vectorSize<T>(), raiX);
		mipp::store<T>(this->accelerations.y + iVec * mipp::vectorSize<T>(), raiY);
		mipp::store<T>(this->accelerations.z + iVec * mipp::vectorSize<T>(), raiZ);
		if(!this->dtConstant)
			mipp::store<T>(this->closestNeighborDist + iVec * mipp::vectorSize<T>(), rclosNeighi);
	}
}

template <typename T>
void SimulationNBodyV3IntrinsicsBH<T>::computeLocalBodiesAccelerationWithBlackHole()
{
	const T *masses = this->getBodies()->getMasses();

	const T *positionsX = this->getBodies()->getPositionsX();
	const T *positionsY = this->getBodies()->getPositionsY();
	const T *positionsZ = this->getBodies()->getPositionsZ();

	T mbhGained = 0;
	for(unsigned long iBody; iBody < this->bodies->getN(); iBody++)
	{
		T dist = SimulationNBodyV3IntrinsicsBH<T>::computeAccelerationBetweenBodyAndBlackHole(
		                                                                              this->G, this->softeningSquared,
		                                                                              positionsX               [iBody],
		                                                                              positionsY               [iBody],
		                                                                              positionsZ               [iBody],
		                                                                              this->accelerations.x    [iBody],
		                                                                              this->accelerations.y    [iBody],
		                                                                              this->accelerations.z    [iBody],
		                                                                              this->mbh,
		                                                                              this->rbh,
		                                                                              this->qbhX,
		                                                                              this->qbhY,
		                                                                              this->qbhZ);

		if(dist <= 0)
			mbhGained += masses[iBody];
	}

	this->mbh += mbhGained;
}

template <typename T>
void SimulationNBodyV3IntrinsicsBH<T>::computeLocalBodiesAcceleration()
{
	this->_computeLocalBodiesAcceleration();
	this->computeLocalBodiesAccelerationWithBlackHole();
}

template <>
void SimulationNBodyV3IntrinsicsBH<float>::computeLocalBodiesAcceleration()
{
	this->_computeLocalBodiesAcceleration();
	this->computeLocalBodiesAccelerationWithBlackHole();

	for(unsigned long iBody = 0; iBody < this->bodies->getN(); iBody++)
		this->closestNeighborDist[iBody] = 1.0 / this->closestNeighborDist[iBody];
}


template <typename T>
T SimulationNBodyV3IntrinsicsBH<T>::computeAccelerationBetweenBodyAndBlackHole(const T &G,   const T &softSquared,
                                                                               const T &qiX, const T &qiY, const T &qiZ,
                                                                                     T &aiX,       T &aiY,       T &aiZ,
                                                                               const T &mbh,
                                                                               const T &rbh,
                                                                               const T &qbhX, const T &qbhY, const T &qbhZ)
{
	const T ribhX = qbhX - qiX; // 1 flop
	const T ribhY = qbhY - qiY; // 1 flop
	const T ribhZ = qbhZ - qiZ; // 1 flop

	// compute the distance (|| rij ||² + epsilon²) between body i and body j
	const T ribhSquared = (ribhX * ribhX) + (ribhY * ribhY) + (ribhZ * ribhZ) + softSquared; // 6 flops

	// compute the distance ~= || rij ||
	const T ribh = std::sqrt(ribhSquared); // 1 flop

	// compute the acceleration value between body i and body j: || ai || = G.mj / || rij ||^3
	const T ai = G * mbh / (ribhSquared * ribh); // 3 flops

	// add the acceleration value into the acceleration vector: ai += || ai ||.rij
	aiX += ai * ribhX; // 2 flops
	aiY += ai * ribhY; // 2 flops
	aiZ += ai * ribhZ; // 2 flops

	return ribh - rbh; // 1 flop
}
