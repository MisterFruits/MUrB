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

#include "../../../utils/myIntrinsicsPlusPlus.h"

#include "SimulationNBodyV1Vectors.h"

template <typename T>
SimulationNBodyV1Vectors<T>::SimulationNBodyV1Vectors(const unsigned long nBodies)
	: SimulationNBodyLocal<T>(nBodies)
{
	this->init();
}

template <typename T>
SimulationNBodyV1Vectors<T>::SimulationNBodyV1Vectors(const std::string inputFileName)
	: SimulationNBodyLocal<T>(inputFileName)
{
	this->init();
}

template <typename T>
void SimulationNBodyV1Vectors<T>::init()
{
	this->flopsPerIte = 18 * (this->bodies.getN() -1) * this->bodies.getN();
}

template <typename T>
SimulationNBodyV1Vectors<T>::~SimulationNBodyV1Vectors()
{
}

template <typename T>
void SimulationNBodyV1Vectors<T>::initIteration()
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
void SimulationNBodyV1Vectors<T>::computeLocalBodiesAcceleration()
{
	const T *masses = this->getBodies().getMasses();

	const T *positionsX = this->getBodies().getPositionsX();
	const T *positionsY = this->getBodies().getPositionsY();
	const T *positionsZ = this->getBodies().getPositionsZ();

#pragma omp parallel for schedule(runtime)
	for(unsigned long iVec = 0; iVec < this->bodies.getNVecs(); iVec++)
		for(unsigned long jVec = 0; jVec < this->bodies.getNVecs(); jVec++)
			if(iVec != jVec)
				for(unsigned short iVecPos = 0; iVecPos < mipp::vectorSize<T>(); iVecPos++)
				{
					const unsigned long iBody = iVecPos + iVec * mipp::vectorSize<T>();
					for(unsigned short jVecPos = 0; jVecPos < mipp::vectorSize<T>(); jVecPos++)
					{
						const unsigned long jBody = jVecPos + jVec * mipp::vectorSize<T>();
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
				}
			else
				for(unsigned short iVecPos = 0; iVecPos < mipp::vectorSize<T>(); iVecPos++)
				{
					const unsigned long iBody = iVecPos + iVec * mipp::vectorSize<T>();
					for(unsigned short jVecPos = 0; jVecPos < mipp::vectorSize<T>(); jVecPos++)
					{
						const unsigned long jBody = jVecPos + jVec * mipp::vectorSize<T>();
						if(iVecPos != jVecPos)
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
				}
}

// 18 flops
template <typename T>
void SimulationNBodyV1Vectors<T>::computeAccelerationBetweenTwoBodies(const T &iPosX, const T &iPosY, const T &iPosZ,
                                                                            T &iAccsX,      T &iAccsY,      T &iAccsZ,
                                                                            T &iClosNeiDist,
                                                                      const T &jMasses,
                                                                      const T &jPosX, const T &jPosY, const T &jPosZ)
{
	const T diffPosX = jPosX - iPosX; // 1 flop
	const T diffPosY = jPosY - iPosY; // 1 flop
	const T diffPosZ = jPosZ - iPosZ; // 1 flop
	const T squareDist = (diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ); // 5 flops
	const T dist = std::sqrt(squareDist); // 1 flop

	const T acc = this->G * jMasses / (squareDist * dist); // 3 flops
	iAccsX += acc * diffPosX; // 2 flop
	iAccsY += acc * diffPosY; // 2 flop
	iAccsZ += acc * diffPosZ; // 2 flop

	/*
	if(!this->dtConstant)
		// min
		if(dist < iClosNeiDist)
			iClosNeiDist = dist;
	*/
	iClosNeiDist = std::min(iClosNeiDist, dist);
}
