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

#include "../utils/myIntrinsicsPlusPlus.h"

#include "SimulationNBodyCollisionLocal.h"

template <typename T>
SimulationNBodyCollisionLocal<T>::SimulationNBodyCollisionLocal(const unsigned long nBodies)
	: SimulationNBody<T>(nBodies), collisions(this->bodies.getN())
{
}

template <typename T>
SimulationNBodyCollisionLocal<T>::SimulationNBodyCollisionLocal(const std::string inputFileName)
	: SimulationNBody<T>(inputFileName), collisions(this->bodies.getN())
{
}

template <typename T>
SimulationNBodyCollisionLocal<T>::~SimulationNBodyCollisionLocal()
{
}

template <typename T>
void SimulationNBodyCollisionLocal<T>::computeOneIteration()
{
	this->initIteration();
	this->computeLocalBodiesAcceleration();
	if(!this->dtConstant)
		this->findTimeStep();

	this->bodies.applyCollisions(this->collisions);

	this->bodies.updatePositionsAndVelocities(this->accelerations, this->dt);
}

template <typename T>
void SimulationNBodyCollisionLocal<T>::findTimeStep()
{
	// TODO: be careful with the V1Intrinsics version: with fake bodies added at the end of the last vector, the
	//       dynamic time step is broken.
	//       It is necessary to launch the simulation with a number of bodies multiple of mipp::vectorSize<T>()!
	if(!this->dtConstant)
	{
		this->dt = std::numeric_limits<T>::infinity();
		for(unsigned long iBody = 0; iBody < this->bodies.getN(); iBody++)
		{
			const T newDt = this->computeTimeStep(iBody);

			if(newDt < this->dt)
				this->dt = newDt;
		}
	}
}