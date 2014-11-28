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

#include "SimulationNBody.h"

template <typename T>
SimulationNBody<T>::SimulationNBody(Bodies<T> *bodies)
	: bodies             (bodies),
	  closestNeighborDist(NULL),
	  dtConstant         (false),
	  dt                 (std::numeric_limits<T>::infinity()),
	  minDt              (0),
	  flopsPerIte        (0),
	  allocatedBytes     (bodies->getAllocatedBytes()),
	  nMaxThreads        (omp_get_max_threads())
{
	assert(bodies != nullptr);

	this->allocateBuffers();
}

template <typename T>
SimulationNBody<T>::~SimulationNBody()
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

	if(this->closestNeighborDist != nullptr) {
		delete[] this->closestNeighborDist;
		this->closestNeighborDist = nullptr;
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

	if(this->closestNeighborDist != nullptr) {
		_mm_free(this->closestNeighborDist);
		this->closestNeighborDist = nullptr;
	}
#endif
}

template <typename T>
void SimulationNBody<T>::allocateBuffers()
{
#ifdef __ARM_NEON__
	this->accelerations.x = new T[this->bodies->getN() + this->bodies->getPadding()];
	this->accelerations.y = new T[this->bodies->getN() + this->bodies->getPadding()];
	this->accelerations.z = new T[this->bodies->getN() + this->bodies->getPadding()];

	this->closestNeighborDist = new T[this->bodies->getN() + this->bodies->getPadding()];
#else
	this->accelerations.x = (T*)_mm_malloc((this->bodies->getN() + this->bodies->getPadding()) * sizeof(T),
	                                       mipp::RequiredAlignement);
	this->accelerations.y = (T*)_mm_malloc((this->bodies->getN() + this->bodies->getPadding()) * sizeof(T),
	                                       mipp::RequiredAlignement);
	this->accelerations.z = (T*)_mm_malloc((this->bodies->getN() + this->bodies->getPadding()) * sizeof(T),
	                                       mipp::RequiredAlignement);

	this->closestNeighborDist = (T*)_mm_malloc((this->bodies->getN() + this->bodies->getPadding()) * sizeof(T),
	                                           mipp::RequiredAlignement);
#endif

	this->allocatedBytes += (this->bodies->getN() + this->bodies->getPadding()) * sizeof(T) * 4;
}

template <typename T>
const Bodies<T>* SimulationNBody<T>::getBodies() const
{
	return const_cast<const Bodies<T>*>(this->bodies);
}

template <typename T>
void SimulationNBody<T>::setDtConstant(T dtVal)
{
	this->dtConstant = true;
	this->dt = dtVal;
}

template <typename T>
void SimulationNBody<T>::setDtVariable(T minDt)
{
	this->dtConstant = false;
	this->dt = std::numeric_limits<T>::infinity();
	this->minDt = minDt;
}

template <typename T>
const T& SimulationNBody<T>::getDt() const
{
	return this->dt;
}

template <typename T>
const float& SimulationNBody<T>::getFlopsPerIte() const
{
	return this->flopsPerIte;
}

template <typename T>
const float& SimulationNBody<T>::getAllocatedBytes() const
{
	return this->allocatedBytes;
}

template <typename T>
T SimulationNBody<T>::computeTimeStep(const unsigned long iBody)
{
	const T *velocitiesX = this->bodies->getVelocitiesX();
	const T *velocitiesY = this->bodies->getVelocitiesY();
	const T *velocitiesZ = this->bodies->getVelocitiesZ();

	// || velocity[iBody] ||
	const T vi = std::sqrt((velocitiesX[iBody] * velocitiesX[iBody]) +
	                       (velocitiesY[iBody] * velocitiesY[iBody]) +
	                       (velocitiesZ[iBody] * velocitiesZ[iBody]));

	// || acceleration[iBody] ||
	const T ai = std::sqrt((this->accelerations.x[iBody] * this->accelerations.x[iBody]) +
	                       (this->accelerations.y[iBody] * this->accelerations.y[iBody]) +
	                       (this->accelerations.z[iBody] * this->accelerations.z[iBody]));

	/*
	 * compute dt
	 * solve:  (a/2)*dt^2 + v*dt + (-0.1)*ClosestNeighborDist = 0
	 * <=>     dt = [ (-v) +/-  sqrt( v^2 - 4 * (a/2) * (-0.1)*ClosestNeighborDist ) ] / [ 2 (a/2) ]
	 *
	 * dt should be positive (+/- becomes + because result of sqrt is positive)
	 * <=>     dt = [ -v + sqrt( v^2 + 0.2*ClosestNeighborDist*a) ] / a
	 */
	T dt = (std::sqrt(vi * vi + 0.2 * ai * this->closestNeighborDist[iBody]) - vi) / ai;

	if(dt == 0)
		dt = std::numeric_limits<T>::epsilon() / ai;

	return dt;
}
