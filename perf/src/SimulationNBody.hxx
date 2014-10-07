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

#include "SimulationNBody.h"

template <typename T>
SimulationNBody<T>::SimulationNBody(const unsigned long nBodies)
	: bodies(nBodies),
	  closestNeighborDist(NULL),
	  dt                 (std::numeric_limits<T>::infinity()),
	  dtConstant         (false)
{
	assert(nBodies > 0);
}

template <typename T>
SimulationNBody<T>::SimulationNBody(const std::string inputFileName)
	: bodies             (inputFileName),
	  closestNeighborDist(NULL),
	  dt                 (std::numeric_limits<T>::infinity()),
	  dtConstant         (false)
{
}

template <typename T>
SimulationNBody<T>::~SimulationNBody()
{
	if(this->accelerations.x)
		delete[] this->accelerations.x;
	if(this->accelerations.y)
		delete[] this->accelerations.y;
	if(this->accelerations.z)
		delete[] this->accelerations.z;

	if(this->closestNeighborDist)
		delete[] this->closestNeighborDist;
}

template <typename T>
Bodies<T>& SimulationNBody<T>::getBodies()
{
	return this->bodies;
}

template <typename T>
void SimulationNBody<T>::setDtConstant(T dtVal)
{
	this->dtConstant = true;
	this->dt = dtVal;
}

template <typename T>
void SimulationNBody<T>::setDtVariable()
{
	this->dtConstant = false;
	this->dt = std::numeric_limits<T>::infinity();
}

template <typename T>
T SimulationNBody<T>:: getDt()
{
	return this->dt;
}

template <typename T>
void SimulationNBody<T>::computeOneIteration()
{
	this->initIteration();
	this->computeBodiesAcceleration();
	if(!this->dtConstant)
		this->findTimeStep();
	this->bodies.updatePositionsAndVelocities(this->accelerations, this->dt);
}

template <typename T>
void SimulationNBody<T>::findTimeStep()
{
	if(!this->dtConstant)
	{
		this->dt = std::numeric_limits<T>::infinity();
		for(unsigned long iBody = 0; iBody < this->bodies.getN(); iBody++)
		{
			const T newDt = computeTimeStep(iBody);

			if(newDt < this->dt)
				this->dt = newDt;
		}
	}
}

template <typename T>
T SimulationNBody<T>::computeTimeStep(const unsigned long iBody)
{
	const T *velocitiesX = this->bodies.getVelocitiesX();
	const T *velocitiesY = this->bodies.getVelocitiesY();
	const T *velocitiesZ = this->bodies.getVelocitiesZ();

	// || velocity[iBody] ||
	const T v = std::sqrt((velocitiesX[iBody] * velocitiesX[iBody]) +
	                      (velocitiesY[iBody] * velocitiesY[iBody]) +
	                      (velocitiesZ[iBody] * velocitiesZ[iBody]));

	// || acceleration[iBody] ||
	const T a = std::sqrt((this->accelerations.x[iBody] * this->accelerations.x[iBody]) +
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
	T dt = (std::sqrt(v * v + 0.2 * a * this->closestNeighborDist[iBody]) - v) / a;

	if(dt == 0)
		dt = std::numeric_limits<T>::epsilon() / a;

	return dt;
}
