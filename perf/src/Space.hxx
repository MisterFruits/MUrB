/*
 * Do not remove.
 * MPI/OpenMP training courses
 * Adrien Cassagne, ASA - CINES, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#include <cmath>
#include <limits>
#include <string>
#include <cassert>
#include <fstream>
#include <iostream>

#include "Space.h"

template <typename T>
Space<T>::Space(const unsigned long nBodies)
	: nBodies(nBodies), masses(NULL), closestNeighborLen(NULL), dt(std::numeric_limits<T>::infinity())
{
	assert(nBodies > 0);
	this->initBodiesRandomly();
}

template <typename T>
Space<T>::Space(const std::string inputFileName)
	: nBodies(0), masses(NULL), closestNeighborLen(NULL), dt(std::numeric_limits<T>::infinity())
{
	this->initBodiesWithFile(inputFileName);
}

template <typename T>
Space<T>::~Space() {
	if(this->masses)
		delete[] this->masses;

	if(this->positions.x)
		delete[] this->positions.x;
	if(this->positions.y)
		delete[] this->positions.y;
	if(this->positions.z)
		delete[] this->positions.z;

	if(this->speeds.x)
		delete[] this->speeds.x;
	if(this->speeds.y)
		delete[] this->speeds.y;
	if(this->speeds.z)
		delete[] this->speeds.z;

	if(this->accelerations.x)
		delete[] this->accelerations.x;
	if(this->accelerations.y)
		delete[] this->accelerations.y;
	if(this->accelerations.z)
		delete[] this->accelerations.z;

	if(this->closestNeighborLen)
		delete[] this->closestNeighborLen;
}

template <typename T>
void Space<T>::initBuffers()
{
	this->masses = new T[this->nBodies];

	this->positions.x = new T[this->nBodies];
	this->positions.y = new T[this->nBodies];
	this->positions.z = new T[this->nBodies];

	this->speeds.x = new T[this->nBodies];
	this->speeds.y = new T[this->nBodies];
	this->speeds.z = new T[this->nBodies];

	this->accelerations.x = new T[this->nBodies];
	this->accelerations.y = new T[this->nBodies];
	this->accelerations.z = new T[this->nBodies];

	this->closestNeighborLen = new T[this->nBodies];
}

template <typename T>
void Space<T>::initBodiesRandomly()
{
	this->initBuffers();

	srand(123);
	for(unsigned long iBody = 0; iBody < this->nBodies; iBody++)
	{
		const double mass   = (rand() / (double) RAND_MAX) * 2000000;
		const double posX   = (rand() / (double) RAND_MAX) * 800;
		const double posY   = (rand() / (double) RAND_MAX) * 600;
		const double speedX = ((rand() - RAND_MAX/2) / (double) (RAND_MAX/2)) * 0.02;
		const double speedY = ((rand() - RAND_MAX/2) / (double) (RAND_MAX/2)) * 0.02;

		this->initBody(iBody, mass, posX, posY, 0, speedX, speedY, 0);
	}
}

template <typename T>
void Space<T>::initBodiesWithFile(const std::string inputFileName)
{
	std::ifstream bodiesFile;
	bodiesFile.open(inputFileName.c_str(), std::ios::in);

	if(!bodiesFile.is_open())
		std::cout << "Can't open \"" << inputFileName << "\" file (reading)." << std::endl;

	bodiesFile >> this->nBodies;

	if(this->nBodies)
		this->initBuffers();
	else
	{
		std::cout << "Empty \"" << inputFileName << "\" file (reading)... exiting." << std::endl;
		exit(-1);
	}

	double mass, posX, posY, speedX, speedY;
	for(unsigned long iBody = 0; iBody < this->nBodies; iBody++)
	{
		bodiesFile >> mass;
		bodiesFile >> posX;
		bodiesFile >> posY;
		bodiesFile >> speedX;
		bodiesFile >> speedY;

		this->initBody(iBody, mass, posX, posY, 0, speedX, speedY, 0);

		if(!bodiesFile.good())
		{
			std::cout << "Something bad occurred during the reading of \"" << inputFileName
			          << "\" file... exiting." << std::endl;
			bodiesFile.close();
			exit(-1);
		}
	}

	bodiesFile.close();
}

template <typename T>
void Space<T>::initBody(unsigned long iBody, T mass, T posX, T posY, T posZ, T speedX, T speedY, T speedZ)
{
	this->masses[iBody] = mass * G;

	this->positions.x[iBody] = posX;
	this->positions.y[iBody] = posY;
	this->positions.z[iBody] = posZ;

	this->speeds.x[iBody] = speedX;
	this->speeds.y[iBody] = speedY;
	this->speeds.z[iBody] = speedZ;

	this->accelerations.x[iBody] = 0;
	this->accelerations.y[iBody] = 0;
	this->accelerations.z[iBody] = 0;

	this->closestNeighborLen[iBody] = std::numeric_limits<double>::infinity();
}

template <typename T>
void Space<T>::computeBodiesAcceleration()
{
	// flops ~= nBody^2 * 12
	for(unsigned long iBody = 0; iBody < this->nBodies; iBody++)
		// flops ~= nBody * 12
		for(unsigned long jBody = 0; jBody < this->nBodies; jBody++)
			if(iBody != jBody)
				this->computeAccelerationBetweenTwoBodies(iBody, jBody); // 12 flops
}

template <typename T>
void Space<T>::computeAccelerationBetweenTwoBodies(const unsigned long iBody, const unsigned long jBody)
{
	const T vecX       = this->positions.x[jBody] - this->positions.x[iBody]; // 1 flop
	const T vecY       = this->positions.y[jBody] - this->positions.y[iBody]; // 1 flop
	const T vecLen     = (vecX * vecX) + (vecY * vecY);                       // 3 flops
	const T sqrtVecLen = sqrt(vecLen);                                        

	if(vecLen == 0)
		std::cout << "Collision at {" << this->positions.x[jBody] << ", "
		                              << this->positions.y[jBody] << "}" << std::endl;
	assert(vecLen != 0);

	const T acc  = this->masses[jBody] / (vecLen * sqrtVecLen); // 2 flops
	const T accX = acc * vecX;                                  // 1 flop
	const T accY = acc * vecY;                                  // 1 flop

	this->accelerations.x[iBody] += accX; // 1 flop
	this->accelerations.y[iBody] += accY; // 1 flop

	if(sqrtVecLen < this->closestNeighborLen[iBody])
		this->closestNeighborLen[iBody] = sqrtVecLen;
}

template <typename T>
void Space<T>::findTimeStep()
{
	this->dt = std::numeric_limits<double>::infinity();
	// flops = nBody * 12
	for(unsigned long iBody = 0; iBody < this->nBodies; iBody++)
	{
		const T newDt = computeTimeStep(iBody); // 12 flops

		if(newDt < this->dt)
			this->dt = newDt;
	}
}

template <typename T>
T Space<T>::computeTimeStep(const unsigned long iBody)
{
	/* || lb.speed ||        */
	const T s = sqrt((this->speeds.x[iBody] * this->speeds.x[iBody]) +
	                 (this->speeds.y[iBody] * this->speeds.y[iBody])); // 3 flops

	/* || lb.acceleration || */
	const T a = sqrt((this->accelerations.x[iBody] * this->accelerations.x[iBody]) +
	                 (this->accelerations.y[iBody] * this->accelerations.y[iBody])); // 3 flops

	/*
	 * compute dt
	 * solve:  (a/2)*dt^2 + s*dt + (-0.1)*ClosestNeighborLen = 0
	 * <=>     dt = [ (-s) +/-  sqrt( s^2 - 4 * (a/2) * (-0.1)*ClosestNeighborLen ) ] / [ 2 (a/2) ]
	 *
	 * dt should be positive (+/- becomes + because result of sqrt is positive)
	 * <=>     dt = [ -s + sqrt( s^2 + 0.2*ClosestNeighborLen*a) ] / a
	 */

	T dt = (sqrt(s * s + 0.2 * a * this->closestNeighborLen[iBody]) - s) / a; // 6 flops

	if(dt == 0)
		dt = std::numeric_limits<T>::epsilon() / a;

	return dt;
}

template <typename T>
void Space<T>::updateBodiesPositionAndSpeed()
{
	// flops = nBody * 12
	for(unsigned long iBody = 0; iBody < this->nBodies; iBody++)
	{
		T accXMultDt = this->accelerations.x[iBody] * this->dt; // 1 flop
		T accYMultDt = this->accelerations.y[iBody] * this->dt; // 1 flop

		this->positions.x[iBody] += (this->speeds.x[iBody] + accXMultDt / 2.0) * this->dt; // 4 flops
		this->positions.y[iBody] += (this->speeds.y[iBody] + accYMultDt / 2.0) * this->dt; // 4 flops

		this->speeds.x[iBody] += accXMultDt; // 1 flop
		this->speeds.y[iBody] += accYMultDt; // 1 flop

		this->accelerations.x[iBody] = 0;
		this->accelerations.y[iBody] = 0;

		this->closestNeighborLen[iBody] = std::numeric_limits<double>::infinity();
	}
}

template <typename T>
void Space<T>::write(std::ostream& stream)
{
	assert(this->nBodies > 0);

	stream << this->nBodies << std::endl;

	for(unsigned long iBody = 0; iBody < this->nBodies; iBody++)
		stream << this->masses     [iBody] / G << " "
		       << this->positions.x[iBody]     << " "
		       << this->positions.y[iBody]     << " "
		       << this->speeds.x   [iBody]     << " "
		       << this->speeds.y   [iBody]     << std::endl;
}

template <typename T>
void Space<T>::writeIntoFile(const std::string outputFileName)
{
	assert(this->nBodies > 0);

	std::fstream bodiesFile(outputFileName.c_str(), std::ios_base::out);
	if(!bodiesFile.is_open())
	{
		std::cout << "Can't open \"" << outputFileName << "\" file (writing). Exiting..." << std::endl;
		exit(-1);
	}

	this->write(bodiesFile);

	bodiesFile.close();
}
