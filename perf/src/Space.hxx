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
unsigned long Space<T>::getNBodies()
{
	return this->nBodies;
}

template <typename T>
void Space<T>::allocateBuffers()
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
	this->allocateBuffers();

	srand(123);
	for(unsigned long iBody = 0; iBody < this->nBodies; iBody++)
	{
		this->masses[iBody] = ((rand() / (T) RAND_MAX) * 2000000) * G;

		this->positions.x[iBody] = (rand() / (T) RAND_MAX) * 800;
		this->positions.y[iBody] = (rand() / (T) RAND_MAX) * 600;
		this->positions.z[iBody] = 0;

		this->speeds.x[iBody] = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 0.02;
		this->speeds.y[iBody] = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 0.02;
		this->speeds.z[iBody] = 0;

		this->accelerations.x[iBody] = 0;
		this->accelerations.y[iBody] = 0;
		this->accelerations.z[iBody] = 0;

		this->closestNeighborLen[iBody] = std::numeric_limits<T>::infinity();
	}
}

template <typename T>
void Space<T>::initBodiesWithFile(const std::string inputFileName)
{
	std::ifstream bodiesFile;
	bodiesFile.open(inputFileName.c_str(), std::ios::in);

	if(!bodiesFile.is_open())
		std::cout << "Can't open \"" << inputFileName << "\" file (reading)." << std::endl;

	bool isOk = this->read(bodiesFile);
	bodiesFile.close();

	if(!isOk)
	{
		std::cout << "Something bad occurred during the reading of \"" << inputFileName
		          << "\" file... exiting." << std::endl;
		exit(-1);
	}
}

template <typename T>
void Space<T>::computeBodiesAcceleration()
{
	// flops ~= nBody^2 * 17
	for(unsigned long iBody = 0; iBody < this->nBodies; iBody++)
		// flops ~= nBody * 17
		for(unsigned long jBody = 0; jBody < this->nBodies; jBody++)
			if(iBody != jBody)
				this->computeAccelerationBetweenTwoBodies(iBody, jBody); // 17 flops
}

template <typename T>
void Space<T>::computeAccelerationBetweenTwoBodies(const unsigned long iBody, const unsigned long jBody)
{
	const T vecX   = this->positions.x[jBody] - this->positions.x[iBody]; // 1 flop
	const T vecY   = this->positions.y[jBody] - this->positions.y[iBody]; // 1 flop
	const T vecZ   = this->positions.z[jBody] - this->positions.z[iBody]; // 1 flop
	const T vecLen = sqrt((vecX * vecX) + (vecY * vecY) + (vecZ * vecZ)); // 5 flops

	if(vecLen == 0)
		std::cout << "Collision at {" << this->positions.x[jBody] << ", "
		                              << this->positions.y[jBody] << ", "
		                              << this->positions.z[jBody] << "}" << std::endl;
	assert(vecLen != 0);

	const T acc  = this->masses[jBody] / (vecLen * vecLen * vecLen); // 3 flops
	this->accelerations.x[iBody] += acc * vecX;                      // 2 flop
	this->accelerations.y[iBody] += acc * vecY;                      // 2 flop
	this->accelerations.z[iBody] += acc * vecZ;                      // 2 flop

	if(vecLen < this->closestNeighborLen[iBody])
		this->closestNeighborLen[iBody] = vecLen;
}

template <typename T>
void Space<T>::findTimeStep()
{
	this->dt = std::numeric_limits<T>::infinity();
	// flops = nBodies * 16
	for(unsigned long iBody = 0; iBody < this->nBodies; iBody++)
	{
		const T newDt = computeTimeStep(iBody); // 16 flops

		if(newDt < this->dt)
			this->dt = newDt;
	}
}

template <typename T>
T Space<T>::computeTimeStep(const unsigned long iBody)
{
	/* || lb.speed ||        */
	const T s = sqrt((this->speeds.x[iBody] * this->speeds.x[iBody]) +
	                 (this->speeds.y[iBody] * this->speeds.y[iBody]) +
	                 (this->speeds.z[iBody] * this->speeds.z[iBody])); // 5 flops

	/* || lb.acceleration || */
	const T a = sqrt((this->accelerations.x[iBody] * this->accelerations.x[iBody]) +
	                 (this->accelerations.y[iBody] * this->accelerations.y[iBody]) +
	                 (this->accelerations.z[iBody] * this->accelerations.z[iBody])); // 5 flops

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
	// flops = nBodies * 18
	for(unsigned long iBody = 0; iBody < this->nBodies; iBody++)
	{
		T accXMultDt = this->accelerations.x[iBody] * this->dt; // 1 flop
		T accYMultDt = this->accelerations.y[iBody] * this->dt; // 1 flop
		T accZMultDt = this->accelerations.z[iBody] * this->dt; // 1 flop

		this->positions.x[iBody] += (this->speeds.x[iBody] + accXMultDt * 0.5) * this->dt; // 4 flops
		this->positions.y[iBody] += (this->speeds.y[iBody] + accYMultDt * 0.5) * this->dt; // 4 flops
		this->positions.z[iBody] += (this->speeds.z[iBody] + accZMultDt * 0.5) * this->dt; // 4 flops;

		this->speeds.x[iBody] += accXMultDt; // 1 flop
		this->speeds.y[iBody] += accYMultDt; // 1 flop
		this->speeds.z[iBody] += accZMultDt; // 1 flop

		this->accelerations.x[iBody] = 0;
		this->accelerations.y[iBody] = 0;
		this->accelerations.z[iBody] = 0;

		this->closestNeighborLen[iBody] = std::numeric_limits<T>::infinity();
	}
}

template <typename T>
bool Space<T>::read(std::istream& stream)
{
	this->nBodies = 0;
	stream >> this->nBodies;

	if(this->nBodies)
		this->allocateBuffers();
	else
		return false;

	for(unsigned long iBody = 0; iBody < this->nBodies; iBody++)
	{
		stream >> this->masses[iBody];
		this->masses[iBody] *= G;

		stream >> this->positions.x[iBody];
		stream >> this->positions.y[iBody];
		stream >> this->positions.z[iBody];

		stream >> this->speeds.x[iBody];
		stream >> this->speeds.y[iBody];
		stream >> this->speeds.z[iBody];

		this->accelerations.x[iBody] = 0;
		this->accelerations.y[iBody] = 0;
		this->accelerations.z[iBody] = 0;

		this->closestNeighborLen[iBody] = std::numeric_limits<T>::infinity();

		if(!stream.good())
			return false;
	}

	return true;
}

template <typename T>
void Space<T>::write(std::ostream& stream)
{
	stream << this->nBodies << std::endl;

	for(unsigned long iBody = 0; iBody < this->nBodies; iBody++)
		stream << this->masses     [iBody] / G << " "
		       << this->positions.x[iBody]     << " "
		       << this->positions.y[iBody]     << " "
		       << this->positions.z[iBody]     << " "
		       << this->speeds.x   [iBody]     << " "
		       << this->speeds.y   [iBody]     << " "
		       << this->speeds.z   [iBody]     << std::endl;
}

template <typename T>
void Space<T>::writeIntoFile(const std::string outputFileName)
{
	std::fstream bodiesFile(outputFileName.c_str(), std::ios_base::out);
	if(!bodiesFile.is_open())
	{
		std::cout << "Can't open \"" << outputFileName << "\" file (writing). Exiting..." << std::endl;
		exit(-1);
	}

	this->write(bodiesFile);

	bodiesFile.close();
}

template <typename T>
std::ostream& operator<<(std::ostream &o, const Space<T>& s)
{
	s.write(o);
	return o;
}
