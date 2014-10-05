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

#include "Space.h"

template <typename T>
Space<T>::Space(const unsigned long nBodies)
	: nBodies            (nBodies),
	  masses             (NULL),
	  radiuses           (NULL),
	  closestNeighborDist(NULL),
	  dt                 (std::numeric_limits<T>::infinity()),
	  dtConstant         (false)
{
	assert(nBodies > 0);
	this->initBodiesRandomly();
}

template <typename T>
Space<T>::Space(const std::string inputFileName)
	: nBodies            (0),
	  masses             (NULL),
	  radiuses           (NULL),
	  closestNeighborDist(NULL),
	  dt                 (std::numeric_limits<T>::infinity()),
	  dtConstant         (false)
{
	this->initBodiesFromFile(inputFileName);
}

template <typename T>
void Space<T>::allocateBuffers()
{
	this->masses = new T[this->nBodies];

	this->radiuses = new T[this->nBodies];

	this->positions.x = new T[this->nBodies];
	this->positions.y = new T[this->nBodies];
	this->positions.z = new T[this->nBodies];

	this->speeds.x = new T[this->nBodies];
	this->speeds.y = new T[this->nBodies];
	this->speeds.z = new T[this->nBodies];

	this->accelerations.x = new T[this->nBodies];
	this->accelerations.y = new T[this->nBodies];
	this->accelerations.z = new T[this->nBodies];

	this->closestNeighborDist = new T[this->nBodies];
}

template <typename T>
Space<T>::~Space() {
	if(this->masses)
		delete[] this->masses;

	if(this->radiuses)
		delete[] this->radiuses;

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

	if(this->closestNeighborDist)
		delete[] this->closestNeighborDist;
}

template <typename T>
unsigned long Space<T>::getNBodies()
{
	return this->nBodies;
}

template <typename T>
void Space<T>::setDtConstant(T dtVal)
{
	this->dtConstant = true;
	this->dt = dtVal;
}

template <typename T>
void Space<T>::setDtVariable()
{
	this->dtConstant = false;
	this->dt = std::numeric_limits<T>::infinity();
}

template <typename T>
T Space<T>:: getDt()
{
	return this->dt;
}

template <typename T>
void Space<T>::setBody(const unsigned long &iBody,
                       const T &mass, const T &radius,
                       const T &posX, const T &posY, const T &posZ,
                       const T &speedX, const T &speedY, const T &speedZ)
{
	this->masses[iBody] = mass;

	this->radiuses[iBody] = radius;

	this->positions.x[iBody] = posX;
	this->positions.y[iBody] = posY;
	this->positions.z[iBody] = posZ;

	this->speeds.x[iBody] = speedX;
	this->speeds.y[iBody] = speedY;
	this->speeds.z[iBody] = speedZ;

	this->accelerations.x[iBody] = 0;
	this->accelerations.y[iBody] = 0;
	this->accelerations.z[iBody] = 0;

	this->closestNeighborDist[iBody] = std::numeric_limits<T>::infinity();
}

template <typename T>
void Space<T>::initBodiesRandomly()
{
	this->allocateBuffers();

	srand(123);
	for(unsigned long iBody = 0; iBody < this->nBodies; iBody++)
	{
		T mass, radius, posX, posY, posZ, speedX, speedY, speedZ;

		mass = ((rand() / (T) RAND_MAX) * 5.0e21);

		radius = mass * 0.6e-15;

		posX = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * (5.0e8 * 1.33);
		posY = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 5.0e8;
		posZ = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 5.0e8 -10.0e8;

		speedX = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
		speedY = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
		speedZ = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;

		this->setBody(iBody, mass, radius, posX, posY, posZ, speedX, speedY, speedZ);
	}
}

template <typename T>
void Space<T>::initBodiesFromFile(const std::string inputFileName)
{
	std::ifstream bodiesFile;
	bodiesFile.open(inputFileName.c_str(), std::ios::in);

	if(!bodiesFile.is_open())
	{
		std::cout << "Can't open \"" << inputFileName << "\" file (reading)." << std::endl;
		exit(-1);
	}

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
#pragma omp parallel for schedule(runtime)
	for(unsigned long iBody = 0; iBody < this->nBodies; iBody++)
		for(unsigned long jBody = 0; jBody < this->nBodies; jBody++)
			if(iBody != jBody)
				//this->computeAccelerationBetweenTwoBodies(iBody, jBody);
				this->computeAccelerationBetweenTwoBodiesNaive(iBody, jBody);
}

template <typename T>
void Space<T>::computeBodiesAccelerationV2()
{
#pragma omp parallel for schedule(runtime)
	for(unsigned long iBody = 0; iBody < this->nBodies; iBody++)
		for(unsigned long jBody = iBody +1; jBody < this->nBodies; jBody++)
			//this->computeAccelerationBetweenTwoBodies(iBody, jBody);
			this->computeAccelerationBetweenTwoBodiesNaiveV2(iBody, jBody);
}

/* 
	AI  = (23 * blockSize * nBodies * nBlocks)  / ((4 * blockSize + 7 * nBodies) * nBlocks) <=>
	AI  = (23 * blockSize * nBodies)            /  (4 * blockSize + 7 * nBodies)            <=>
	AI  = (23 * blockSize * nBlock * blockSize) /  (4 * blockSize + 7 * nBlock * blockSize) <=>
	AI  = (23 * nBlock * blockSize²)            / ((4 + 7 * nBlock) * blockSize)            <=>
	AI  = (23 * nBlock * blockSize)             /  (4 + 7 * nBlock)                         <=>
	AI ~= (23 * nBlock * blockSize)             /      (7 * nBlock)                         <=>
	AI ~= (23 * blockSize)                      /       7
	-------------------------------------------------------------------------------------------
	OI  = AI                                    /      sizeof(T)                            <=>
	OI  = (23 * blockSize)                      / (7 * sizeof(T))
*/
template <typename T>
void Space<T>::computeBodiesAccelerationCB()
{
	unsigned long blockSize = 512;
	// flops  = 23 * blockSize * nBodies      * nBlocks
	// memops = (4 * blockSize + 7 * nBodies) * nBlocks
	for(unsigned long jOff = 0; jOff < this->nBodies; jOff += blockSize)
	{
		blockSize = std::min(blockSize, this->nBodies - jOff);
		// flops  = 23 * blockSize * nBodies
		// memops =  4 * blockSize + 7 * nBodies 
#pragma omp parallel for schedule(runtime)
		for(unsigned long iBody = 0; iBody < this->nBodies; iBody++)
			// flops  = 23 * blockSize
			// memops =  4 * blockSize + 7 
			for(unsigned long jBody = jOff; jBody < jOff + blockSize; jBody++)
				if(iBody != jBody)
					this->computeAccelerationBetweenTwoBodies(iBody, jBody);
					//this->computeAccelerationBetweenTwoBodiesNaive(iBody, jBody);
	}
}

// 18 flops
template <typename T>
void Space<T>::computeAccelerationBetweenTwoBodies(const unsigned long iBody, const unsigned long jBody)
{
	assert(iBody != jBody);

	const T diffPosX = this->positions.x[jBody] - this->positions.x[iBody]; // 1 flop
	const T diffPosY = this->positions.y[jBody] - this->positions.y[iBody]; // 1 flop
	const T diffPosZ = this->positions.z[jBody] - this->positions.z[iBody]; // 1 flop
	const T squareDist = (diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ); // 5 flops
	const T dist = std::sqrt(squareDist); // 1 flop
	assert(dist != 0);

	const T acc = G * this->masses[jBody] / (squareDist * dist); // 3 flops
	this->accelerations.x[iBody] += acc * diffPosX; // 2 flop
	this->accelerations.y[iBody] += acc * diffPosY; // 2 flop
	this->accelerations.z[iBody] += acc * diffPosZ; // 2 flop

	if(!this->dtConstant)
		if(dist < this->closestNeighborDist[iBody])
			this->closestNeighborDist[iBody] = dist;
}

// 23 flops
template <typename T>
void Space<T>::computeAccelerationBetweenTwoBodiesNaive(const unsigned long iBody, const unsigned long jBody)
{
	assert(iBody != jBody);

	const T diffPosX  = this->positions.x[jBody] - this->positions.x[iBody]; // 1 flop
	const T diffPosY  = this->positions.y[jBody] - this->positions.y[iBody]; // 1 flop
	const T diffPosZ  = this->positions.z[jBody] - this->positions.z[iBody]; // 1 flop

	// compute distance between iBody and jBody: Dij
	const T dist = std::sqrt((diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ)); // 6 flops

	// we cannot divide by 0
	if(dist == 0)
	{
		std::cout << "Collision at {" << this->positions.x[jBody] << ", "
		                              << this->positions.y[jBody] << ", "
		                              << this->positions.z[jBody] << "}" << std::endl;
		assert(dist != 0);
	}

	// compute the force value between iBody and jBody: || F || = G.mi.mj / Dij²
	const T force = G * this->masses[iBody] * this->masses[jBody] / (dist * dist); // 4 flops

	// compute the acceleration value: || a || = || F || / mi
	const T acc = force / this->masses[iBody]; // 1 flop

	// normalize and add acceleration value into acceleration vector: a += || a ||.u
	this->accelerations.x[iBody] += acc * (diffPosX / dist); // 3 flops
	this->accelerations.y[iBody] += acc * (diffPosY / dist); // 3 flops
	this->accelerations.z[iBody] += acc * (diffPosZ / dist); // 3 flops

	if(!this->dtConstant)
		if(dist < this->closestNeighborDist[iBody])
			this->closestNeighborDist[iBody] = dist;
}

// 32 flops
template <typename T>
void Space<T>::computeAccelerationBetweenTwoBodiesNaiveV2(const unsigned long iBody, const unsigned long jBody)
{
	assert(iBody != jBody);

	const T diffPosX  = this->positions.x[jBody] - this->positions.x[iBody]; // 1 flop
	const T diffPosY  = this->positions.y[jBody] - this->positions.y[iBody]; // 1 flop
	const T diffPosZ  = this->positions.z[jBody] - this->positions.z[iBody]; // 1 flop

	// compute distance between iBody and jBody: Dij
	const T dist = std::sqrt((diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ)); // 6 flops

	// we cannot divide by 0
	if(dist == 0)
	{
		std::cout << "Collision at {" << this->positions.x[jBody] << ", "
		                              << this->positions.y[jBody] << ", "
		                              << this->positions.z[jBody] << "}" << std::endl;
		assert(dist != 0);
	}

	// compute the force value between iBody and jBody: || F || = G.mi.mj / Dij²
	const T force = G * this->masses[iBody] * this->masses[jBody] / (dist * dist); // 4 flops

	// compute the acceleration value: || a || = || F || / mi
	const T acc = force / this->masses[iBody]; // 1 flop

	// normalize and add acceleration value into acceleration vector: a += || a ||.u
	this->accelerations.x[iBody] += acc * (diffPosX / dist); // 3 flops
	this->accelerations.y[iBody] += acc * (diffPosY / dist); // 3 flops
	this->accelerations.z[iBody] += acc * (diffPosZ / dist); // 3 flops

	this->accelerations.x[jBody] += acc * ((-diffPosX) / dist); // 3 flops
	this->accelerations.y[jBody] += acc * ((-diffPosY) / dist); // 3 flops
	this->accelerations.z[jBody] += acc * ((-diffPosZ) / dist); // 3 flops

	if(!this->dtConstant)
	{
		if(dist < this->closestNeighborDist[iBody])
			this->closestNeighborDist[iBody] = dist;
		if(dist < this->closestNeighborDist[jBody])
			this->closestNeighborDist[jBody] = dist;
	}
}

template <typename T>
void Space<T>::findTimeStep()
{
	if(!this->dtConstant)
	{
		this->dt = std::numeric_limits<T>::infinity();
		for(unsigned long iBody = 0; iBody < this->nBodies; iBody++)
		{
			const T newDt = computeTimeStep(iBody);

			if(newDt < this->dt)
				this->dt = newDt;
		}
	}
}

template <typename T>
T Space<T>::computeTimeStep(const unsigned long iBody)
{
	// || lb.speed ||
	const T s = std::sqrt((this->speeds.x[iBody] * this->speeds.x[iBody]) +
	                      (this->speeds.y[iBody] * this->speeds.y[iBody]) +
	                      (this->speeds.z[iBody] * this->speeds.z[iBody]));

	// || lb.acceleration ||
	const T a = std::sqrt((this->accelerations.x[iBody] * this->accelerations.x[iBody]) +
	                      (this->accelerations.y[iBody] * this->accelerations.y[iBody]) +
	                      (this->accelerations.z[iBody] * this->accelerations.z[iBody]));

	/*
	 * compute dt
	 * solve:  (a/2)*dt^2 + s*dt + (-0.1)*ClosestNeighborDist = 0
	 * <=>     dt = [ (-s) +/-  sqrt( s^2 - 4 * (a/2) * (-0.1)*ClosestNeighborDist ) ] / [ 2 (a/2) ]
	 *
	 * dt should be positive (+/- becomes + because result of sqrt is positive)
	 * <=>     dt = [ -s + sqrt( s^2 + 0.2*ClosestNeighborDist*a) ] / a
	 */
	T dt = (std::sqrt(s * s + 0.2 * a * this->closestNeighborDist[iBody]) - s) / a;

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
		T mass, radius, posX, posY, posZ, speedX, speedY, speedZ;

		mass = this->masses[iBody];

		radius = this->radiuses[iBody];

		T accXMultDt = this->accelerations.x[iBody] * this->dt;
		T accYMultDt = this->accelerations.y[iBody] * this->dt;
		T accZMultDt = this->accelerations.z[iBody] * this->dt;

		posX = this->positions.x[iBody] + (this->speeds.x[iBody] + accXMultDt * 0.5) * this->dt;
		posY = this->positions.y[iBody] + (this->speeds.y[iBody] + accYMultDt * 0.5) * this->dt;
		posZ = this->positions.z[iBody] + (this->speeds.z[iBody] + accZMultDt * 0.5) * this->dt;

		speedX = this->speeds.x[iBody] + accXMultDt;
		speedY = this->speeds.y[iBody] + accYMultDt;
		speedZ = this->speeds.z[iBody] + accZMultDt;

		this->setBody(iBody, mass, radius, posX, posY, posZ, speedX, speedY, speedZ);
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
		T mass, radius, posX, posY, posZ, speedX, speedY, speedZ;

		stream >> mass;

		stream >> radius;

		stream >> posX;
		stream >> posY;
		stream >> posZ;

		stream >> speedX;
		stream >> speedY;
		stream >> speedZ;

		this->setBody(iBody, mass, radius, posX, posY, posZ, speedX, speedY, speedZ);

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
		stream << this->masses     [iBody] << " "
		       << this->radiuses   [iBody] << " "
		       << this->positions.x[iBody] << " "
		       << this->positions.y[iBody] << " "
		       << this->positions.z[iBody] << " "
		       << this->speeds.x   [iBody] << " "
		       << this->speeds.y   [iBody] << " "
		       << this->speeds.z   [iBody] << std::endl;
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
