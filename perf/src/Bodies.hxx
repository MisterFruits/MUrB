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

#include "utils/myIntrinsicsPlusPlus.h"

#include "Bodies.h"

template <typename T>
Bodies<T>::Bodies(const unsigned long n)
	: n             (n),
	  masses        (NULL),
	  radiuses      (NULL),
	  nVecs         (ceil((T) n / (T) mipp::vectorSize<T>())),
	  allocatedBytes(0)
{
	assert(n > 0);
	this->initRandomly();
}

template <typename T>
Bodies<T>::Bodies(const std::string inputFileName)
	: n             (0),
	  masses        (NULL),
	  radiuses      (NULL),
	  nVecs         (0),
	  allocatedBytes(0)
{
	this->initFromFile(inputFileName);
}

template <typename T>
Bodies<T>::Bodies(const Bodies<T>& bodies)
	: n             (bodies.n),
	  masses        (bodies.masses),
	  radiuses      (bodies.radiuses),
	  nVecs         (bodies.nVecs),
	  allocatedBytes(bodies.allocatedBytes)
{
	this->positions.x = bodies.positions.x;
	this->positions.y = bodies.positions.y;
	this->positions.z = bodies.positions.z;

	this->velocities.x = bodies.velocities.x;
	this->velocities.y = bodies.velocities.y;
	this->velocities.z = bodies.velocities.z;
}

template <typename T>
Bodies<T>& Bodies<T>::operator=(const Bodies<T>& bodies)
{
	return Bodies<T>(bodies);
}

template <typename T>
void Bodies<T>::allocateBuffers()
{
	const unsigned long padding = (this->nVecs * mipp::vectorSize<T>()) - this->n;

#ifdef __ARM_NEON__
	this->masses = new T[this->n + padding];

	this->radiuses = new T[this->n + padding];

	this->positions.x = new T[this->n + padding];
	this->positions.y = new T[this->n + padding];
	this->positions.z = new T[this->n + padding];

	this->velocities.x = new T[this->n + padding];
	this->velocities.y = new T[this->n + padding];
	this->velocities.z = new T[this->n + padding];
#else
	this->masses = (T*)_mm_malloc((this->n + padding) * sizeof(T), mipp::RequiredAlignement);

	this->radiuses = (T*)_mm_malloc((this->n + padding) * sizeof(T), mipp::RequiredAlignement);

	this->positions.x = (T*)_mm_malloc((this->n + padding) * sizeof(T), mipp::RequiredAlignement);
	this->positions.y = (T*)_mm_malloc((this->n + padding) * sizeof(T), mipp::RequiredAlignement);
	this->positions.z = (T*)_mm_malloc((this->n + padding) * sizeof(T), mipp::RequiredAlignement);

	this->velocities.x = (T*)_mm_malloc((this->n + padding) * sizeof(T), mipp::RequiredAlignement);
	this->velocities.y = (T*)_mm_malloc((this->n + padding) * sizeof(T), mipp::RequiredAlignement);
	this->velocities.z = (T*)_mm_malloc((this->n + padding) * sizeof(T), mipp::RequiredAlignement);
#endif

	this->allocatedBytes = (this->n + padding) * sizeof(T) * 8;
}

template <typename T>
Bodies<T>::~Bodies() {
#ifdef __ARM_NEON__
	if(this->masses != nullptr) {
		delete[] this->masses;
		this->masses = nullptr;
	}

	if(this->radiuses != nullptr) {
		delete[] this->radiuses;
		this->radiuses = nullptr;
	}

	if(this->positions.x != nullptr) {
		delete[] this->positions.x;
		this->positions.x = nullptr;
	}
	if(this->positions.y != nullptr) {
		delete[] this->positions.y;
		this->positions.y = nullptr;
	}
	if(this->positions.z != nullptr) {
		delete[] this->positions.z;
		this->positions.z = nullptr;
	}

	if(this->velocities.x != nullptr) {
		delete[] this->velocities.x;
		this->velocities.x = nullptr;
	}
	if(this->velocities.y != nullptr) {
		delete[] this->velocities.y;
		this->velocities.y = nullptr;
	}
	if(this->velocities.z != nullptr) {
		delete[] this->velocities.z;
		this->velocities.z = nullptr;
	}
#else
	if(this->masses != nullptr) {
		_mm_free(this->masses);
		this->masses = nullptr;
	}

	if(this->radiuses != nullptr) {
		_mm_free(this->radiuses);
		this->radiuses = nullptr;
	}

	if(this->positions.x != nullptr) {
		_mm_free(this->positions.x);
		this->positions.x = nullptr;
	}
	if(this->positions.y != nullptr) {
		_mm_free(this->positions.y);
		this->positions.y = nullptr;
	}
	if(this->positions.z != nullptr) {
		_mm_free(this->positions.z);
		this->positions.z = nullptr;
	}

	if(this->velocities.x != nullptr) {
		_mm_free(this->velocities.x);
		this->velocities.x = nullptr;
	}
	if(this->velocities.y != nullptr) {
		_mm_free(this->velocities.y);
		this->velocities.y = nullptr;
	}
	if(this->velocities.z != nullptr) {
		_mm_free(this->velocities.z);
		this->velocities.z = nullptr;
	}
#endif
}

template <typename T>
const unsigned long& Bodies<T>::getN()
{
	return const_cast<const unsigned long&>(this->n);
}

template <typename T>
const unsigned long& Bodies<T>::getNVecs()
{
	return const_cast<const unsigned long&>(this->nVecs);
}

template <typename T>
const T* Bodies<T>::getMasses()
{
	return const_cast<const T*>(this->masses);
}

template <typename T>
const T* Bodies<T>::getRadiuses()
{
	return const_cast<const T*>(this->radiuses);
}

template <typename T>
const T* Bodies<T>::getPositionsX()
{
	return const_cast<const T*>(this->positions.x);
}

template <typename T>
const T* Bodies<T>::getPositionsY()
{
	return const_cast<const T*>(this->positions.y);
}

template <typename T>
const T* Bodies<T>::getPositionsZ()
{
	return const_cast<const T*>(this->positions.z);
}

template <typename T>
const T* Bodies<T>::getVelocitiesX()
{
	return const_cast<const T*>(this->velocities.x);
}

template <typename T>
const T* Bodies<T>::getVelocitiesY()
{
	return const_cast<const T*>(this->velocities.y);
}

template <typename T>
const T* Bodies<T>::getVelocitiesZ()
{
	return const_cast<const T*>(this->velocities.z);
}

template <typename T>
const float& Bodies<T>::getAllocatedBytes()
{
	return this->allocatedBytes;
}

template <typename T>
void Bodies<T>::setBody(const unsigned long &iBody,
                        const T &mass, const T &radius,
                        const T &posX, const T &posY, const T &posZ,
                        const T &velocityX, const T &velocityY, const T &velocityZ)
{
	this->masses[iBody] = mass;

	this->radiuses[iBody] = radius;

	this->positions.x[iBody] = posX;
	this->positions.y[iBody] = posY;
	this->positions.z[iBody] = posZ;

	this->velocities.x[iBody] = velocityX;
	this->velocities.y[iBody] = velocityY;
	this->velocities.z[iBody] = velocityZ;
}

template <typename T>
void Bodies<T>::initRandomly()
{
	this->allocateBuffers();

	srand(123);
	for(unsigned long iBody = 0; iBody < this->n; iBody++)
	{
		T mass, radius, posX, posY, posZ, velocityX, velocityY, velocityZ;

		mass = ((rand() / (T) RAND_MAX) * 5.0e21);

		radius = mass * 0.6e-15;

		posX = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * (5.0e8 * 1.33);
		posY = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 5.0e8;
		posZ = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 5.0e8 -10.0e8;

		velocityX = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
		velocityY = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
		velocityZ = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;

		this->setBody(iBody, mass, radius, posX, posY, posZ, velocityX, velocityY, velocityZ);
	}

	// fill the bodies in the padding zone
	const unsigned long padding = (this->nVecs * mipp::vectorSize<T>()) - this->n;
	for(unsigned long iBody = this->n; iBody < this->n + padding; iBody++)
	{
		T posX, posY, posZ, velocityX, velocityY, velocityZ;

		posX = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * (5.0e8 * 1.33);
		posY = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 5.0e8;
		posZ = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 5.0e8 -10.0e8;

		velocityX = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
		velocityY = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
		velocityZ = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;

		this->setBody(iBody, 0, 0, posX, posY, posZ, velocityX, velocityY, velocityZ);
	}
}

template <typename T>
void Bodies<T>::initFromFile(const std::string inputFileName)
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
void Bodies<T>::updatePositionsAndVelocities(const vector3<T> &accelerations, T &dt)
{
	// flops = n * 18
	for(unsigned long iBody = 0; iBody < this->n; iBody++)
	{
		T mass, radius, posX, posY, posZ, velocityX, velocityY, velocityZ;

		mass = this->masses[iBody];

		radius = this->radiuses[iBody];

		T accXMultDt = accelerations.x[iBody] * dt;
		T accYMultDt = accelerations.y[iBody] * dt;
		T accZMultDt = accelerations.z[iBody] * dt;

		posX = this->positions.x[iBody] + (this->velocities.x[iBody] + accXMultDt * 0.5) * dt;
		posY = this->positions.y[iBody] + (this->velocities.y[iBody] + accYMultDt * 0.5) * dt;
		posZ = this->positions.z[iBody] + (this->velocities.z[iBody] + accZMultDt * 0.5) * dt;

		velocityX = this->velocities.x[iBody] + accXMultDt;
		velocityY = this->velocities.y[iBody] + accYMultDt;
		velocityZ = this->velocities.z[iBody] + accZMultDt;

		this->setBody(iBody, mass, radius, posX, posY, posZ, velocityX, velocityY, velocityZ);
	}
}

template <typename T>
bool Bodies<T>::read(std::istream& stream)
{
	this->n = 0;
	stream >> this->n;

	this->nVecs = ceil((T) this->n / (T) mipp::vectorSize<T>());

	if(this->n)
		this->allocateBuffers();
	else
		return false;

	for(unsigned long iBody = 0; iBody < this->n; iBody++)
	{
		T mass, radius, posX, posY, posZ, velocityX, velocityY, velocityZ;

		stream >> mass;

		stream >> radius;

		stream >> posX;
		stream >> posY;
		stream >> posZ;

		stream >> velocityX;
		stream >> velocityY;
		stream >> velocityZ;

		this->setBody(iBody, mass, radius, posX, posY, posZ, velocityX, velocityY, velocityZ);

		if(!stream.good())
			return false;
	}

	// fill the bodies in the padding zone
	const unsigned long padding = (this->nVecs * mipp::vectorSize<T>()) - this->n;
	for(unsigned long iBody = this->n; iBody < this->n + padding; iBody++)
	{
		T posX, posY, posZ, velocityX, velocityY, velocityZ;

		posX = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * (5.0e8 * 1.33);
		posY = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 5.0e8;
		posZ = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 5.0e8 -10.0e8;

		velocityX = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
		velocityY = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
		velocityZ = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;

		this->setBody(iBody, 0, 0, posX, posY, posZ, velocityX, velocityY, velocityZ);
	}

	return true;
}

template <typename T>
void Bodies<T>::write(std::ostream& stream)
{
	stream << this->n << std::endl;

	for(unsigned long iBody = 0; iBody < this->n; iBody++)
		stream << this->masses      [iBody] << " "
		       << this->radiuses    [iBody] << " "
		       << this->positions.x [iBody] << " "
		       << this->positions.y [iBody] << " "
		       << this->positions.z [iBody] << " "
		       << this->velocities.x[iBody] << " "
		       << this->velocities.y[iBody] << " "
		       << this->velocities.z[iBody] << std::endl;
}

template <typename T>
void Bodies<T>::writeIntoFile(const std::string outputFileName)
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
std::ostream& operator<<(std::ostream &o, const Bodies<T>& s)
{
	s.write(o);
	return o;
}
