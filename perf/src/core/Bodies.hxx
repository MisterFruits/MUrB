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
#include <sys/stat.h>

#include "../utils/myIntrinsicsPlusPlus.h"

#include "Bodies.h"

template <typename T>
Bodies<T>::Bodies()
	: n             (0),
	  masses        (nullptr),
	  radiuses      (nullptr),
	  nVecs         (0),
	  padding       (0),
	  allocatedBytes(0)
{
}

template <typename T>
Bodies<T>::Bodies(const unsigned long n, const unsigned long randInit)
	: n             (n),
	  masses        (nullptr),
	  radiuses      (nullptr),
	  nVecs         (ceil((T) n / (T) mipp::vectorSize<T>())),
	  padding       ((this->nVecs * mipp::vectorSize<T>()) - this->n),
	  allocatedBytes(0)
{
	assert(n > 0);
	this->initRandomly(randInit);
}

template <typename T>
Bodies<T>::Bodies(const std::string inputFileName)
	: n             (0),
	  masses        (nullptr),
	  radiuses      (nullptr),
	  nVecs         (0),
	  padding       (0),
	  allocatedBytes(0)
{
	if(!this->initFromFile(inputFileName))
		exit(-1);
}

template <typename T>
Bodies<T>::Bodies(const Bodies<T>& bodies)
	: n             (bodies.n),
	  masses        (bodies.masses),
	  radiuses      (bodies.radiuses),
	  nVecs         (bodies.nVecs),
	  padding       (bodies.padding),
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
void Bodies<T>::hardCopy(const Bodies<T>& bodies)
{
	assert((this->n + this->padding) == (bodies.n + bodies.padding));

	this->n              = bodies.n;
	this->padding        = bodies.padding;
	this->nVecs          = bodies.nVecs;
	this->allocatedBytes = bodies.allocatedBytes;

	for(unsigned long iBody = 0; iBody < this->n + this->padding; iBody++)
	{
		this->masses      [iBody] = bodies.masses      [iBody];
		this->radiuses    [iBody] = bodies.radiuses    [iBody];
		this->positions.x [iBody] = bodies.positions.x [iBody];
		this->positions.y [iBody] = bodies.positions.y [iBody];
		this->positions.z [iBody] = bodies.positions.z [iBody];
		this->velocities.x[iBody] = bodies.velocities.x[iBody];
		this->velocities.y[iBody] = bodies.velocities.y[iBody];
		this->velocities.z[iBody] = bodies.velocities.z[iBody];
	}
}

template <typename T>
void Bodies<T>::allocateBuffers()
{
#ifdef __ARM_NEON__
	this->masses = new T[this->n + this->padding];

	this->radiuses = new T[this->n + this->padding];

	this->positions.x = new T[this->n + this->padding];
	this->positions.y = new T[this->n + this->padding];
	this->positions.z = new T[this->n + this->padding];

	this->velocities.x = new T[this->n + this->padding];
	this->velocities.y = new T[this->n + this->padding];
	this->velocities.z = new T[this->n + this->padding];
#else
	this->masses = (T*)_mm_malloc((this->n + this->padding) * sizeof(T), mipp::RequiredAlignement);

	this->radiuses = (T*)_mm_malloc((this->n + this->padding) * sizeof(T), mipp::RequiredAlignement);

	this->positions.x = (T*)_mm_malloc((this->n + this->padding) * sizeof(T), mipp::RequiredAlignement);
	this->positions.y = (T*)_mm_malloc((this->n + this->padding) * sizeof(T), mipp::RequiredAlignement);
	this->positions.z = (T*)_mm_malloc((this->n + this->padding) * sizeof(T), mipp::RequiredAlignement);

	this->velocities.x = (T*)_mm_malloc((this->n + this->padding) * sizeof(T), mipp::RequiredAlignement);
	this->velocities.y = (T*)_mm_malloc((this->n + this->padding) * sizeof(T), mipp::RequiredAlignement);
	this->velocities.z = (T*)_mm_malloc((this->n + this->padding) * sizeof(T), mipp::RequiredAlignement);
#endif

	this->allocatedBytes = (this->n + this->padding) * sizeof(T) * 11;
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
const unsigned long& Bodies<T>::getN() const
{
	return const_cast<const unsigned long&>(this->n);
}

template <typename T>
const unsigned long& Bodies<T>::getNVecs() const
{
	return const_cast<const unsigned long&>(this->nVecs);
}

template <typename T>
const unsigned short& Bodies<T>::getPadding() const
{
	return const_cast<const unsigned short&>(this->padding);
}

template <typename T>
const T* Bodies<T>::getMasses() const
{
	return const_cast<const T*>(this->masses);
}

template <typename T>
const T* Bodies<T>::getRadiuses() const
{
	return const_cast<const T*>(this->radiuses);
}

template <typename T>
const T* Bodies<T>::getPositionsX() const
{
	return const_cast<const T*>(this->positions.x);
}

template <typename T>
const T* Bodies<T>::getPositionsY() const
{
	return const_cast<const T*>(this->positions.y);
}

template <typename T>
const T* Bodies<T>::getPositionsZ() const
{
	return const_cast<const T*>(this->positions.z);
}

template <typename T>
const T* Bodies<T>::getVelocitiesX() const
{
	return const_cast<const T*>(this->velocities.x);
}

template <typename T>
const T* Bodies<T>::getVelocitiesY() const
{
	return const_cast<const T*>(this->velocities.y);
}

template <typename T>
const T* Bodies<T>::getVelocitiesZ() const
{
	return const_cast<const T*>(this->velocities.z);
}

template <typename T>
const float& Bodies<T>::getAllocatedBytes() const
{
	return this->allocatedBytes;
}

template <typename T>
void Bodies<T>::setBody(const unsigned long &iBody,
                        const T &mi, const T &radi,
                        const T &qiX, const T &qiY, const T &qiZ,
                        const T &viX, const T &viY, const T &viZ)
{
	this->masses[iBody] = mi;

	this->radiuses[iBody] = radi;

	this->positions.x[iBody] = qiX;
	this->positions.y[iBody] = qiY;
	this->positions.z[iBody] = qiZ;

	this->velocities.x[iBody] = viX;
	this->velocities.y[iBody] = viY;
	this->velocities.z[iBody] = viZ;
}

template <typename T>
void Bodies<T>::initRandomly(const unsigned long randInit)
{
	this->allocateBuffers();

	srand(randInit);
	for(unsigned long iBody = 0; iBody < this->n; iBody++)
	{
		T mi, radi, qiX, qiY, qiZ, viX, viY, viZ;

		mi = ((rand() / (T) RAND_MAX) * 5.0e21);

		radi = mi * 0.5e-14;
		//radius = mass * 0.5e-15;

		qiX = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * (5.0e8 * 1.33);
		qiY = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 5.0e8;
		qiZ = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 5.0e8 -10.0e8;

		viX = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
		viY = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
		viZ = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;

		this->setBody(iBody, mi, radi, qiX, qiY, qiZ, viX, viY, viZ);
	}

	// fill the bodies in the padding zone
	for(unsigned long iBody = this->n; iBody < this->n + this->padding; iBody++)
	{
		T qiX, qiY, qiZ, viX, viY, viZ;

		qiX = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * (5.0e8 * 1.33);
		qiY = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 5.0e8;
		qiZ = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 5.0e8 -10.0e8;

		viX = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
		viY = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
		viZ = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;

		this->setBody(iBody, 0, 0, qiX, qiY, qiZ, viX, viY, viZ);
	}
}

template <typename T>
bool Bodies<T>::initFromFile(const std::string inputFileName)
{
	std::ifstream bodiesFile;
	bodiesFile.open(inputFileName.c_str(), std::ios::in);

	if(!bodiesFile.is_open())
	{
		std::cout << "Can't open \"" << inputFileName << "\" file (reading)." << std::endl;
		return false;
	}

	bool isOk = this->read(bodiesFile);
	bodiesFile.close();

	if(!isOk)
	{
		std::cout << "Something bad occurred during the reading of \"" << inputFileName
		          << "\" file... exiting." << std::endl;
		return false;
	}

	return true;
}

template <typename T>
void Bodies<T>::updatePositionsAndVelocities(const vector3<T> &accelerations, T &dt)
{
	// flops = n * 18
	for(unsigned long iBody = 0; iBody < this->n; iBody++)
	{
		T mi, radi, qiX, qiY, qiZ, viX, viY, viZ;

		mi = this->masses[iBody];

		radi = this->radiuses[iBody];

		T aiXDt = accelerations.x[iBody] * dt;
		T aiYDt = accelerations.y[iBody] * dt;
		T aiZSt = accelerations.z[iBody] * dt;

		qiX = this->positions.x[iBody] + (this->velocities.x[iBody] + aiXDt * 0.5) * dt;
		qiY = this->positions.y[iBody] + (this->velocities.y[iBody] + aiYDt * 0.5) * dt;
		qiZ = this->positions.z[iBody] + (this->velocities.z[iBody] + aiZSt * 0.5) * dt;

		viX = this->velocities.x[iBody] + aiXDt;
		viY = this->velocities.y[iBody] + aiYDt;
		viZ = this->velocities.z[iBody] + aiZSt;

		this->setBody(iBody, mi, radi, qiX, qiY, qiZ, viX, viY, viZ);
	}
}

template <typename T>
bool Bodies<T>::readFromFile(const std::string inputFileName)
{
	assert(this->n == 0);
	return this->initFromFile(inputFileName);
}

template <typename T>
bool Bodies<T>::read(std::istream& stream)
{
	this->n = 0;
	stream >> this->n;

	this->nVecs   = ceil((T) this->n / (T) mipp::vectorSize<T>());
	this->padding = (this->nVecs * mipp::vectorSize<T>()) - this->n;

	if(this->n)
		this->allocateBuffers();
	else
		return false;

	for(unsigned long iBody = 0; iBody < this->n; iBody++)
	{
		T mi, radi, qiX, qiY, qiZ, viX, viY, viZ;

		stream >> mi;

		stream >> radi;

		stream >> qiX;
		stream >> qiY;
		stream >> qiZ;

		stream >> viX;
		stream >> viY;
		stream >> viZ;

		this->setBody(iBody, mi, radi, qiX, qiY, qiZ, viX, viY, viZ);

		if(!stream.good())
			return false;
	}

	// fill the bodies in the padding zone
	for(unsigned long iBody = this->n; iBody < this->n + this->padding; iBody++)
	{
		T qiX, qiY, qiZ, viX, viY, viZ;

		qiX = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * (5.0e8 * 1.33);
		qiY = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 5.0e8;
		qiZ = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 5.0e8 -10.0e8;

		viX = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
		viY = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
		viZ = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;

		this->setBody(iBody, 0, 0, qiX, qiY, qiZ, viX, viY, viZ);
	}

	return true;
}

template <typename T>
void Bodies<T>::write(std::ostream& stream, bool writeN) const
{
	if(writeN)
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
bool Bodies<T>::writeIntoFile(const std::string outputFileName) const
{
	std::fstream bodiesFile(outputFileName.c_str(), std::ios_base::out);
	if(!bodiesFile.is_open())
	{
		std::cout << "Can't open \"" << outputFileName << "\" file (writing). Exiting..." << std::endl;
		return false;
	}

	this->write(bodiesFile);

	bodiesFile.close();

	return true;
}

template <typename T>
bool Bodies<T>::writeIntoFileMPI(const std::string outputFileName, const unsigned long MPINBodies) const
{
	std::fstream bodiesFile;

	if(MPINBodies)
		bodiesFile.open(outputFileName.c_str(), std::fstream::out | std::fstream::trunc);
	else
		bodiesFile.open(outputFileName.c_str(), std::fstream::out | std::fstream::app);

	if(!bodiesFile.is_open())
	{
		std::cout << "Can't open \"" << outputFileName << "\" file (writing). Exiting..." << std::endl;
		return false;
	}

	if(MPINBodies)
		bodiesFile << MPINBodies << std::endl;

	bool writeN = false;
	this->write(bodiesFile, writeN);

	bodiesFile.close();

	return true;
}

template <typename T>
std::ostream& operator<<(std::ostream &o, const Bodies<T>& s)
{
	s.write(o);
	return o;
}
