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

#include "Bodies.h"

template <typename T>
Bodies<T>::Bodies(const unsigned long n)
	: n       (n),
	  masses  (NULL),
	  radiuses(NULL)
{
	assert(n > 0);
	this->initRandomly();
}

template <typename T>
Bodies<T>::Bodies(const std::string inputFileName)
	: n       (0),
	  masses  (NULL),
	  radiuses(NULL)
{
	this->initFromFile(inputFileName);
}

template <typename T>
void Bodies<T>::allocateBuffers()
{
	this->masses = new T[this->n];

	this->radiuses = new T[this->n];

	this->positions.x = new T[this->n];
	this->positions.y = new T[this->n];
	this->positions.z = new T[this->n];

	this->velocities.x = new T[this->n];
	this->velocities.y = new T[this->n];
	this->velocities.z = new T[this->n];
}

template <typename T>
Bodies<T>::~Bodies() {
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

	if(this->velocities.x)
		delete[] this->velocities.x;
	if(this->velocities.y)
		delete[] this->velocities.y;
	if(this->velocities.z)
		delete[] this->velocities.z;
}

template <typename T>
unsigned long Bodies<T>::getNBodies()
{
	return this->n;
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
void Bodies<T>::updatePositionsAndVelocities(vector3<T> &accelerations, T &dt)
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
