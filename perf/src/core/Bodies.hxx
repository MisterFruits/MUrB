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

	this->newVelocities.x = new T[this->n + this->padding];
	this->newVelocities.y = new T[this->n + this->padding];
	this->newVelocities.z = new T[this->n + this->padding];
#else
	this->masses = (T*)_mm_malloc((this->n + this->padding) * sizeof(T), mipp::RequiredAlignement);

	this->radiuses = (T*)_mm_malloc((this->n + this->padding) * sizeof(T), mipp::RequiredAlignement);

	this->positions.x = (T*)_mm_malloc((this->n + this->padding) * sizeof(T), mipp::RequiredAlignement);
	this->positions.y = (T*)_mm_malloc((this->n + this->padding) * sizeof(T), mipp::RequiredAlignement);
	this->positions.z = (T*)_mm_malloc((this->n + this->padding) * sizeof(T), mipp::RequiredAlignement);

	this->velocities.x = (T*)_mm_malloc((this->n + this->padding) * sizeof(T), mipp::RequiredAlignement);
	this->velocities.y = (T*)_mm_malloc((this->n + this->padding) * sizeof(T), mipp::RequiredAlignement);
	this->velocities.z = (T*)_mm_malloc((this->n + this->padding) * sizeof(T), mipp::RequiredAlignement);

	this->newVelocities.x = (T*)_mm_malloc((this->n + this->padding) * sizeof(T), mipp::RequiredAlignement);
	this->newVelocities.y = (T*)_mm_malloc((this->n + this->padding) * sizeof(T), mipp::RequiredAlignement);
	this->newVelocities.z = (T*)_mm_malloc((this->n + this->padding) * sizeof(T), mipp::RequiredAlignement);
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

	if(this->newVelocities.x != nullptr) {
		delete[] this->newVelocities.x;
		this->newVelocities.x = nullptr;
	}
	if(this->newVelocities.y != nullptr) {
		delete[] this->newVelocities.y;
		this->newVelocities.y = nullptr;
	}
	if(this->newVelocities.z != nullptr) {
		delete[] this->newVelocities.z;
		this->newVelocities.z = nullptr;
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

	if(this->newVelocities.x != nullptr) {
		_mm_free(this->newVelocities.x);
		this->newVelocities.x = nullptr;
	}
	if(this->newVelocities.y != nullptr) {
		_mm_free(this->newVelocities.y);
		this->newVelocities.y = nullptr;
	}
	if(this->newVelocities.z != nullptr) {
		_mm_free(this->newVelocities.z);
		this->newVelocities.z = nullptr;
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
const unsigned short& Bodies<T>::getPadding()
{
	return const_cast<const unsigned short&>(this->padding);
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
void Bodies<T>::initRandomly(const unsigned long randInit)
{
	this->allocateBuffers();

	srand(randInit);
	for(unsigned long iBody = 0; iBody < this->n; iBody++)
	{
		T mass, radius, posX, posY, posZ, velocityX, velocityY, velocityZ;

		mass = ((rand() / (T) RAND_MAX) * 5.0e21);

		//radius = mass * 0.5e-14;
		radius = mass * 0.5e-15;

		posX = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * (5.0e8 * 1.33);
		posY = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 5.0e8;
		posZ = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 5.0e8 -10.0e8;

		velocityX = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
		velocityY = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
		velocityZ = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;

		this->setBody(iBody, mass, radius, posX, posY, posZ, velocityX, velocityY, velocityZ);
	}

	// fill the bodies in the padding zone
	for(unsigned long iBody = this->n; iBody < this->n + this->padding; iBody++)
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

/* collisions 2D: work only if mi == mj (pool elastics collisions)
template <typename T>
void Bodies<T>::applyCollisions(std::vector<std::vector<unsigned long>> collisions)
{
	for(unsigned long iBody = 0; iBody < this->n; iBody++)
	{
		if(collisions[iBody].size())
		{
			std::cout << "Collisions for body n°" << iBody << ":" << std::endl;

			for(unsigned long jBody = 0; jBody < collisions[iBody].size(); jBody++)
			{
				std::cout << "  - body n°" << collisions[iBody][jBody] << std::endl;
			}
		}
	}


	for(unsigned long iBody = 0; iBody < this->n; iBody++)
		if(collisions[iBody].size())
		{
			assert(collisions[iBody].size() < 2);
			unsigned long jBody = collisions[iBody][0];

			// compute iBody ------------------------------------------------------------------------------------------
			T diffPosX = this->positions.x[jBody] - this->positions.x[iBody];
			T diffPosY = this->positions.y[jBody] - this->positions.y[iBody];

			T squareDist = (diffPosX * diffPosX) + (diffPosY * diffPosY);
			T dist = std::sqrt(squareDist);

			T normX = diffPosX / dist;
			T normY = diffPosY / dist;

			T perpendicularX = -normY;
			T perpendicularY =  normX;

			const T iVelocityX = this->velocities.x[iBody];
			const T iVelocityY = this->velocities.y[iBody];

			const T iSquareVelocityDist = (iVelocityX * iVelocityX) + (iVelocityY * iVelocityY);
			const T iVelocityDist = std::sqrt(iSquareVelocityDist);

			const T iCosTeta = ((perpendicularX * iVelocityX) + (perpendicularY * iVelocityY)) / iVelocityDist;
			const T iTeta = std::acos(iCosTeta);

			const T tmpIVelX = perpendicularX * (iCosTeta * iVelocityDist);
			const T tmpIVelY = perpendicularY * (iCosTeta * iVelocityDist);

			const T tmpJVelX = normX * (std::sin(iTeta) * iVelocityDist);
			const T tmpJVelY = normY * (std::sin(iTeta) * iVelocityDist);

			// compute jBody ------------------------------------------------------------------------------------------
			normX = -normX;
			normY = -normY;

			perpendicularX = -normY;
			perpendicularY =  normX;

			const T jVelocityX = this->velocities.x[jBody];
			const T jVelocityY = this->velocities.y[jBody];

			const T jSquareVelocityDist = (jVelocityX * jVelocityX) + (jVelocityY * jVelocityY);
			const T jVelocityDist = std::sqrt(jSquareVelocityDist);

			const T jCosTeta = ((perpendicularX * jVelocityX) + (perpendicularY * jVelocityY)) / jVelocityDist;
			const T jTeta = std::acos(jCosTeta);

			this->velocities.x[jBody] = perpendicularX * (jCosTeta * jVelocityDist) + tmpJVelX;
			this->velocities.y[jBody] = perpendicularY * (jCosTeta * jVelocityDist) + tmpJVelY;

			this->velocities.x[iBody] = normX * (std::sin(jTeta) * jVelocityDist) + tmpIVelX;
			this->velocities.y[iBody] = normY * (std::sin(jTeta) * jVelocityDist) + tmpIVelY;

			collisions[jBody].clear();
		}
}
*/

/* collisions 2D
template <typename T>
void Bodies<T>::applyCollisions(std::vector<std::vector<unsigned long>> collisions)
{
	for(unsigned long iBody = 0; iBody < this->n; iBody++)
	{
		if(collisions[iBody].size())
		{
			std::cout << "Collisions for body n°" << iBody << ":" << std::endl;

			for(unsigned long jBody = 0; jBody < collisions[iBody].size(); jBody++)
			{
				std::cout << "  - body n°" << collisions[iBody][jBody] << std::endl;
			}
		}
	}

	for(unsigned long iBody = 0; iBody < this->n; iBody++)
		if(collisions[iBody].size())
		{
			assert(collisions[iBody].size() < 2);
			unsigned long jBody = collisions[iBody][0];

			// 1. Find unit normal and unit tangent vectors
			T normX = this->positions.x[jBody] - this->positions.x[iBody];
			T normY = this->positions.y[jBody] - this->positions.y[iBody];

			T dij = std::sqrt((normX * normX) + (normY * normY));

			T uNormX = normX / dij;
			T uNormY = normY / dij;

			T uTangX = -uNormY;
			T uTangY =  uNormX;

			// 3 . Projecting the velocity vectors onto the unit normal and unit tangent vectors
			T viNorm = uNormX * this->velocities.x[iBody] + uNormY * this->velocities.y[iBody];
			T viTang = uTangX * this->velocities.x[iBody] + uTangY * this->velocities.y[iBody];

			T vjNorm = uNormX * this->velocities.x[jBody] + uNormY * this->velocities.y[jBody];
			T vjTang = uTangX * this->velocities.x[jBody] + uTangY * this->velocities.y[jBody];

			// 4. Find the new tangential velocities
			T viNewTang = viTang;
			T vjNewTang = vjTang;

			T iMass = this->masses[iBody];
			T jMass = this->masses[jBody];

			// 5. Find the new normal velocities
			T viNewNorm = (viNorm * (iMass - jMass) + (2.0 * jMass * vjNorm)) / (iMass + jMass);
			T vjNewNorm = (vjNorm * (jMass - iMass) + (2.0 * iMass * viNorm)) / (iMass + jMass);

			// 6. Convert the scalar normal and tangential velocities into vectors
			T viNewNormX = viNewNorm * uNormX;
			T viNewNormY = viNewNorm * uNormY;

			T viNewTangX = viNewTang * uTangX;
			T viNewTangY = viNewTang * uTangY;

			T vjNewNormX = vjNewNorm * uNormX;
			T vjNewNormY = vjNewNorm * uNormY;

			T vjNewTangX = vjNewTang * uTangX;
			T vjNewTangY = vjNewTang * uTangY;

			// 7. Find the final velocity vectors by adding the normal and tangential components
			this->velocities.x[iBody] = viNewNormX + viNewTangX;
			this->velocities.y[iBody] = viNewNormY + viNewTangY;

			this->velocities.x[jBody] = vjNewNormX + vjNewTangX;
			this->velocities.y[jBody] = vjNewNormY + vjNewTangY;

			collisions[jBody].clear();
		}
}
*/

/* collisions 3D */
template <typename T>
void Bodies<T>::applyCollisions(std::vector<std::vector<unsigned long>> collisions)
{
	/*
	for(unsigned long iBody = 0; iBody < this->n; iBody++)
		if(collisions[iBody].size())
		{
			std::cout << "Collisions for body n°" << iBody << ":" << std::endl;

			for(unsigned long jBody = 0; jBody < collisions[iBody].size(); jBody++)
				std::cout << "  - with body n°" << collisions[iBody][jBody] << std::endl;
		}
	*/

	for(unsigned long iBody = 0; iBody < this->n; iBody++)
		if(collisions[iBody].size())
		{
			assert(collisions[iBody].size() < 2);
			unsigned long jBody = collisions[iBody][0];

			// 1. Find unit normal and unit tangent 1 and tangent 2 vectors
			T normX = this->positions.x[jBody] - this->positions.x[iBody];
			T normY = this->positions.y[jBody] - this->positions.y[iBody];
			T normZ = this->positions.z[jBody] - this->positions.z[iBody];

			T dNorm = std::sqrt((normX * normX) + (normY * normY) + (normZ * normZ));

			T uNormX = normX / dNorm;
			T uNormY = normY / dNorm;
			T uNormZ = normZ / dNorm;
			//std::cout << "uNorm = {" << uNormX << ", " << uNormY << ", " << uNormZ << "}" << std::endl;

			// uTan1 and uTan2 vectors define the collision tangent plan (for 3D collisions)
			//T uTan1X = -uNormY;
			//T uTan1Y =  uNormZ;
			//T uTan1Z =  uNormX;

			T tan1X, tan1Y, tan1Z;
			if(uNormZ != 0)
			{
				tan1X = 1.0;
				tan1Y = 0.0;
				tan1Z = -(uNormX / uNormZ);
			}
			else if (uNormY != 0)
			{
				tan1X = 1.0;
				tan1Y = -(uNormX / uNormY);
				tan1Z = 0.0;
			}
			else if (uNormX != 0)
			{
				tan1X = -(uNormY / uNormX);
				tan1Y = 1.0;
				tan1Z = 0.0;
			}
			else
			{
				std::cout << "Error: uNorm = {" << uNormX << ", " << uNormY << ", " << uNormZ << "}" << std::endl;
				exit(-1);
			}

			T dTan1 = std::sqrt((tan1X * tan1X) + (tan1Y * tan1Y) + (tan1Z * tan1Z));

			T uTan1X = tan1X / dTan1;
			T uTan1Y = tan1Y / dTan1;
			T uTan1Z = tan1Z / dTan1;

			//std::cout << "uTan1 = {" << uTan1X << ", " << uTan1Y << ", " << uTan1Z << "}" << std::endl;

			//std::cout << "uNorm.uTan1 = " << ((uNormX * uTan1X) + (uNormY * uTan1Y) + (uNormZ * uTan1Z)) << std::endl;

			// uTan2 vector = (uNormX vector) ^ (uTan1 vector); (cross product or vector product)
			T uTan2X = (uNormY * uTan1Z) - (uNormZ * uTan1Y);
			T uTan2Y = (uNormZ * uTan1X) - (uNormX * uTan1Z);
			T uTan2Z = (uNormX * uTan1Y) - (uNormY * uTan1X);
			//std::cout << "uTan2 = {" << uTan2X << ", " << uTan2Y << ", " << uTan2Z << "}" << std::endl;

			//std::cout << "uNorm.uTan2 = " << ((uNormX * uTan2X) + (uNormY * uTan2Y) + (uNormZ * uTan2Z)) << std::endl;
			//std::cout << "uTan1.uTan2 = " << ((uTan1X * uTan2X) + (uTan1Y * uTan2Y) + (uTan1Z * uTan2Z)) << std::endl;

			// 2. Create the initial (before the collision) velocity vectors, iVel and jVel
			T iVelX = this->velocities.x[iBody];
			T iVelY = this->velocities.y[iBody];
			T iVelZ = this->velocities.z[iBody];

			T jVelX = this->velocities.x[jBody];
			T jVelY = this->velocities.y[jBody];
			T jVelZ = this->velocities.z[jBody];

			// 3. Projecting the velocity vectors onto the unit normal and unit tangent vectors
			// (scalar product or dot product)
			T viNorm = (uNormX * iVelX) + (uNormY * iVelY) + (uNormZ * iVelZ);
			T viTan1 = (uTan1X * iVelX) + (uTan1Y * iVelY) + (uTan1Z * iVelZ);
			T viTan2 = (uTan2X * iVelX) + (uTan2Y * iVelY) + (uTan2Z * iVelZ);

			T vjNorm = (uNormX * jVelX) + (uNormY * jVelY) + (uNormZ * jVelZ);
			T vjTan1 = (uTan1X * jVelX) + (uTan1Y * jVelY) + (uTan1Z * jVelZ);
			T vjTan2 = (uTan2X * jVelX) + (uTan2Y * jVelY) + (uTan2Z * jVelZ);

			// 4. Find the new tangential velocities
			T viNewTan1 = viTan1;
			T viNewTan2 = viTan2;

			T vjNewTan1 = vjTan1;
			T vjNewTan2 = vjTan2;

			T iMass = this->masses[iBody];
			T jMass = this->masses[jBody];

			// 5. Find the new normal velocities (apply elastic collision based on momentum and kinetic energy conservation)
			T viNewNorm = (viNorm * (iMass - jMass) + (2.0 * jMass * vjNorm)) / (iMass + jMass);
			T vjNewNorm = (vjNorm * (jMass - iMass) + (2.0 * iMass * viNorm)) / (iMass + jMass);

			// 6. Convert the scalar normal and tangential velocities into vectors
			T viNewNormX = viNewNorm * uNormX;
			T viNewNormY = viNewNorm * uNormY;
			T viNewNormZ = viNewNorm * uNormZ;

			T viNewTan1X = viNewTan1 * uTan1X;
			T viNewTan1Y = viNewTan1 * uTan1Y;
			T viNewTan1Z = viNewTan1 * uTan1Z;

			T viNewTan2X = viNewTan2 * uTan2X;
			T viNewTan2Y = viNewTan2 * uTan2Y;
			T viNewTan2Z = viNewTan2 * uTan2Z;

			T vjNewNormX = vjNewNorm * uNormX;
			T vjNewNormY = vjNewNorm * uNormY;
			T vjNewNormZ = vjNewNorm * uNormZ;

			T vjNewTan1X = vjNewTan1 * uTan1X;
			T vjNewTan1Y = vjNewTan1 * uTan1Y;
			T vjNewTan1Z = vjNewTan1 * uTan1Z;

			T vjNewTan2X = vjNewTan2 * uTan2X;
			T vjNewTan2Y = vjNewTan2 * uTan2Y;
			T vjNewTan2Z = vjNewTan2 * uTan2Z;

			// 7. Find the final velocity vectors by adding the normal and tangential components
			this->velocities.x[iBody] = viNewNormX + viNewTan1X + viNewTan2X;
			this->velocities.y[iBody] = viNewNormY + viNewTan1Y + viNewTan2Y;
			this->velocities.z[iBody] = viNewNormZ + viNewTan1Z + viNewTan2Z;

			this->velocities.x[jBody] = vjNewNormX + vjNewTan1X + vjNewTan2X;
			this->velocities.y[jBody] = vjNewNormY + vjNewTan1Y + vjNewTan2Y;
			this->velocities.z[jBody] = vjNewNormZ + vjNewTan1Z + vjNewTan2Z;

			collisions[jBody].clear();
		}
}

/* multi collisions 3D: not working */
template <typename T>
void Bodies<T>::applyMultiCollisions(std::vector<std::vector<unsigned long>> collisions)
{
	for(unsigned long iBody = 0; iBody < this->n; iBody++)
		if(collisions[iBody].size())
		{
			std::cout << "Collisions for body n°" << iBody << ":" << std::endl;

			for(unsigned long jBody = 0; jBody < collisions[iBody].size(); jBody++)
				std::cout << "  - with body n°" << collisions[iBody][jBody] << std::endl;
		}

	for(unsigned long iBody = 0; iBody < this->n; iBody++)
	{
		this->newVelocities.x[iBody] = 0.0;
		this->newVelocities.y[iBody] = 0.0;
		this->newVelocities.z[iBody] = 0.0;
	}

	for(unsigned long iBody = 0; iBody < this->n; iBody++)
	{
		for(unsigned long iCollision = 0; iCollision < collisions[iBody].size(); iCollision++)
		{
			unsigned long jBody = collisions[iBody][0];

			if(jBody > iBody)
			{
				// 1. Find unit normal and unit tangent 1 and tangent 2 vectors
				T normX = this->positions.x[jBody] - this->positions.x[iBody];
				T normY = this->positions.y[jBody] - this->positions.y[iBody];
				T normZ = this->positions.z[jBody] - this->positions.z[iBody];

				T dNorm = std::sqrt((normX * normX) + (normY * normY) + (normZ * normZ));

				T uNormX = normX / dNorm;
				T uNormY = normY / dNorm;
				T uNormZ = normZ / dNorm;

				// uTan1 and uTan2 vectors define the collision tangent plan (for 3D collisions)
				T tan1X, tan1Y, tan1Z;
				if(uNormZ != 0)
				{
					tan1X = 1.0;
					tan1Y = 0.0;
					tan1Z = -(uNormX / uNormZ);
				}
				else if (uNormY != 0)
				{
					tan1X = 1.0;
					tan1Y = -(uNormX / uNormY);
					tan1Z = 0.0;
				}
				else if (uNormX != 0)
				{
					tan1X = -(uNormY / uNormX);
					tan1Y = 1.0;
					tan1Z = 0.0;
				}
				else
				{
					std::cout << "Error: uNorm = {" << uNormX << ", " << uNormY << ", " << uNormZ << "}" << std::endl;
					exit(-1);
				}

				T dTan1 = std::sqrt((tan1X * tan1X) + (tan1Y * tan1Y) + (tan1Z * tan1Z));

				T uTan1X = tan1X / dTan1;
				T uTan1Y = tan1Y / dTan1;
				T uTan1Z = tan1Z / dTan1;

				// uTan2 vector = (uNormX vector) ^ (uTan1 vector); (cross product or vector product)
				T uTan2X = (uNormY * uTan1Z) - (uNormZ * uTan1Y);
				T uTan2Y = (uNormZ * uTan1X) - (uNormX * uTan1Z);
				T uTan2Z = (uNormX * uTan1Y) - (uNormY * uTan1X);

				// 2. Create the initial (before the collision) velocity vectors, iVel and jVel
				T iVelX = this->velocities.x[iBody];
				T iVelY = this->velocities.y[iBody];
				T iVelZ = this->velocities.z[iBody];

				T jVelX = this->velocities.x[jBody];
				T jVelY = this->velocities.y[jBody];
				T jVelZ = this->velocities.z[jBody];

				// 3. Projecting the velocity vectors onto the unit normal and unit tangent vectors
				// (scalar product or dot product)
				T viNorm = (uNormX * iVelX) + (uNormY * iVelY) + (uNormZ * iVelZ);
				T viTan1 = (uTan1X * iVelX) + (uTan1Y * iVelY) + (uTan1Z * iVelZ);
				T viTan2 = (uTan2X * iVelX) + (uTan2Y * iVelY) + (uTan2Z * iVelZ);

				T vjNorm = (uNormX * jVelX) + (uNormY * jVelY) + (uNormZ * jVelZ);
				T vjTan1 = (uTan1X * jVelX) + (uTan1Y * jVelY) + (uTan1Z * jVelZ);
				T vjTan2 = (uTan2X * jVelX) + (uTan2Y * jVelY) + (uTan2Z * jVelZ);

				// 4. Find the new tangential velocities
				T viNewTan1 = viTan1;
				T viNewTan2 = viTan2;

				T vjNewTan1 = vjTan1;
				T vjNewTan2 = vjTan2;

				T iMass = this->masses[iBody];
				T jMass = this->masses[jBody];

				// 5. Find the new normal velocities (apply elastic collision based on momentum and kinetic energy conservation)
				T viNewNorm = (viNorm * (iMass - jMass) + (2.0 * jMass * vjNorm)) /
				              ((iMass + jMass) * collisions[iBody].size() * collisions[jBody].size());
				T vjNewNorm = (vjNorm * (jMass - iMass) + (2.0 * iMass * viNorm)) /
				              ((iMass + jMass) * collisions[iBody].size() * collisions[jBody].size());

				// 6. Convert the scalar normal and tangential velocities into vectors
				T viNewNormX = viNewNorm * uNormX;
				T viNewNormY = viNewNorm * uNormY;
				T viNewNormZ = viNewNorm * uNormZ;

				T viNewTan1X = viNewTan1 * uTan1X;
				T viNewTan1Y = viNewTan1 * uTan1Y;
				T viNewTan1Z = viNewTan1 * uTan1Z;

				T viNewTan2X = viNewTan2 * uTan2X;
				T viNewTan2Y = viNewTan2 * uTan2Y;
				T viNewTan2Z = viNewTan2 * uTan2Z;

				T vjNewNormX = vjNewNorm * uNormX;
				T vjNewNormY = vjNewNorm * uNormY;
				T vjNewNormZ = vjNewNorm * uNormZ;

				T vjNewTan1X = vjNewTan1 * uTan1X;
				T vjNewTan1Y = vjNewTan1 * uTan1Y;
				T vjNewTan1Z = vjNewTan1 * uTan1Z;

				T vjNewTan2X = vjNewTan2 * uTan2X;
				T vjNewTan2Y = vjNewTan2 * uTan2Y;
				T vjNewTan2Z = vjNewTan2 * uTan2Z;

				// 7. Find the final velocity vectors by adding the normal and tangential components
				this->newVelocities.x[iBody] += viNewNormX + viNewTan1X + viNewTan2X;
				this->newVelocities.y[iBody] += viNewNormY + viNewTan1Y + viNewTan2Y;
				this->newVelocities.z[iBody] += viNewNormZ + viNewTan1Z + viNewTan2Z;

				this->newVelocities.x[jBody] += vjNewNormX + vjNewTan1X + vjNewTan2X;
				this->newVelocities.y[jBody] += vjNewNormY + vjNewTan1Y + vjNewTan2Y;
				this->newVelocities.z[jBody] += vjNewNormZ + vjNewTan1Z + vjNewTan2Z;
			}
		}
	}

	for(unsigned long iBody = 0; iBody < this->n; iBody++)
	{
		for(unsigned long iCollision = 0; iCollision < collisions[iBody].size(); iCollision++)
		{
			unsigned long jBody = collisions[iBody][0];

			this->velocities.x[iBody] = this->newVelocities.x[iBody];
			this->velocities.y[iBody] = this->newVelocities.y[iBody];
			this->velocities.z[iBody] = this->newVelocities.z[iBody];

			this->velocities.x[jBody] = this->newVelocities.x[jBody];
			this->velocities.y[jBody] = this->newVelocities.y[jBody];
			this->velocities.z[jBody] = this->newVelocities.z[jBody];
		}
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
	for(unsigned long iBody = this->n; iBody < this->n + this->padding; iBody++)
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
void Bodies<T>::write(std::ostream& stream, bool writeN)
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
bool Bodies<T>::writeIntoFile(const std::string outputFileName)
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
bool Bodies<T>::writeIntoFileMPI(const std::string outputFileName, const unsigned long MPINBodies)
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
