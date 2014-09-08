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
	: nBodies(nBodies),
	  nVecs(ceil(nBodies * 1.0 / VECTOR_SIZE)),
	  masses(NULL),
	  radiuses(NULL),
	  closestNeighborDist(NULL),
	  dt(std::numeric_limits<T>::infinity()),
	  dtConstant(false)
{
	assert(nBodies > 0);
	this->initBodiesRandomly();
}

template <typename T>
Space<T>::Space(const std::string inputFileName)
	: nBodies(0),
	  nVecs(0),
	  masses(NULL),
	  radiuses(NULL),
	  closestNeighborDist(NULL),
	  dt(std::numeric_limits<T>::infinity()),
	  dtConstant(false)
{
	this->initBodiesFromFile(inputFileName);
}

template <typename T>
void Space<T>::allocateBuffers()
{
	this->masses = (vec_t<T>*)_mm_malloc(this->nVecs * sizeof(vec_t<T>), REQUIRED_ALIGNEMENT);

	this->radiuses = (vec_t<T>*)_mm_malloc(this->nVecs * sizeof(vec_t<T>), REQUIRED_ALIGNEMENT);

	this->positions.x = (vec_t<T>*)_mm_malloc(this->nVecs * sizeof(vec_t<T>), REQUIRED_ALIGNEMENT);
	this->positions.y = (vec_t<T>*)_mm_malloc(this->nVecs * sizeof(vec_t<T>), REQUIRED_ALIGNEMENT);
	this->positions.z = (vec_t<T>*)_mm_malloc(this->nVecs * sizeof(vec_t<T>), REQUIRED_ALIGNEMENT);

	this->speeds.x = (vec_t<T>*)_mm_malloc(this->nVecs * sizeof(vec_t<T>), REQUIRED_ALIGNEMENT);
	this->speeds.y = (vec_t<T>*)_mm_malloc(this->nVecs * sizeof(vec_t<T>), REQUIRED_ALIGNEMENT);
	this->speeds.z = (vec_t<T>*)_mm_malloc(this->nVecs * sizeof(vec_t<T>), REQUIRED_ALIGNEMENT);

	this->accelerations.x = (vec_t<T>*)_mm_malloc(this->nVecs * sizeof(vec_t<T>), REQUIRED_ALIGNEMENT);
	this->accelerations.y = (vec_t<T>*)_mm_malloc(this->nVecs * sizeof(vec_t<T>), REQUIRED_ALIGNEMENT);
	this->accelerations.z = (vec_t<T>*)_mm_malloc(this->nVecs * sizeof(vec_t<T>), REQUIRED_ALIGNEMENT);

	this->closestNeighborDist = (vec_t<T>*)_mm_malloc(this->nVecs * sizeof(vec_t<T>), REQUIRED_ALIGNEMENT);
}

template <typename T>
Space<T>::~Space() {
	if(this->masses)
		_mm_free(this->masses);

	if(this->radiuses)
		_mm_free(this->radiuses);

	if(this->positions.x)
		_mm_free(this->positions.x);
	if(this->positions.y)
		_mm_free(this->positions.y);
	if(this->positions.z)
		_mm_free(this->positions.z);

	if(this->speeds.x)
		_mm_free(this->speeds.x);
	if(this->speeds.y)
		_mm_free(this->speeds.y);
	if(this->speeds.z)
		_mm_free(this->speeds.z);

	if(this->accelerations.x)
		_mm_free(this->accelerations.x);
	if(this->accelerations.y)
		_mm_free(this->accelerations.y);
	if(this->accelerations.z)
		_mm_free(this->accelerations.z);

	if(this->closestNeighborDist)
		_mm_free(this->closestNeighborDist);
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
void Space<T>::initBodiesRandomly()
{
	this->allocateBuffers();

	srand(123);
	for(unsigned long iVec = 0; iVec < this->nVecs; iVec++)
	{
		for(unsigned short iBody = 0; iBody < VECTOR_SIZE; iBody++)
		{
			unsigned long realBody = iBody + iVec * VECTOR_SIZE;
			if(realBody < this->nBodies)
				this->masses[iVec].vec_data[iBody] = ((rand() / (T) RAND_MAX) * 5.0e21);
			else // fake body just to fit into the last vector
				this->masses[iVec].vec_data[iBody] = 0;

			this->radiuses[iVec].vec_data[iBody] = this->masses[iVec].vec_data[iBody] * 0.6e-15;

			this->positions.x[iVec].vec_data[iBody] = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * (5.0e8 * 1.33);
			this->positions.y[iVec].vec_data[iBody] = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 5.0e8;
			this->positions.z[iVec].vec_data[iBody] = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 5.0e8 -10.0e8;

			this->speeds.x[iVec].vec_data[iBody] = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
			this->speeds.y[iVec].vec_data[iBody] = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
			this->speeds.z[iVec].vec_data[iBody] = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;

			this->accelerations.x[iVec].vec_data[iBody] = 0;
			this->accelerations.y[iVec].vec_data[iBody] = 0;
			this->accelerations.z[iVec].vec_data[iBody] = 0;

			this->closestNeighborDist[iVec].vec_data[iBody] = std::numeric_limits<T>::infinity();

		}
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
	for(unsigned long iVec = 0; iVec < this->nVecs; iVec++)
		for(unsigned long jVec = 0; jVec < this->nVecs; jVec++)
			if(iVec != jVec)
				this->iComputeAccelerationBetweenTwoVectorOfBodies(this->positions.x        [iVec].vec_data,
				                                                   this->positions.y        [iVec].vec_data,
				                                                   this->positions.z        [iVec].vec_data,
				                                                   this->accelerations.x    [iVec].vec_data,
				                                                   this->accelerations.y    [iVec].vec_data,
				                                                   this->accelerations.z    [iVec].vec_data,
				                                                   this->closestNeighborDist[iVec].vec_data,
				                                                   this->masses             [jVec].vec_data,
				                                                   this->positions.x        [jVec].vec_data,
				                                                   this->positions.y        [jVec].vec_data,
				                                                   this->positions.z        [jVec].vec_data);
			else
				for(unsigned short iBody = 0; iBody < VECTOR_SIZE; iBody++)
					for(unsigned short jBody = 0; jBody < VECTOR_SIZE; jBody++)
						if(iBody != jBody)
							this->computeAccelerationBetweenTwoBodies(this->positions.x        [iVec].vec_data[iBody],
							                                          this->positions.y        [iVec].vec_data[iBody],
							                                          this->positions.z        [iVec].vec_data[iBody],
							                                          this->accelerations.x    [iVec].vec_data[iBody],
							                                          this->accelerations.y    [iVec].vec_data[iBody],
							                                          this->accelerations.z    [iVec].vec_data[iBody],
							                                          this->closestNeighborDist[iVec].vec_data[iBody],
							                                          this->masses             [jVec].vec_data[jBody],
							                                          this->positions.x        [jVec].vec_data[jBody],
							                                          this->positions.y        [jVec].vec_data[jBody],
							                                          this->positions.z        [jVec].vec_data[jBody]);
}

template <typename T>
void Space<T>::iComputeBodiesAcceleration()
{
	const vec rG = vec_set1(G);

#pragma omp parallel for schedule(runtime) firstprivate(rG)
	for(unsigned long iVec = 0; iVec < this->nVecs; iVec++)
	{
		// load vectors
		const vec rIPosX = vec_load(this->positions.x[iVec].vec_data);
		const vec rIPosY = vec_load(this->positions.y[iVec].vec_data);
		const vec rIPosZ = vec_load(this->positions.z[iVec].vec_data);

		vec rIAccX = vec_load(this->accelerations.x[iVec].vec_data);
		vec rIAccY = vec_load(this->accelerations.y[iVec].vec_data);
		vec rIAccZ = vec_load(this->accelerations.z[iVec].vec_data);

		vec rIClosNeiDist = vec_load(this->closestNeighborDist[iVec].vec_data);

		for(unsigned long jVec = 0; jVec < this->nVecs; jVec++)
		{
			// load vectors
			vec rJMass = vec_load(this->masses[jVec].vec_data);
			vec rJPosX = vec_load(this->positions.x[jVec].vec_data);
			vec rJPosY = vec_load(this->positions.y[jVec].vec_data);
			vec rJPosZ = vec_load(this->positions.z[jVec].vec_data);

			if(iVec != jVec)
			{
				for(unsigned short iRot = 0; iRot < VECTOR_SIZE; iRot++)
				{
					this->iComputeAccelerationBetweenTwoBodies(rG,
					                                           rIPosX, rIPosY, rIPosZ,
					                                           rIAccX, rIAccY, rIAccZ,
					                                           rIClosNeiDist,
					                                           rJMass,
					                                           rJPosX, rJPosY, rJPosZ);

					// we make one useless rotate in the last iteration...
					rJMass = vec_rot(rJMass);
					rJPosX = vec_rot(rJPosX); rJPosY = vec_rot(rJPosY); rJPosZ = vec_rot(rJPosZ);
				}
			}
			else
			{
				for(unsigned short iRot = 1; iRot < VECTOR_SIZE; iRot++)
				{
					rJMass = vec_rot(rJMass);
					rJPosX = vec_rot(rJPosX); rJPosY = vec_rot(rJPosY); rJPosZ = vec_rot(rJPosZ);

					this->iComputeAccelerationBetweenTwoBodies(rG,
					                                           rIPosX, rIPosY, rIPosZ,
					                                           rIAccX, rIAccY, rIAccZ,
					                                           rIClosNeiDist,
					                                           rJMass,
					                                           rJPosX, rJPosY, rJPosZ);
				}
			}
		}

		// store vectors
		vec_store(this->accelerations.x[iVec].vec_data, rIAccX);
		vec_store(this->accelerations.y[iVec].vec_data, rIAccY);
		vec_store(this->accelerations.z[iVec].vec_data, rIAccZ);
		vec_store(this->closestNeighborDist[iVec].vec_data, rIClosNeiDist);
	}
}

template <typename T>
void Space<T>::computeAccelerationBetweenTwoVectorOfBodies(const T* __restrict iVecPosX,
                                                           const T* __restrict iVecPosY,
                                                           const T* __restrict iVecPosZ,
                                                                 T* __restrict iVecAccsX,
                                                                 T* __restrict iVecAccsY,
                                                                 T* __restrict iVecAccsZ,
                                                                 T* __restrict iClosNeiDist,
                                                           const T* __restrict jVecMasses,
                                                           const T* __restrict jVecPosX,
                                                           const T* __restrict jVecPosY,
                                                           const T* __restrict jVecPosZ)
{
	for(unsigned short iBody = 0; iBody < VECTOR_SIZE; iBody++)
		for(unsigned short jBody = 0; jBody < VECTOR_SIZE; jBody++)
			this->computeAccelerationBetweenTwoBodies(iVecPosX    [iBody],
			                                          iVecPosY    [iBody],
			                                          iVecPosZ    [iBody],
			                                          iVecAccsX   [iBody],
			                                          iVecAccsY   [iBody],
			                                          iVecAccsZ   [iBody],
			                                          iClosNeiDist[iBody],
			                                          jVecMasses  [jBody],
			                                          jVecPosX    [jBody],
			                                          jVecPosY    [jBody],
			                                          jVecPosZ    [jBody]);
}

template <typename T>
void Space<T>::computeAccelerationBetweenTwoBodies(const T &iPosX, const T &iPosY, const T &iPosZ,
                                                         T &iAccsX,      T &iAccsY,      T &iAccsZ,
                                                         T &iClosNeiDist,
                                                   const T &jMasses,
                                                   const T &jPosX, const T &jPosY, const T &jPosZ)
{
	const T diffPosX = jPosX - iPosX; // 1 flop
	const T diffPosY = jPosY - iPosY; // 1 flop
	const T diffPosZ = jPosZ - iPosZ; // 1 flop
	const T squareDist = (diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ); // 5 flops
	const T dist = std::sqrt(squareDist);
	assert(dist != 0);

	const T acc = G * jMasses / (squareDist * dist); // 3 flops
	iAccsX += acc * diffPosX; // 2 flop
	iAccsY += acc * diffPosY; // 2 flop
	iAccsZ += acc * diffPosZ; // 2 flop

	if(!this->dtConstant)
		if(dist < iClosNeiDist)
			iClosNeiDist = dist;
}

template <typename T>
void Space<T>::iComputeAccelerationBetweenTwoVectorOfBodies(const T* __restrict iVecPosX,
                                                            const T* __restrict iVecPosY,
                                                            const T* __restrict iVecPosZ,
                                                                  T* __restrict iVecAccsX,
                                                                  T* __restrict iVecAccsY,
                                                                  T* __restrict iVecAccsZ,
                                                                  T* __restrict iClosNeiDist,
                                                            const T* __restrict jVecMasses,
                                                            const T* __restrict jVecPosX,
                                                            const T* __restrict jVecPosY,
                                                            const T* __restrict jVecPosZ)
{
	assert(alignof(iVecPosX)     == REQUIRED_ALIGNEMENT);
	assert(alignof(iVecPosY)     == REQUIRED_ALIGNEMENT);
	assert(alignof(iVecPosZ)     == REQUIRED_ALIGNEMENT);
	assert(alignof(iVecAccsX)    == REQUIRED_ALIGNEMENT);
	assert(alignof(iVecAccsY)    == REQUIRED_ALIGNEMENT);
	assert(alignof(iVecAccsZ)    == REQUIRED_ALIGNEMENT);
	assert(alignof(iClosNeiDist) == REQUIRED_ALIGNEMENT);
	assert(alignof(jVecMasses)   == REQUIRED_ALIGNEMENT);
	assert(alignof(jVecPosX)     == REQUIRED_ALIGNEMENT);
	assert(alignof(jVecPosY)     == REQUIRED_ALIGNEMENT);
	assert(alignof(jVecPosZ)     == REQUIRED_ALIGNEMENT);

	// load vectors
	vec rG     = vec_set1(G);

	vec rIPosX = vec_load(iVecPosX);
	vec rIPosY = vec_load(iVecPosY);
	vec rIPosZ = vec_load(iVecPosZ);

	vec rJPosX = vec_load(jVecPosX);
	vec rJPosY = vec_load(jVecPosY);
	vec rJPosZ = vec_load(jVecPosZ);


	vec rJMass = vec_load(jVecMasses);

	vec rIAccX = vec_load(iVecAccsX);
	vec rIAccY = vec_load(iVecAccsY);
	vec rIAccZ = vec_load(iVecAccsZ);

	vec rIClosNeiDist = vec_load(iClosNeiDist);

	for(unsigned short iBody = 0; iBody < VECTOR_SIZE; iBody++)
	{
		this->iComputeAccelerationBetweenTwoBodies(rG,
		                                           rIPosX, rIPosY, rIPosZ,
		                                           rIAccX, rIAccY, rIAccZ,
		                                           rIClosNeiDist,
		                                           rJMass,
		                                           rJPosX, rJPosY, rJPosZ);

		// we make one useless rotate in the last iteration...
		rJPosX = vec_rot(rJPosX);
		rJPosY = vec_rot(rJPosY);
		rJPosZ = vec_rot(rJPosZ);
		rJMass = vec_rot(rJMass);
	}

	// stores vectors
	vec_store(iVecAccsX,    rIAccX);
	vec_store(iVecAccsY,    rIAccY);
	vec_store(iVecAccsZ,    rIAccZ);
	vec_store(iClosNeiDist, rIClosNeiDist);
}

template <typename T>
void Space<T>::iComputeAccelerationBetweenTwoBodies(const vec &rG,
                                                    const vec &rIPosX,
                                                    const vec &rIPosY,
                                                    const vec &rIPosZ,
                                                          vec &rIAccX,
                                                          vec &rIAccY,
                                                          vec &rIAccZ,
                                                          vec &rIClosNeiDist,
                                                          vec &rJMass,
                                                          vec &rJPosX,
                                                          vec &rJPosY,
                                                          vec &rJPosZ)
{
	//const T diffPosX = jPosX - iPosX;
	vec rDiffPosX = vec_sub(rJPosX, rIPosX);
	//const T diffPosY = jPosY - iPosY;
	vec rDiffPosY = vec_sub(rJPosY, rIPosY);
	//const T diffPosZ = jPosZ - iPosZ;
	vec rDiffPosZ = vec_sub(rJPosZ, rIPosZ);

	//const T squareDist = (diffPosX * diffPosX) + (diffPosY * diffPosY) + (diffPosZ * diffPosZ);
	vec rSquareDist = vec_set1(0);
	rSquareDist = vec_fmadd(rDiffPosX, rDiffPosX, rSquareDist);
	rSquareDist = vec_fmadd(rDiffPosY, rDiffPosY, rSquareDist);
	rSquareDist = vec_fmadd(rDiffPosZ, rDiffPosZ, rSquareDist);

	//const T dist = std::sqrt(squareDist);
	vec rDist = vec_sqrt(rSquareDist);

	//const T acc = G * jMasses / (squareDist * dist);
	vec rAcc = vec_div(vec_mul(rG, rJMass), vec_mul(rDist, rSquareDist));

	//iAccsX += acc * diffPosX;
	rIAccX = vec_fmadd(rAcc, rDiffPosX, rIAccX);
	//iAccsY += acc * diffPosY;
	rIAccY = vec_fmadd(rAcc, rDiffPosY, rIAccY);
	//iAccsZ += acc * diffPosZ;
	rIAccZ = vec_fmadd(rAcc, rDiffPosZ, rIAccZ);

	//if(!this->dtConstant)
	//	if(dist < iClosNeiDist)
	//		iClosNeiDist = dist;
	rIClosNeiDist = vec_min(rDist, rIClosNeiDist);
}

template <typename T>
void Space<T>::findTimeStep()
{
	// TODO: be careful, with fake bodies added at the end of the last vector, the dynamic time step is broken.
	//       It is necessary to launch the simulation with a number of bodies multiple of VECTOR_SIZE!
	if(!this->dtConstant)
	{
		this->dt = std::numeric_limits<T>::infinity();

		// flops = nBodies * 16
		for(unsigned long iVec = 0; iVec < this->nVecs; iVec++)
		{
			const T newDt = computeTimeStep(iVec); // 16 flops

			if(newDt < this->dt)
				this->dt = newDt;
		}
	}
}

template <typename T>
T Space<T>::computeTimeStep(const unsigned long iVec)
{
	T newDt;
	vec_t<T> vecNewDt;

	for(unsigned short iBody = 0; iBody < VECTOR_SIZE; iBody++)
	{
		// || lb.speed ||
		const T s = std::sqrt((this->speeds.x[iVec].vec_data[iBody] * this->speeds.x[iVec].vec_data[iBody]) +
		                      (this->speeds.y[iVec].vec_data[iBody] * this->speeds.y[iVec].vec_data[iBody]) +
		                      (this->speeds.z[iVec].vec_data[iBody] * this->speeds.z[iVec].vec_data[iBody])); // 5 flops

		// || lb.acceleration ||
		const T a = std::sqrt((this->accelerations.x[iVec].vec_data[iBody] * this->accelerations.x[iVec].vec_data[iBody]) +
		                      (this->accelerations.y[iVec].vec_data[iBody] * this->accelerations.y[iVec].vec_data[iBody]) +
		                      (this->accelerations.z[iVec].vec_data[iBody] * this->accelerations.z[iVec].vec_data[iBody])); // 5 flops

		/*
		 * compute dt
		 * solve:  (a/2)*dt^2 + s*dt + (-0.1)*ClosestNeighborDist = 0
		 * <=>     dt = [ (-s) +/-  sqrt( s^2 - 4 * (a/2) * (-0.1)*ClosestNeighborDist ) ] / [ 2 (a/2) ]
		 *
		 * dt should be positive (+/- becomes + because result of sqrt is positive)
		 * <=>     dt = [ -s + sqrt( s^2 + 0.2*ClosestNeighborDist*a) ] / a
		 */
		vecNewDt.vec_data[iBody] = (std::sqrt(s * s + 0.2 * a * this->closestNeighborDist[iVec].vec_data[iBody]) - s) / a; // 6 flops

		if(vecNewDt.vec_data[iBody] == 0)
			vecNewDt.vec_data[iBody] = std::numeric_limits<T>::epsilon() / a;
	}

	// looking for min dt in the vector of dt
	newDt = vecNewDt.vec_data[0];
	for(unsigned short iBody = 1; iBody < VECTOR_SIZE; iBody++)
	{
		if(vecNewDt.vec_data[iBody] < newDt)
			newDt = vecNewDt.vec_data[iBody];
	}

	return newDt;
}

template <typename T>
void Space<T>::updateBodiesPositionAndSpeed()
{
	// flops = nBodies * 18
	for(unsigned long iVec = 0; iVec < this->nVecs; iVec++)
	{
		for(unsigned short iBody = 0; iBody < VECTOR_SIZE; iBody++)
		{
			T accXMultDt = this->accelerations.x[iVec].vec_data[iBody] * this->dt; // 1 flop
			T accYMultDt = this->accelerations.y[iVec].vec_data[iBody] * this->dt; // 1 flop
			T accZMultDt = this->accelerations.z[iVec].vec_data[iBody] * this->dt; // 1 flop

			this->positions.x[iVec].vec_data[iBody] += (this->speeds.x[iVec].vec_data[iBody] + accXMultDt * 0.5) * this->dt; // 4 flops
			this->positions.y[iVec].vec_data[iBody] += (this->speeds.y[iVec].vec_data[iBody] + accYMultDt * 0.5) * this->dt; // 4 flops
			this->positions.z[iVec].vec_data[iBody] += (this->speeds.z[iVec].vec_data[iBody] + accZMultDt * 0.5) * this->dt; // 4 flops;

			this->speeds.x[iVec].vec_data[iBody] += accXMultDt; // 1 flop
			this->speeds.y[iVec].vec_data[iBody] += accYMultDt; // 1 flop
			this->speeds.z[iVec].vec_data[iBody] += accZMultDt; // 1 flop

			this->accelerations.x[iVec].vec_data[iBody] = 0;
			this->accelerations.y[iVec].vec_data[iBody] = 0;
			this->accelerations.z[iVec].vec_data[iBody] = 0;

			this->closestNeighborDist[iVec].vec_data[iBody] = std::numeric_limits<T>::infinity();
		}
	}
}

template <typename T>
bool Space<T>::read(std::istream& stream)
{
	this->nBodies = 0;
	stream >> this->nBodies;

	if(this->nBodies)
	{
		this->nVecs = ceil(this->nBodies * 1.0 / VECTOR_SIZE);
		this->allocateBuffers();
	}
	else
		return false;

	for(unsigned long iVec = 0; iVec < this->nVecs; iVec++)
	{
		for(unsigned short iBody = 0; iBody < VECTOR_SIZE; iBody++)
		{
			unsigned long realBody = iBody + iVec * VECTOR_SIZE;
			if(realBody < this->nBodies) // read from file
			{
				stream >> this->masses[iVec].vec_data[iBody];

				stream >> this->radiuses[iVec].vec_data[iBody];

				stream >> this->positions.x[iVec].vec_data[iBody];
				stream >> this->positions.y[iVec].vec_data[iBody];
				stream >> this->positions.z[iVec].vec_data[iBody];

				stream >> this->speeds.x[iVec].vec_data[iBody];
				stream >> this->speeds.y[iVec].vec_data[iBody];
				stream >> this->speeds.z[iVec].vec_data[iBody];
			}
			else // fake body just to fit into the last vector
			{
				this->masses[iVec].vec_data[iBody] = 0;

				this->radiuses[iVec].vec_data[iBody] = this->masses[iVec].vec_data[iBody] * 0.6e-15;

				this->positions.x[iVec].vec_data[iBody] = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * (5.0e8 * 1.33);
				this->positions.y[iVec].vec_data[iBody] = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 5.0e8;
				this->positions.z[iVec].vec_data[iBody] = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 5.0e8 -10.0e8;

				this->speeds.x[iVec].vec_data[iBody] = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
				this->speeds.y[iVec].vec_data[iBody] = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
				this->speeds.z[iVec].vec_data[iBody] = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 1.0e2;
			}

			this->accelerations.x[iVec].vec_data[iBody] = 0;
			this->accelerations.y[iVec].vec_data[iBody] = 0;
			this->accelerations.z[iVec].vec_data[iBody] = 0;

			this->closestNeighborDist[iVec].vec_data[iBody] = std::numeric_limits<T>::infinity();

			if(!stream.good())
				return false;
		}
	}

	return true;
}

template <typename T>
void Space<T>::write(std::ostream& stream)
{
	stream << this->nBodies << std::endl;

	for(unsigned long iVec = 0; iVec < this->nVecs; iVec++)
		for(unsigned short iBody = 0; iBody < VECTOR_SIZE; iBody++)
			if((iBody + iVec * VECTOR_SIZE) < this->nBodies) // do not write fake bodies
				stream << this->masses     [iVec].vec_data[iBody] << " "
				       << this->radiuses   [iVec].vec_data[iBody] << " "
				       << this->positions.x[iVec].vec_data[iBody] << " "
				       << this->positions.y[iVec].vec_data[iBody] << " "
				       << this->positions.z[iVec].vec_data[iBody] << " "
				       << this->speeds.x   [iVec].vec_data[iBody] << " "
				       << this->speeds.y   [iVec].vec_data[iBody] << " "
				       << this->speeds.z   [iVec].vec_data[iBody] << std::endl;
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
