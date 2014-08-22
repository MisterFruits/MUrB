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

#include "./utils/myIntrinsics.h"
#include "Space.h"

#define dim 4

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
	this->masses = (T*)_mm_malloc(this->nBodies*sizeof(T),64);

	this->positions.x = (T*)_mm_malloc(this->nBodies*sizeof(T),64);
	this->positions.y = (T*)_mm_malloc(this->nBodies*sizeof(T),64);
	this->positions.z = (T*)_mm_malloc(this->nBodies*sizeof(T),64);

	this->speeds.x = (T*)_mm_malloc(this->nBodies*sizeof(T),64);
	this->speeds.y = (T*)_mm_malloc(this->nBodies*sizeof(T),64);
	this->speeds.z = (T*)_mm_malloc(this->nBodies*sizeof(T),64);

	this->accelerations.x = (T*)_mm_malloc(this->nBodies*sizeof(T),64);
	this->accelerations.y = (T*)_mm_malloc(this->nBodies*sizeof(T),64);
	this->accelerations.z = (T*)_mm_malloc(this->nBodies*sizeof(T),64);

	this->closestNeighborLen = (T*)_mm_malloc(this->nBodies*sizeof(T),64);
}

template <typename T>
void Space<T>::initBodiesRandomly()
{
	this->initBuffers();

	srand(123);
	for(unsigned long iBody = 0; iBody < this->nBodies; iBody++)
	{
		const T mass   = (rand() / (T) RAND_MAX) * 2000000;
		const T posX   = (rand() / (T) RAND_MAX) * 800;
		const T posY   = (rand() / (T) RAND_MAX) * 600;
		const T speedX = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 0.02;
		const T speedY = ((rand() - RAND_MAX/2) / (T) (RAND_MAX/2)) * 0.02;

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

	T mass, posX, posY, posZ, speedX, speedY, speedZ;
	for(unsigned long iBody = 0; iBody < this->nBodies; iBody++)
	{
		bodiesFile >> mass;
		bodiesFile >> posX;
		bodiesFile >> posY;
		bodiesFile >> posZ;
		bodiesFile >> speedX;
		bodiesFile >> speedY;
		bodiesFile >> speedZ;

		this->initBody(iBody, mass, posX, posY, posZ, speedX, speedY, speedZ);

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

	this->closestNeighborLen[iBody] = std::numeric_limits<T>::infinity();
}

template <typename T>
void Space<T>::computeBodiesAcceleration()
{
	// flops ~= nBody^2 * 12
	for(unsigned long iBody = 0; iBody < this->nBodies; iBody+=dim)
	{
		//pos1 = this->positions.x[iBody];
		// flops ~= nBody * 12
		for(unsigned long jBody = 0; jBody < this->nBodies; jBody+=dim)
			if(iBody != jBody)
//				this->computeAccelerationBetweenTwoBodies(iBody, jBody, dim); // 12 flops
				//this->vectorComputeAccelerationBetweenBodies(iBody, jBody, dim); // 12 flops
				this->intrinComputeAccelerationBetweenBodies(iBody, jBody, dim); // 12 flops
			else
				this->selfVectorComputeAccelerationBetweenBodies(iBody, dim);
	}
}

template <typename T>
void Space<T>::selfVectorComputeAccelerationBetweenBodies(const unsigned long iBody, const int vecDim)
{
	T vecX[dim],vecY[dim], vecLen[dim], acc[dim], accX[dim], accY[dim], sqrtVecLen[dim];
        for(int j=0; j<vecDim; j++)
        {
		for(int i=0; i<vecDim; i++)
		{
                vecX[i]       = this->positions.x[iBody+(iBody+j)%vecDim] - this->positions.x[iBody+i]; // 1 flop
                vecY[i]       = this->positions.y[iBody+(iBody+j)%vecDim] - this->positions.y[iBody+i]; // 1 flop
                vecLen[i]     = (vecX[i] * vecX[i]) + (vecY[i] * vecY[i]);                       // 3 flops
                if(vecLen[i]!=0){
			sqrtVecLen[i] = sqrt(vecLen[i]);    

			acc[i]  = this->masses[iBody+(iBody+j)%vecDim] / (vecLen[i] * sqrtVecLen[i]); // 2 flops
			accX[i] = acc[i] * vecX[i];                                  // 1 flop
			accY[i] = acc[i] * vecY[i];                                  // 1 flop

			this->accelerations.x[iBody+i] += accX[i]; // 1 flop
			this->accelerations.y[iBody+i] += accY[i]; // 1 flop
			
			if(sqrtVecLen[i] < this->closestNeighborLen[iBody+i])
				this->closestNeighborLen[iBody+i] = sqrtVecLen[i];
			}
		}
	}
}

template <typename T>
void Space<T>::vectorComputeAccelerationBetweenBodies(const unsigned long iBody, const unsigned long jBody, const int vecDim)
{
	T vecX[dim], vecY[dim], vecLen[dim], acc[dim], accX[dim], accY[dim], sqrtVecLen[dim];
	
	int jShuff[dim];
        for(int j=0; j<vecDim; j++)
        {
		for(int i=0; i<vecDim; i++)
		{
		jShuff[i]     = jBody+(j+i)%vecDim;

                vecX[i]       = this->positions.x[jShuff[i]] - this->positions.x[iBody+i]; // 1 flop
                vecY[i]       = this->positions.y[jShuff[i]] - this->positions.y[iBody+i]; // 1 flop
                vecLen[i]     = (vecX[i] * vecX[i]) + (vecY[i] * vecY[i]);                       // 3 flops
		sqrtVecLen[i] = sqrt(vecLen[i]);    

		acc[i]  = this->masses[jShuff[i]] / (vecLen[i] * sqrtVecLen[i]); // 2 flops
		accX[i] = acc[i] * vecX[i];                                  // 1 flop
		accY[i] = acc[i] * vecY[i];                                  // 1 flop

		this->accelerations.x[iBody+i] += accX[i]; // 1 flop
		this->accelerations.y[iBody+i] += accY[i]; // 1 flop
		
		this->closestNeighborLen[iBody+i] = std::min(sqrtVecLen[i],this->closestNeighborLen[iBody+i]);
		}
	}
}


template <typename T>
void Space<T>::intrinComputeAccelerationBetweenBodies(const unsigned long iBody, const unsigned long jBody, const int vecDim)
{
	vec px, py, rpx, rpy, masses, accx, accy, closest;
	vec vecX, vecY, vecLen, acc, sqrtVecLen;

	px       =  vec_load(&(this->positions.x[iBody]));
	py       =  vec_load(&(this->positions.y[iBody]));

	rpx      =  vec_load(&(this->positions.x[jBody]));
	rpy      =  vec_load(&(this->positions.y[jBody]));

	masses   =  vec_load(&(this->masses[jBody]));

	accx     =  vec_load(&(this->accelerations.x[iBody]));
	accy     =  vec_load(&(this->accelerations.y[iBody]));

	closest  =  vec_load(&(this->closestNeighborLen[iBody]));


	for(int i=0; i<vecDim; i++)
	{	

	vecX       = vec_sub(rpx,px); // 1 flop
	vecY       = vec_sub(rpy,py); // 1 flop
	vecLen     = vec_add( vec_mul(vecX,vecX) , vec_mul(vecY,vecY) ); // 3 flops

	sqrtVecLen = vec_sqrt(vecLen);    

	acc  = vec_div( masses , vec_mul(vecLen,sqrtVecLen) ); // 2 flops

	accx = vec_fmadd(acc, vecX, accx); // 2 flop
	accy = vec_fmadd(acc, vecY, accy); // 2 flop
	
	closest = vec_min(sqrtVecLen,closest);

	rpx = vec_permute(rpx,_MM_SHUFFLE(0,3,2,1));
	rpy = vec_permute(rpy,_MM_SHUFFLE(0,3,2,1));

	masses = vec_permute(masses,_MM_SHUFFLE(0,3,2,1));
	}

	vec_store(&(this->positions.x[iBody]), px);
	vec_store(&(this->positions.y[iBody]), py);

        vec_store(&(this->accelerations.x[iBody]), accx);
        vec_store(&(this->accelerations.y[iBody]), accy);

        vec_store(&(this->closestNeighborLen[iBody]), closest);

}

template <typename T>
void Space<T>::computeAccelerationBetweenTwoBodies(const unsigned long iBody, const unsigned long jBody, const int vecDim)
{
	T vecX[dim],vecY[dim], vecLen[dim], acc[dim], accX[dim], accY[dim], sqrtVecLen[dim];

	for(int i=0; i<vecDim; i++)
	{
		vecX[i]       = this->positions.x[jBody+i] - this->positions.x[iBody]; // 1 flop
		vecY[i]       = this->positions.y[jBody+i] - this->positions.y[iBody]; // 1 flop
		vecLen[i]     = (vecX[i] * vecX[i]) + (vecY[i] * vecY[i]);                       // 3 flops
		sqrtVecLen[i] = sqrt(vecLen[i]);                                        

		//if(vecLen[i] != 0){
		//	std::cout << "Collision at {" << this->positions.x[iBody] << ", "
		//				      << this->positions.y[iBody] << "}" << std::endl;
		//assert(vecLen[i] != 0);

		acc[i]  = this->masses[jBody+i] / (vecLen[i] * sqrtVecLen[i]); // 2 flops
		accX[i] = acc[i] * vecX[i];                                  // 1 flop
		accY[i] = acc[i] * vecY[i];                                  // 1 flop

		this->accelerations.x[iBody] += accX[i]; // 1 flop
		this->accelerations.y[iBody] += accY[i]; // 1 flop

		if(sqrtVecLen[i] < this->closestNeighborLen[iBody])
			this->closestNeighborLen[iBody] = sqrtVecLen[i];
		//this->closestNeighborLen[iBody] = std::min(sqrtVecLen[i],this->closestNeighborLen[iBody]);
		//}
	}
}

template <typename T>
void Space<T>::findTimeStep()
{
	this->dt = std::numeric_limits<T>::infinity();
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

		this->closestNeighborLen[iBody] = std::numeric_limits<T>::infinity();
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
		       << this->positions.z[iBody]     << " "
		       << this->speeds.x   [iBody]     << " "
		       << this->speeds.y   [iBody]     << " " 
		       << this->speeds.z   [iBody]     << std::endl;
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
