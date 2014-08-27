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

#define dim 4


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
void Space<T>::vectorComputeBodiesAcceleration()
{
	//nBodies²*16 flops
        for(unsigned long iBody = 0; iBody < this->nBodies; iBody+=dim)
        {
                for(unsigned long jBody = 0; jBody < this->nBodies; jBody+=dim)
			// 16 flops
                        if(iBody != jBody)
                                this->vectorComputeAccelerationBetweenBodies(iBody, jBody, dim);
                        else
                                this->selfVectorComputeAccelerationBetweenBodies(iBody, dim);
        }   
}

template <typename T>
void Space<T>::selfVectorComputeAccelerationBetweenBodies(const unsigned long iBody, const int vecDim)
{
        T vecX[dim], vecY[dim], vecZ[dim];
	T vecLen[dim], acc[dim], accX[dim], accY[dim], accZ[dim], sqrtVecLen[dim];
        for(int j=0; j<vecDim; j++)
        {   
                for(int i=0; i<vecDim; i++)
                {   
                vecX[i]       = this->positions.x[iBody+(iBody+j)%vecDim] - this->positions.x[iBody+i]; 
                vecY[i]       = this->positions.y[iBody+(iBody+j)%vecDim] - this->positions.y[iBody+i];
                vecZ[i]       = this->positions.z[iBody+(iBody+j)%vecDim] - this->positions.z[iBody+i];
                vecLen[i]     = vecX[i] * vecX[i] + vecY[i] * vecY[i] + vecZ[i]*vecZ[i];                       
                if(vecLen[i]!=0){
                        sqrtVecLen[i] = sqrt(vecLen[i]);    

                        acc[i]  = this->masses[iBody+(iBody+j)%vecDim] / (vecLen[i] * sqrtVecLen[i]); 
                        accX[i] = acc[i] * vecX[i];                                 
                        accY[i] = acc[i] * vecY[i];                                  
                        accZ[i] = acc[i] * vecZ[i];                                 

                        this->accelerations.x[iBody+i] += accX[i]; 
                        this->accelerations.y[iBody+i] += accY[i];
                        this->accelerations.z[iBody+i] += accZ[i];
    
                        if(sqrtVecLen[i] < this->closestNeighborLen[iBody+i])
                                this->closestNeighborLen[iBody+i] = sqrtVecLen[i];
                        }   
                }   
        }   
}

template <typename T>
void Space<T>::vectorComputeAccelerationBetweenBodies(const unsigned long iBody, const unsigned long jBody, const int vecDim)
{
        T vecX[dim], vecY[dim], vecZ[dim], vecLen[dim], acc[dim], accX[dim], accY[dim], accZ[dim], sqrtVecLen[dim];

        int jShuff[dim];
        for(int j=0; j<vecDim; j++)
        {
                for(int i=0; i<vecDim; i++)
                {
                jShuff[i]     = jBody+(j+i)%vecDim;

                vecX[i]       = this->positions.x[jShuff[i]] - this->positions.x[iBody+i]; // 1 flop
                vecY[i]       = this->positions.y[jShuff[i]] - this->positions.y[iBody+i]; // 1 flop
                vecZ[i]       = this->positions.z[jShuff[i]] - this->positions.z[iBody+i]; // 1 flop
                vecLen[i]     = vecX[i] * vecX[i] + vecY[i] * vecY[i] + vecZ[i]*vecZ[i];   // 5 flop                    
                sqrtVecLen[i] = sqrt(vecLen[i]);

                acc[i]  = this->masses[jShuff[i]] / (vecLen[i] * sqrtVecLen[i]); // 2 flop
                accX[i] = acc[i] * vecX[i];                                      // 1 flop
                accY[i] = acc[i] * vecY[i];                                      // 1 flop
                accZ[i] = acc[i] * vecZ[i];                                      // 1 flop

                this->accelerations.x[iBody+i] += accX[i]; // 1 flop
                this->accelerations.y[iBody+i] += accY[i]; // 1 flop
                this->accelerations.z[iBody+i] += accZ[i]; // 1 flop

                this->closestNeighborLen[iBody+i] = std::min(sqrtVecLen[i],this->closestNeighborLen[iBody+i]);
                }
        }
}

template <typename T>
void Space<T>::intrinComputeBodiesAcceleration()
{
        for(unsigned long iBody = 0; iBody < this->nBodies; iBody+=dim)
        {
		vec px, py, pz, accx, accy, accz, closest;

		px       =  vec_load(&(this->positions.x[iBody]));
		py       =  vec_load(&(this->positions.y[iBody]));
		pz       =  vec_load(&(this->positions.z[iBody]));

		accx     =  vec_load(&(this->accelerations.x[iBody]));
		accy     =  vec_load(&(this->accelerations.y[iBody]));
		accz     =  vec_load(&(this->accelerations.z[iBody]));
		
		closest  =  vec_load(&(this->closestNeighborLen[iBody]));
		
	         for(unsigned long jBody = 0; jBody < this->nBodies; jBody+=dim)
		 {
                        if(iBody != jBody)
                                this->intrinComputeAccelerationBetweenBodies(iBody, jBody, dim, px, py, pz, &accx, &accy, &accz, &closest);
                        else
                                this->selfIntrinComputeAccelerationBetweenBodies(iBody, dim, px, py, pz, &accx, &accy, &accz, &closest);
		 }        

		vec_store(&(this->accelerations.x[iBody]), accx);
		vec_store(&(this->accelerations.y[iBody]), accy);
		vec_store(&(this->accelerations.z[iBody]), accz);

		vec_store(&(this->closestNeighborLen[iBody]), closest);

         }
}

template <typename T>
void Space<T>::intrinComputeAccelerationBetweenBodies(const unsigned long iBody, const unsigned long jBody, const int vecDim, 
						           vec px, vec py, vec pz, vec *accx, vec *accy, vec *accz, vec *closest)
{
	vec rpx, rpy, rpz, masses;
        
	vec vecX, vecY, vecZ, vecLen, acc, sqrtVecLen;

        rpx      =  vec_load(&(this->positions.x[jBody]));
        rpy      =  vec_load(&(this->positions.y[jBody]));
        rpz      =  vec_load(&(this->positions.z[jBody]));

        masses   =  vec_load(&(this->masses[jBody]));

        for(int i=0; i<vecDim; i++)
        {

        vecX       = vec_sub(rpx,px);
        vecY       = vec_sub(rpy,py); 
        vecZ       = vec_sub(rpz,pz); 
        vecLen     = vec_add( vec_mul(vecZ,vecZ), vec_add( vec_mul(vecX,vecX) , vec_mul(vecY,vecY) )); 

        sqrtVecLen = vec_sqrt(vecLen);

        acc  = vec_div( masses , vec_mul(vecLen,sqrtVecLen) ); 

        *accx = vec_fmadd(acc, vecX, *accx); 
        *accy = vec_fmadd(acc, vecY, *accy); 
        *accz = vec_fmadd(acc, vecZ, *accz); 

        *closest = vec_min(sqrtVecLen,*closest);


		/*_MM_SHUFFLE(z, y, x, w)*/
		/* expands to the following value */
		/*(z<<6) | (y<<4) | (x<<2) | w*/
        
	//avant: [ a , b , c , d ]
        rpx = vec_permute(rpx,_MM_SHUFFLE(0,3,2,1));
	//après: [ b , c , d , a ]
        rpy = vec_permute(rpy,_MM_SHUFFLE(0,3,2,1));
        rpz = vec_permute(rpz,_MM_SHUFFLE(0,3,2,1));

        masses = vec_permute(masses,_MM_SHUFFLE(0,3,2,1));
        }
}

template <typename T>
void Space<T>::selfIntrinComputeAccelerationBetweenBodies(const unsigned long iBody, const int vecDim, 
						           vec px, vec py, vec pz, vec *accx, vec *accy, vec *accz, vec *closest)
{
	vec rpx, rpy, rpz, masses;
        
	vec vecX, vecY, vecZ, vecLen, acc, sqrtVecLen;

	T result[8];

        rpx      =  vec_load(&(this->positions.x[iBody]));
        rpy      =  vec_load(&(this->positions.y[iBody]));
        rpz      =  vec_load(&(this->positions.z[iBody]));

        masses   =  vec_load(&(this->masses[iBody]));

        rpx = vec_permute(rpx,_MM_SHUFFLE(0,3,2,1));
	rpy = vec_permute(rpy,_MM_SHUFFLE(0,3,2,1));
        rpz = vec_permute(rpz,_MM_SHUFFLE(0,3,2,1));

        masses = vec_permute(masses,_MM_SHUFFLE(0,3,2,1));

        for(int i=0; i<vecDim-1; i++)
        {

        
	vecX       = vec_sub(rpx,px);
        vecY       = vec_sub(rpy,py); 
        vecZ       = vec_sub(rpz,pz);
        
	vecLen     = vec_add( vec_mul(vecZ,vecZ), vec_add( vec_mul(vecX,vecX) , vec_mul(vecY,vecY) )); 
        sqrtVecLen = vec_sqrt(vecLen);

        acc  = vec_div( masses , vec_mul(vecLen,sqrtVecLen) ); 
        
        *accx = vec_fmadd(acc, vecX, *accx); 
        *accy = vec_fmadd(acc, vecY, *accy); 
        *accz = vec_fmadd(acc, vecZ, *accz); 

        *closest = vec_min(sqrtVecLen,*closest);

	//avant: [ a , b , c , d ]
        rpx = vec_permute(rpx,_MM_SHUFFLE(0,3,2,1));
	//après: [ b , c , d , a ]
        
	rpy = vec_permute(rpy,_MM_SHUFFLE(0,3,2,1));
        rpz = vec_permute(rpz,_MM_SHUFFLE(0,3,2,1));

        masses = vec_permute(masses,_MM_SHUFFLE(0,3,2,1));
        }
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
