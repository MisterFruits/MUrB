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

#ifdef _OPENMP
#include <omp.h>
#else
#ifndef _MYOPENMP
#define _MYOPENMP
inline void omp_set_num_threads(int) {           }
inline int  omp_get_num_threads(   ) { return 1; }
inline int  omp_get_max_threads(   ) { return 1; }
inline int  omp_get_thread_num (   ) { return 0; }
#endif
#endif

#include "SimulationNBodyV1CB.h"

template <typename T>
SimulationNBodyV1CB<T>::SimulationNBodyV1CB(const unsigned long nBodies)
	: SimulationNBodyV1<T>(nBodies)
{
}

template <typename T>
SimulationNBodyV1CB<T>::SimulationNBodyV1CB(const std::string inputFileName)
	: SimulationNBodyV1<T>(inputFileName)
{
}

template <typename T>
SimulationNBodyV1CB<T>::~SimulationNBodyV1CB()
{
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
void SimulationNBodyV1CB<T>::computeBodiesAcceleration()
{
	unsigned long blockSize = 512;
	// flops  = 23 * blockSize * nBodies      * nBlocks
	// memops = (4 * blockSize + 7 * nBodies) * nBlocks
	for(unsigned long jOff = 0; jOff < this->bodies.getN(); jOff += blockSize)
	{
		blockSize = std::min(blockSize, this->bodies.getN() - jOff);
		// flops  = 23 * blockSize * nBodies
		// memops =  4 * blockSize + 7 * nBodies
#pragma omp parallel for schedule(runtime)
		for(unsigned long iBody = 0; iBody < this->bodies.getN(); iBody++)
			// flops  = 23 * blockSize
			// memops =  4 * blockSize + 7
			for(unsigned long jBody = jOff; jBody < jOff + blockSize; jBody++)
				if(iBody != jBody)
					//this->computeAccelerationBetweenTwoBodiesNaive(iBody, jBody);
					this->computeAccelerationBetweenTwoBodies(iBody, jBody);
	}
}
