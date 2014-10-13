/*
 * Do not remove.
 * Gabriel Hautreux, gabriel.hautreux@gmail.com
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

#include "SimulationNBodyV2CB.h"

template <typename T>
SimulationNBodyV2CB<T>::SimulationNBodyV2CB(const unsigned long nBodies)
	: SimulationNBodyV2<T>(nBodies)
{
}

template <typename T>
SimulationNBodyV2CB<T>::SimulationNBodyV2CB(const std::string inputFileName)
	: SimulationNBodyV2<T>(inputFileName)
{
}

template <typename T>
SimulationNBodyV2CB<T>::~SimulationNBodyV2CB()
{
}

template <typename T>
void SimulationNBodyV2CB<T>::computeBodiesAcceleration()
{
	unsigned long blockSize = 512;
	for(unsigned long jOff = 0; jOff < this->bodies.getN(); jOff += blockSize)
	{
		blockSize = std::min(blockSize, this->bodies.getN() - jOff);
#pragma omp parallel
{
		const unsigned tid = omp_get_thread_num();

#pragma omp for schedule(runtime)
		for(unsigned long iBody = jOff + 1 ; iBody < this->bodies.getN(); iBody++)
			for(unsigned long jBody = jOff ; jBody < jOff + blockSize; jBody++)
				if(iBody > jBody)
					//this->computeAccelerationBetweenTwoBodiesNaive(iBody, jBody, tid);
					this->computeAccelerationBetweenTwoBodies(iBody, jBody, tid);
}
	}

	if(this->nMaxThreads > 1)
	{
		for(unsigned long iBody = 0; iBody < this->bodies.getN(); iBody++)
			for(unsigned iThread = 1; iThread < this->nMaxThreads; iThread++)
			{
				this->accelerations.x[iBody] += this->accelerations.x[iBody + iThread * this->bodies.getN()];
				this->accelerations.y[iBody] += this->accelerations.y[iBody + iThread * this->bodies.getN()];
				this->accelerations.z[iBody] += this->accelerations.z[iBody + iThread * this->bodies.getN()];
			}
	}

}
