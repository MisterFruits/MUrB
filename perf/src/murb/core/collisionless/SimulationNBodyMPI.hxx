/*!
 * \file    SimulationNBodyMPI.hxx
 * \brief   Abstract n-body collision less simulation class for MPI computations (between many nodes).
 * \author  A. Cassagne
 * \date    2014
 *
 * \section LICENSE
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode).
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
#ifndef NO_OMP
#define NO_OMP
inline void omp_set_num_threads(int) {           }
inline int  omp_get_num_threads(   ) { return 1; }
inline int  omp_get_max_threads(   ) { return 1; }
inline int  omp_get_thread_num (   ) { return 0; }
#endif
#endif

#include "../../utils/myIntrinsicsPlusPlus.h"
#include "../../utils/ToMPIDatatype.h"

#include "SimulationNBodyMPI.h"

template <typename T>
SimulationNBodyMPI<T>::SimulationNBodyMPI(const unsigned long nBodies)
	: SimulationNBody<T>(new Bodies<T>(nBodies, MPI::COMM_WORLD.Get_rank() * nBodies +1)),
	  neighborBodies    (NULL),
	  MPISize           (MPI::COMM_WORLD.Get_size()),
	  MPIRank           (MPI::COMM_WORLD.Get_rank()),
	  MPIBodiesBuffers  {this->bodies->getN(), this->bodies->getN()}
{
	this->init();
}

template <typename T>
SimulationNBodyMPI<T>::SimulationNBodyMPI(const std::string inputFileName)
	: SimulationNBody<T>(new Bodies<T>(inputFileName)),
	  neighborBodies    (NULL),
	  MPISize           (MPI::COMM_WORLD.Get_size()),
	  MPIRank           (MPI::COMM_WORLD.Get_rank()),
	  MPIBodiesBuffers  {this->bodies->getN(), this->bodies->getN()}
{
	this->init();
}

template <typename T>
void SimulationNBodyMPI<T>::init()
{
	const unsigned long n = this->bodies->getN() + this->bodies->getPadding();
	MPI::Datatype MPIType = ToMPIDatatype<T>::value();

	// who is next or previous MPI process ?
	int MPIPrevRank = (this->MPIRank == 0) ? this->MPISize -1 : this->MPIRank -1;
	int MPINextRank = (this->MPIRank +1) % this->MPISize;

	Bodies<T> *bBuffs[] = {&this->MPIBodiesBuffers[0], &this->MPIBodiesBuffers[1]};

	// send number of bodies
	this->MPIRequests[0][0] = MPI::COMM_WORLD.Send_init(&bBuffs[0]->n,            1, MPI_LONG, MPINextRank, 0);
	this->MPIRequests[1][0] = MPI::COMM_WORLD.Send_init(&bBuffs[1]->n,            1, MPI_LONG, MPINextRank, 1);
	// receive number of bodies
	this->MPIRequests[0][1] = MPI::COMM_WORLD.Recv_init(&bBuffs[1]->n,            1, MPI_LONG, MPIPrevRank, 0);
	this->MPIRequests[1][1] = MPI::COMM_WORLD.Recv_init(&bBuffs[0]->n,            1, MPI_LONG, MPIPrevRank, 1);

	// send masses
	this->MPIRequests[0][2] = MPI::COMM_WORLD.Send_init(bBuffs[0]->masses,        n, MPIType,  MPINextRank, 2);
	this->MPIRequests[1][2] = MPI::COMM_WORLD.Send_init(bBuffs[1]->masses,        n, MPIType,  MPINextRank, 3);
	// receive masses
	this->MPIRequests[0][3] = MPI::COMM_WORLD.Recv_init(bBuffs[1]->masses,        n, MPIType,  MPIPrevRank, 2);
	this->MPIRequests[1][3] = MPI::COMM_WORLD.Recv_init(bBuffs[0]->masses,        n, MPIType,  MPIPrevRank, 3);

	// send radiuses
	this->MPIRequests[0][4] = MPI::COMM_WORLD.Send_init(bBuffs[0]->radiuses,      n, MPIType,  MPINextRank, 4);
	this->MPIRequests[1][4] = MPI::COMM_WORLD.Send_init(bBuffs[1]->radiuses,      n, MPIType,  MPINextRank, 5);
	// receive radiuses
	this->MPIRequests[0][5] = MPI::COMM_WORLD.Recv_init(bBuffs[1]->radiuses,      n, MPIType,  MPIPrevRank, 4);
	this->MPIRequests[1][5] = MPI::COMM_WORLD.Recv_init(bBuffs[0]->radiuses,      n, MPIType,  MPIPrevRank, 5);

	// send positions.x
	this->MPIRequests[0][6] = MPI::COMM_WORLD.Send_init(bBuffs[0]->positions.x,   n, MPIType,  MPINextRank, 6);
	this->MPIRequests[1][6] = MPI::COMM_WORLD.Send_init(bBuffs[1]->positions.x,   n, MPIType,  MPINextRank, 7);
	// receive positions.x
	this->MPIRequests[0][7] = MPI::COMM_WORLD.Recv_init(bBuffs[1]->positions.x,   n, MPIType,  MPIPrevRank, 6);
	this->MPIRequests[1][7] = MPI::COMM_WORLD.Recv_init(bBuffs[0]->positions.x,   n, MPIType,  MPIPrevRank, 7);

	// send positions.y
	this->MPIRequests[0][8] = MPI::COMM_WORLD.Send_init(bBuffs[0]->positions.y,   n, MPIType,  MPINextRank, 8);
	this->MPIRequests[1][8] = MPI::COMM_WORLD.Send_init(bBuffs[1]->positions.y,   n, MPIType,  MPINextRank, 9);
	// receive positions.y
	this->MPIRequests[0][9] = MPI::COMM_WORLD.Recv_init(bBuffs[1]->positions.y,   n, MPIType,  MPIPrevRank, 8);
	this->MPIRequests[1][9] = MPI::COMM_WORLD.Recv_init(bBuffs[0]->positions.y,   n, MPIType,  MPIPrevRank, 9);

	// send positions.z
	this->MPIRequests[0][10] = MPI::COMM_WORLD.Send_init(bBuffs[0]->positions.z,  n, MPIType, MPINextRank, 10);
	this->MPIRequests[1][10] = MPI::COMM_WORLD.Send_init(bBuffs[1]->positions.z,  n, MPIType, MPINextRank, 11);
	// receive positions.z
	this->MPIRequests[0][11] = MPI::COMM_WORLD.Recv_init(bBuffs[1]->positions.z,  n, MPIType, MPIPrevRank, 10);
	this->MPIRequests[1][11] = MPI::COMM_WORLD.Recv_init(bBuffs[0]->positions.z,  n, MPIType, MPIPrevRank, 11);

	// send velocities.x
	this->MPIRequests[0][12] = MPI::COMM_WORLD.Send_init(bBuffs[0]->velocities.x, n, MPIType, MPINextRank, 12);
	this->MPIRequests[1][12] = MPI::COMM_WORLD.Send_init(bBuffs[1]->velocities.x, n, MPIType, MPINextRank, 13);
	// receive velocities.x
	this->MPIRequests[0][13] = MPI::COMM_WORLD.Recv_init(bBuffs[1]->velocities.x, n, MPIType, MPIPrevRank, 12);
	this->MPIRequests[1][13] = MPI::COMM_WORLD.Recv_init(bBuffs[0]->velocities.x, n, MPIType, MPIPrevRank, 13);

	// send velocities.y
	this->MPIRequests[0][14] = MPI::COMM_WORLD.Send_init(bBuffs[0]->velocities.y, n, MPIType, MPINextRank, 14);
	this->MPIRequests[1][14] = MPI::COMM_WORLD.Send_init(bBuffs[1]->velocities.y, n, MPIType, MPINextRank, 15);
	// receive velocities.y
	this->MPIRequests[0][15] = MPI::COMM_WORLD.Recv_init(bBuffs[1]->velocities.y, n, MPIType, MPIPrevRank, 14);
	this->MPIRequests[1][15] = MPI::COMM_WORLD.Recv_init(bBuffs[0]->velocities.y, n, MPIType, MPIPrevRank, 15);

	// send velocities.z
	this->MPIRequests[0][16] = MPI::COMM_WORLD.Send_init(bBuffs[0]->velocities.z, n, MPIType, MPINextRank, 16);
	this->MPIRequests[1][16] = MPI::COMM_WORLD.Send_init(bBuffs[1]->velocities.z, n, MPIType, MPINextRank, 17);
	// receive velocities.z
	this->MPIRequests[0][17] = MPI::COMM_WORLD.Recv_init(bBuffs[1]->velocities.z, n, MPIType, MPIPrevRank, 16);
	this->MPIRequests[1][17] = MPI::COMM_WORLD.Recv_init(bBuffs[0]->velocities.z, n, MPIType, MPIPrevRank, 17);
}

template <typename T>
SimulationNBodyMPI<T>::~SimulationNBodyMPI()
{
	delete this->bodies;
}

template <typename T>
void SimulationNBodyMPI<T>::computeOneIteration()
{
	this->MPIBodiesBuffers[0].hardCopy(*this->bodies);

	this->initIteration();

	// compute bodies acceleration ------------------------------------------------------------------------------------
	MPI::Prequest::Startall(2*9, this->MPIRequests[0]);

	this->computeLocalBodiesAcceleration();
 
	for(int iStep = 1; iStep < this->MPISize; ++iStep)
	{
		MPI::Request::Waitall(2*9, this->MPIRequests[(iStep -1) % 2]);

		this->neighborBodies = &this->MPIBodiesBuffers[iStep % 2];

		this->computeNeighborBodiesAcceleration();

		if(iStep < this->MPISize -1)
			MPI::Prequest::Startall(2*9, this->MPIRequests[iStep % 2]);
	}
	// ----------------------------------------------------------------------------------------------------------------

	// find time step -------------------------------------------------------------------------------------------------
	if(!this->dtConstant)
		this->findTimeStep();
	// ----------------------------------------------------------------------------------------------------------------

	// update positions and velocities --------------------------------------------------------------------------------
	this->bodies->updatePositionsAndVelocities(this->accelerations, this->dt);
	// ----------------------------------------------------------------------------------------------------------------
}

template <typename T>
void SimulationNBodyMPI<T>::findTimeStep()
{
	// TODO: be careful with the V1Intrinsics version: with fake bodies added at the end of the last vector, the
	//       dynamic time step is broken.
	//       It is necessary to launch the simulation with a number of bodies multiple of mipp::vectorSize<T>()!
	if(!this->dtConstant)
	{
		T localDt = std::numeric_limits<T>::infinity();
		for(unsigned long iBody = 0; iBody < this->bodies->getN(); iBody++)
			localDt = std::min(localDt, this->computeTimeStep(iBody));

		MPI::COMM_WORLD.Allreduce(&localDt, &this->dt, 1, ToMPIDatatype<T>::value(), MPI_MIN);

		if(this->dt < this->minDt)
			this->dt = this->minDt;
	}
}
