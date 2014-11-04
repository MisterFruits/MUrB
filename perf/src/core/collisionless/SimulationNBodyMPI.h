/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SIMULATION_N_BODY_MPI_H_
#define SIMULATION_N_BODY_MPI_H_

#include <string>
#include <mpi.h>

#include "../Bodies.h"

#include "../SimulationNBody.h"

template <typename T = double>
class SimulationNBodyMPI : public SimulationNBody<T>
{
protected:
	Bodies<T> *neighborBodies;
	const int MPISize;
	const int MPIRank;

private:
	Bodies<T>     MPIBodiesBuffers[2];
	MPI::Prequest MPIRequests[2][2*9];

protected:
	SimulationNBodyMPI(const unsigned long nBodies);
	SimulationNBodyMPI(const std::string inputFileName);

public:
	virtual ~SimulationNBodyMPI();

	void computeOneIteration();

protected:
	virtual void initIteration()                     = 0;
	virtual void computeLocalBodiesAcceleration()    = 0;
	virtual void computeNeighborBodiesAcceleration() = 0;

private:
	void init();
	void findTimeStep();
};

#include "SimulationNBodyMPI.hxx"

#endif /* SIMULATION_N_BODY_MPI_H_ */
