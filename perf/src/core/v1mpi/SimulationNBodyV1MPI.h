/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SIMULATION_N_BODY_V1_MPI_H_
#define SIMULATION_N_BODY_V1_MPI_H_

#include <string>

#include "../SimulationNBodyMPI.h"

template <typename T = double>
class SimulationNBodyV1MPI : public SimulationNBodyMPI<T>
{
public:
	SimulationNBodyV1MPI(const unsigned long nBodies);
	SimulationNBodyV1MPI(const std::string inputFileName);
	virtual ~SimulationNBodyV1MPI();

protected:
	virtual void initIteration();
	virtual void computeLocalBodiesAcceleration();
	virtual void computeNeighborBodiesAcceleration();

	inline void computeAccelerationBetweenTwoBodiesNaive(const T &iMasses,
	                                                     const T &iPosX, const T &iPosY, const T &iPosZ,
	                                                           T &iAccsX,      T &iAccsY,      T &iAccsZ,
	                                                           T &iClosNeiDist,
	                                                     const T &jMasses,
	                                                     const T &jPosX, const T &jPosY, const T &jPosZ);

	inline void computeAccelerationBetweenTwoBodies(const T &iPosX, const T &iPosY, const T &iPosZ,
	                                                      T &iAccsX,      T &iAccsY,      T &iAccsZ,
	                                                      T &iClosNeiDist,
	                                                const T &jMasses,
	                                                const T &jPosX, const T &jPosY, const T &jPosZ);

private:
	void init();
};

#include "SimulationNBodyV1MPI.hxx"

#endif /* SIMULATION_N_BODY_V1_MPI_H_ */
