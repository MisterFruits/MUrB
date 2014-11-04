/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SIMULATION_N_BODY_MPI_V1_INTRINSICS_H_
#define SIMULATION_N_BODY_MPI_V1_INTRINSICS_H_

#include <string>

#include "../../../../utils/myIntrinsicsPlusPlus.h"

#include "../../SimulationNBodyMPI.h"

template <typename T = double>
class SimulationNBodyMPIV1Intrinsics : public SimulationNBodyMPI<T>
{
public:
	SimulationNBodyMPIV1Intrinsics(const unsigned long nBodies);
	SimulationNBodyMPIV1Intrinsics(const std::string inputFileName);
	virtual ~SimulationNBodyMPIV1Intrinsics();

protected:
	virtual void initIteration();
	virtual void computeLocalBodiesAcceleration();
	virtual void computeNeighborBodiesAcceleration();

	static inline void computeAccelerationBetweenTwoBodies(const mipp::vec &rG,
	                                                       const mipp::vec &rIPosX,
	                                                       const mipp::vec &rIPosY,
	                                                       const mipp::vec &rIPosZ,
	                                                             mipp::vec &rIAccX,
	                                                             mipp::vec &rIAccY,
	                                                             mipp::vec &rIAccZ,
	                                                             mipp::vec &rIClosNeiDist,
	                                                       const mipp::vec &rJMass,
	                                                       const mipp::vec &rJPosX,
	                                                       const mipp::vec &rJPosY,
	                                                       const mipp::vec &rJPosZ);
private:
	void init();
	void _computeLocalBodiesAcceleration();
	void _computeNeighborBodiesAcceleration();
};

#include "SimulationNBodyMPIV1Intrinsics.hxx"

#endif /* SIMULATION_N_BODY_MPI_V1_INTRINSICS_H_ */
