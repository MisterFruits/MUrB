/*!
 * \file    SimulationNBodyMPIV1Intrinsics.h
 * \brief   Implementation of SimulationNBodyMPI with intrinsic function calls (n² computations).
 * \author  A. Cassagne
 * \date    2014
 *
 * \section LICENSE
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode).
 */
#ifndef SIMULATION_N_BODY_MPI_V1_INTRINSICS_H_
#define SIMULATION_N_BODY_MPI_V1_INTRINSICS_H_

#include <string>
#include <mipp.h>

#include "../../SimulationNBodyMPI.h"

/*!
 * \class  SimulationNBodyMPIV1Intrinsics
 * \brief  Implementation of SimulationNBodyMPI with intrinsic function calls (n² computations).
 *
 * \tparam T : Type.
 */
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

	static inline void computeAccelerationBetweenTwoBodies(const mipp::reg &rG,
	                                                       const mipp::reg &rqiX,
	                                                       const mipp::reg &rqiY,
	                                                       const mipp::reg &rqiZ,
	                                                             mipp::reg &raiX,
	                                                             mipp::reg &raiY,
	                                                             mipp::reg &raiZ,
	                                                             mipp::reg &rclosNeighi,
	                                                       const mipp::reg &rmj,
	                                                       const mipp::reg &rqjX,
	                                                       const mipp::reg &rqjY,
	                                                       const mipp::reg &rqjZ);
private:
	void init();
	void _computeLocalBodiesAcceleration();
	void _computeNeighborBodiesAcceleration();
};

#include "SimulationNBodyMPIV1Intrinsics.hxx"

#endif /* SIMULATION_N_BODY_MPI_V1_INTRINSICS_H_ */
