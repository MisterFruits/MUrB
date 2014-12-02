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

#include "../../../../utils/myIntrinsicsPlusPlus.h"

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

	static inline void computeAccelerationBetweenTwoBodies(const mipp::vec &rG,
	                                                       const mipp::vec &rqiX,
	                                                       const mipp::vec &rqiY,
	                                                       const mipp::vec &rqiZ,
	                                                             mipp::vec &raiX,
	                                                             mipp::vec &raiY,
	                                                             mipp::vec &raiZ,
	                                                             mipp::vec &rclosNeighi,
	                                                       const mipp::vec &rmj,
	                                                       const mipp::vec &rqjX,
	                                                       const mipp::vec &rqjY,
	                                                       const mipp::vec &rqjZ);
private:
	void init();
	void _computeLocalBodiesAcceleration();
	void _computeNeighborBodiesAcceleration();
};

#include "SimulationNBodyMPIV1Intrinsics.hxx"

#endif /* SIMULATION_N_BODY_MPI_V1_INTRINSICS_H_ */
