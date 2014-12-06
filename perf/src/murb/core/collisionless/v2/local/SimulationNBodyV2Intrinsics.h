/*!
 * \file    SimulationNBodyV2Intrinsics.h
 * \brief   Implementation of SimulationNBodyLocal with intrinsic function calls (n²/2 computations).
 * \author  G. Hautreux
 * \date    2014
 *
 * \section LICENSE
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode).
 */
#ifndef SIMULATION_N_BODY_V2_INTRINSICS
#define SIMULATION_N_BODY_V2_INTRINSICS

#include <string>

#include "../../../../../common/utils/mipp.h"

#include "SimulationNBodyV2.h"

/*!
 * \class  SimulationNBodyV2Intrinsics
 * \brief  Implementation of SimulationNBodyLocal with intrinsic function calls (n²/2 computations).
 *
 * \tparam T : Type.
 */
template <typename T = double>
class SimulationNBodyV2Intrinsics : public SimulationNBodyV2<T>
{
public:
	SimulationNBodyV2Intrinsics(const unsigned long nBodies);
	SimulationNBodyV2Intrinsics(const std::string inputFileName);
	virtual ~SimulationNBodyV2Intrinsics();

protected:
	virtual void initIteration();
	virtual void computeLocalBodiesAcceleration();

	static inline void computeAccelerationBetweenTwoBodies(const mipp::vec &rG,
	                                                       const mipp::vec &rmi,
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
	                                                       const mipp::vec &rqjZ,
	                                                             mipp::vec &rajX,
	                                                             mipp::vec &rajY,
	                                                             mipp::vec &rajZ,
	                                                             mipp::vec &rclosNeighj);
private:
	void _computeLocalBodiesAcceleration();
	void _reAllocateBuffers();
	void reAllocateBuffers();
	void _initIteration();
};

#include "SimulationNBodyV2Intrinsics.hxx"

#endif /* SIMULATION_N_BODY_V2_INTRINSICS */
