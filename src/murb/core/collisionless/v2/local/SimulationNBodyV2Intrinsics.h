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
#include <mipp.h>

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

	static inline void computeAccelerationBetweenTwoBodies(const mipp::reg &rG,
	                                                       const mipp::reg &rmi,
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
	                                                       const mipp::reg &rqjZ,
	                                                             mipp::reg &rajX,
	                                                             mipp::reg &rajY,
	                                                             mipp::reg &rajZ,
	                                                             mipp::reg &rclosNeighj);
private:
	void _computeLocalBodiesAcceleration();
	void _reAllocateBuffers();
	void reAllocateBuffers();
	void _initIteration();
};

#include "SimulationNBodyV2Intrinsics.hxx"

#endif /* SIMULATION_N_BODY_V2_INTRINSICS */
