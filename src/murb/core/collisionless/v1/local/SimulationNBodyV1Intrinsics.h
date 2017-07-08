/*!
 * \file    SimulationNBodyV1Intrinsics.h
 * \brief   Implementation of SimulationNBodyLocal with intrinsic function calls (n² computations).
 * \author  A. Cassagne
 * \date    2014
 *
 * \section LICENSE
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode).
 */
#ifndef SIMULATION_N_BODY_V1_INTRINSICS_H_
#define SIMULATION_N_BODY_V1_INTRINSICS_H_

#include <string>

#include <mipp.h>

#include "SimulationNBodyV1.h"

/*!
 * \class  SimulationNBodyV1Intrinsics
 * \brief  Implementation of SimulationNBodyLocal with intrinsic function calls (n² computations).
 *
 * \tparam T : Type.
 */
template <typename T = double>
class SimulationNBodyV1Intrinsics : public SimulationNBodyV1<T>
{
public:
	SimulationNBodyV1Intrinsics(const unsigned long nBodies);
	SimulationNBodyV1Intrinsics(const std::string inputFileName);
	virtual ~SimulationNBodyV1Intrinsics();

protected:
	virtual void initIteration();
	virtual void computeLocalBodiesAcceleration();

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
};

#include "SimulationNBodyV1Intrinsics.hxx"

#endif /* SIMULATION_N_BODY_V1_INTRINSICS_H_ */
