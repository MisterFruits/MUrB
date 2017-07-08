/*!
 * \file    SimulationNBodyV3Intrinsics.h
 * \brief   Implementation of SimulationNBodyLocal with intrinsic function calls (n² computations).
 * \author  A. Cassagne
 * \date    2014
 *
 * \section LICENSE
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode).
 */
#ifndef SIMULATION_N_BODY_V3_INTRINSICS_H_
#define SIMULATION_N_BODY_V3_INTRINSICS_H_

#include <string>
#include <mipp.h>

#include "SimulationNBodyV3.h"

/*!
 * \class  SimulationNBodyV3Intrinsics
 * \brief  Implementation of SimulationNBodyLocal with intrinsic function calls (n² computations).
 *
 * \tparam T : Type.
 */
template <typename T = double>
class SimulationNBodyV3Intrinsics : public SimulationNBodyV3<T>
{
public:
	SimulationNBodyV3Intrinsics(const unsigned long nBodies, T softening = 0.035);
	SimulationNBodyV3Intrinsics(const std::string inputFileName, T softening = 0.035);
	virtual ~SimulationNBodyV3Intrinsics();

protected:
	virtual void initIteration();
	virtual void computeLocalBodiesAcceleration();

	static inline void computeAccelerationBetweenTwoBodies(const mipp::reg &rG,
	                                                       const mipp::reg &rSoftSquared,
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

#include "SimulationNBodyV3Intrinsics.hxx"

#endif /* SIMULATION_N_BODY_V3_INTRINSICS_H_ */
