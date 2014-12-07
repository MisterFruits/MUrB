/*!
 * \file    SimulationNBodyV3Intrinsics.h
 * \brief   Implementation of SimulationNBodyLocal with intrinsic function calls (n² computations).
 * \author  A. Cassagne
 * \date    2014
 *
 * \section LICENSE
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode).
 */
#ifndef SIMULATION_N_BODY_V3_INTRINSICS_BH_H_
#define SIMULATION_N_BODY_V3_INTRINSICS_BH_H_

#include <string>

#include "../../../../../common/utils/mipp.h"

#include "SimulationNBodyV3Intrinsics.h"

/*!
 * \class  SimulationNBodyV3IntrinsicsBH
 * \brief  Implementation of SimulationNBodyLocal with intrinsic function calls (n² computations).
 *
 * \tparam T : Type.
 */
template <typename T = double>
class SimulationNBodyV3IntrinsicsBH : public SimulationNBodyV3Intrinsics<T>
{
protected:
	T mbh  = 4000;   /*!< black hole mass */
	T rbh  = 0.5e8;  /*!< black hole radius */
	T qbhX = 0;      /*!< black hole position x */
	T qbhY = 0;      /*!< black hole position y */
	T qbhZ = 0;      /*!< black hole position z */

public:
	SimulationNBodyV3IntrinsicsBH(const unsigned long nBodies, T softening = 0.035);
	SimulationNBodyV3IntrinsicsBH(const std::string inputFileName, T softening = 0.035);
	virtual ~SimulationNBodyV3IntrinsicsBH();

protected:
	virtual void initIteration();
	virtual void computeLocalBodiesAcceleration();
	virtual void computeLocalBodiesAccelerationWithBlackHole();

public:
	static inline T computeAccelerationBetweenBodyAndBlackHole(const T &G,   const T &softSquared,
	                                                           const T &qiX, const T &qiY, const T &qiZ,
	                                                                 T &aiX,       T &aiY,       T &aiZ,
	                                                           const T &mbh,
	                                                           const T &rbh,
	                                                           const T &qbhX, const T &qbhY, const T &qbhZ);

private:
	void init();
	void _computeLocalBodiesAcceleration();
};

#include "SimulationNBodyV3IntrinsicsBH.hxx"

#endif /* SIMULATION_N_BODY_V3_INTRINSICS_BH_H_ */
