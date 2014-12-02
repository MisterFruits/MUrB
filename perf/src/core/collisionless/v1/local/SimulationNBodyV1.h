/*!
 * \file    SimulationNBodyV1.h
 * \brief   Naive implementation of SimulationNBody.
 * \author  A. Cassagne
 * \date    2014
 *
 * \section LICENSE
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode).
 */

#ifndef SIMULATION_N_BODY_V1_H_
#define SIMULATION_N_BODY_V1_H_

#include <string>

#include "../../SimulationNBodyLocal.h"

/*!
 * \class  SimulationNBodyV1
 * \brief  Naive implementation of SimulationNBody (n² computations).
 *
 * \tparam T : Type.
 */
template <typename T = double>
class SimulationNBodyV1 : public SimulationNBodyLocal<T>
{
public:
	SimulationNBodyV1(const unsigned long nBodies);
	SimulationNBodyV1(const std::string inputFileName);
	virtual ~SimulationNBodyV1();

protected:
	virtual void initIteration();
	virtual void computeLocalBodiesAcceleration();

	static inline void computeAccelerationBetweenTwoBodies(const T &G,
	                                                       const T &qiX, const T &qiY, const T &qiZ,
	                                                             T &aiX,       T &aiY,       T &aiZ,
	                                                             T &closNeighi,
	                                                       const T &mj,
	                                                       const T &qjX, const T &qjY, const T &qjZ);
private:
	void init();
};

#include "SimulationNBodyV1.hxx"

#endif /* SIMULATION_N_BODY_V1_H_ */
