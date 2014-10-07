/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SIMULATION_N_BODY_V2_H_
#define SIMULATION_N_BODY_V2_H_

#include <string>

#include "SimulationNBody.h"

template <typename T = double>
class SimulationNBodyV2 : public SimulationNBody<T>
{
public:
	SimulationNBodyV2(const unsigned long nBodies);
	SimulationNBodyV2(const std::string inputFileName);
	virtual ~SimulationNBodyV2();

protected:
	virtual void allocateBuffers();
	virtual void initIteration();
	virtual void computeBodiesAcceleration();

	inline void computeAccelerationBetweenTwoBodiesNaive(const unsigned long iBody, const unsigned long jBody);
	inline void computeAccelerationBetweenTwoBodies     (const unsigned long iBody, const unsigned long jBody);
};

#include "SimulationNBodyV2.hxx"

#endif /* SIMULATION_N_BODY_V2_H_ */
