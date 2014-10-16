/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SIMULATION_N_BODY_V2_H_
#define SIMULATION_N_BODY_V2_H_

#include <string>

#include "../SimulationNBodyLocal.h"

template <typename T = double>
class SimulationNBodyV2 : public SimulationNBodyLocal<T>
{
public:
	SimulationNBodyV2(const unsigned long nBodies);
	SimulationNBodyV2(const std::string inputFileName);
	virtual ~SimulationNBodyV2();

protected:
	virtual void initIteration();
	virtual void computeLocalBodiesAcceleration();

	inline void computeAccelerationBetweenTwoBodiesNaive(const unsigned long &iBody,
	                                                     const unsigned long &jBody,
	                                                     const unsigned      &tid);
	inline void computeAccelerationBetweenTwoBodies(const unsigned long &iBody,
	                                                const unsigned long &jBody,
	                                                const unsigned      &tid);

private:
	void reAllocateBuffers();
};

#include "SimulationNBodyV2.hxx"

#endif /* SIMULATION_N_BODY_V2_H_ */
