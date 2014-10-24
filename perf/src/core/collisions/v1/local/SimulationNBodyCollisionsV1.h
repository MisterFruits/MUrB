/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SIMULATION_N_BODY_COLLISIONS_V1_H_
#define SIMULATION_N_BODY_COLLISIONS_V1_H_

#include <string>

#include "../../SimulationNBodyCollisionsLocal.h"

template <typename T = double>
class SimulationNBodyCollisionsV1 : public SimulationNBodyCollisionsLocal<T>
{
public:
	SimulationNBodyCollisionsV1(const unsigned long nBodies);
	SimulationNBodyCollisionsV1(const std::string inputFileName);
	virtual ~SimulationNBodyCollisionsV1();

protected:
	virtual void initIteration();
	virtual void computeLocalBodiesAcceleration();

	inline void computeAccelerationBetweenTwoBodiesNaive(const unsigned long &iBody, const unsigned long &jBody);
	inline void computeAccelerationBetweenTwoBodies     (const unsigned long &iBody, const unsigned long &jBody);

private:
	void init();
};

#include "SimulationNBodyCollisionsV1.hxx"

#endif /* SIMULATION_N_BODY_COLLISIONS_V1_H_ */
