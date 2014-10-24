/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SIMULATION_N_BODY_COLLISIONS_LOCAL_H_
#define SIMULATION_N_BODY_COLLISIONS_LOCAL_H_

#include <string>

#include "../SimulationNBody.h"

template <typename T = double>
class SimulationNBodyCollisionsLocal : public SimulationNBody<T>
{
protected:
	std::vector<std::vector<unsigned long>> collisions;

	SimulationNBodyCollisionsLocal(const unsigned long nBodies, const unsigned long randInit = 1);
	SimulationNBodyCollisionsLocal(const std::string inputFileName);

public:
	virtual ~SimulationNBodyCollisionsLocal();

	void computeOneIteration();

protected:
	virtual void initIteration()                  = 0;
	virtual void computeLocalBodiesAcceleration() = 0;

private:
	void findTimeStep();
};

#include "SimulationNBodyCollisionsLocal.hxx"

#endif /* SIMULATION_N_BODY_COLLISIONS_LOCAL_H_ */
