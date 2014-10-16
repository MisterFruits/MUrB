/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SIMULATION_N_BODY_LOCAL_H_
#define SIMULATION_N_BODY_LOCAL_H_

#include <string>

#include "Bodies.h"

#include "SimulationNBody.h"

template <typename T = double>
class SimulationNBodyLocal : public SimulationNBody<T>
{
protected:
	SimulationNBodyLocal(const unsigned long nBodies);
	SimulationNBodyLocal(const std::string inputFileName);

public:
	virtual ~SimulationNBodyLocal();

	void computeOneIteration();

protected:
	virtual void initIteration()                  = 0;
	virtual void computeLocalBodiesAcceleration() = 0;

private:
	void findTimeStep();
};

#include "SimulationNBodyLocal.hxx"

#endif /* SIMULATION_N_BODY_LOCAL_H_ */
