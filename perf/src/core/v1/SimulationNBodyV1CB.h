/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SIMULATION_N_BODY_V1_CB_H_
#define SIMULATION_N_BODY_V1_CB_H_

#include <string>

#include "SimulationNBodyV1.h"

template <typename T = double>
class SimulationNBodyV1CB : public SimulationNBodyV1<T>
{
public:
	SimulationNBodyV1CB(const unsigned long nBodies);
	SimulationNBodyV1CB(const std::string inputFileName);
	virtual ~SimulationNBodyV1CB();

protected:
	virtual void computeBodiesAcceleration();
};

#include "SimulationNBodyV1CB.hxx"

#endif /* SIMULATION_N_BODY_V1_CB_H_ */
