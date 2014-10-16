/*
 * Do not remove.
 * Gabriel Hautreux, gabriel.hautreux@gmail.com
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SIMULATION_N_BODY_V2_CB_H_
#define SIMULATION_N_BODY_V2_CB_H_

#include <string>

#include "SimulationNBodyV2.h"

template <typename T = double>
class SimulationNBodyV2CB : public SimulationNBodyV2<T>
{
public:
	SimulationNBodyV2CB(const unsigned long nBodies);
	SimulationNBodyV2CB(const std::string inputFileName);
	virtual ~SimulationNBodyV2CB();

protected:
	virtual void computeLocalBodiesAcceleration();
};

#include "SimulationNBodyV2CB.hxx"

#endif /* SIMULATION_N_BODY_V2_CB_H_ */
