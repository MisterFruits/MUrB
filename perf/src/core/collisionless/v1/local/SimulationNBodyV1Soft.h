/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SIMULATION_N_BODY_V1_SOFT_H_
#define SIMULATION_N_BODY_V1_SOFT_H_

#include <string>

#include "../../SimulationNBodyLocal.h"

template <typename T = double>
class SimulationNBodyV1Soft : public SimulationNBodyLocal<T>
{
private:
	T softeningSquared;
public:
	SimulationNBodyV1Soft(const unsigned long nBodies, T softening = 0.035);
	SimulationNBodyV1Soft(const std::string inputFileName, T softening = 0.035);
	virtual ~SimulationNBodyV1Soft();

protected:
	virtual void initIteration();
	virtual void computeLocalBodiesAcceleration();

	inline void computeAccelerationBetweenTwoBodiesNaive(const unsigned long &iBody, const unsigned long &jBody);
	inline void computeAccelerationBetweenTwoBodies     (const unsigned long &iBody, const unsigned long &jBody);

private:
	void init();
};

#include "SimulationNBodyV1Soft.hxx"

#endif /* SIMULATION_N_BODY_V1_SOFT_H_ */
