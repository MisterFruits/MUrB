/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SIMULATION_N_BODY_V1_H_
#define SIMULATION_N_BODY_V1_H_

#include <string>

#include "../SimulationNBody.h"

template <typename T = double>
class SimulationNBodyV1 : public SimulationNBody<T>
{
public:
	SimulationNBodyV1(const unsigned long nBodies);
	SimulationNBodyV1(const std::string inputFileName);
	virtual ~SimulationNBodyV1();

protected:
	virtual void initIteration();
	virtual void computeBodiesAcceleration();

	inline void computeAccelerationBetweenTwoBodiesNaive(const unsigned long &iBody, const unsigned long &jBody);
	inline void computeAccelerationBetweenTwoBodies     (const unsigned long &iBody, const unsigned long &jBody);

private:
	void init();
};

#include "SimulationNBodyV1.hxx"

#endif /* SIMULATION_N_BODY_V1_H_ */
