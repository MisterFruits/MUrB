/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SIMULATION_N_BODY_V3_H_
#define SIMULATION_N_BODY_V3_H_

#include <string>

#include "../../SimulationNBodyLocal.h"

template <typename T = double>
class SimulationNBodyV3 : public SimulationNBodyLocal<T>
{
private:
	T softeningSquared;

public:
	SimulationNBodyV3(const unsigned long nBodies, T softening = 0.035);
	SimulationNBodyV3(const std::string inputFileName, T softening = 0.035);
	virtual ~SimulationNBodyV3();

protected:
	virtual void initIteration();
	virtual void computeLocalBodiesAcceleration();

	static inline void computeAccelerationBetweenTwoBodies(const T &G,   const T &softSquared,
	                                                       const T &qiX, const T &qiY, const T &qiZ,
	                                                             T &aiX,       T &aiY,       T &aiZ,
	                                                             T &closNeighi,
	                                                       const T &mj,
	                                                       const T &qjX, const T &qjY, const T &qjZ);

private:
	void init();
};

#include "SimulationNBodyV3.hxx"

#endif /* SIMULATION_N_BODY_V3_H_ */
