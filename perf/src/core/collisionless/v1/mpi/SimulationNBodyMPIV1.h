/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SIMULATION_N_BODY_MPI_V1_H_
#define SIMULATION_N_BODY_MPI_V1_H_

#include <string>

#include "../../SimulationNBodyMPI.h"

template <typename T = double>
class SimulationNBodyMPIV1 : public SimulationNBodyMPI<T>
{
public:
	SimulationNBodyMPIV1(const unsigned long nBodies);
	SimulationNBodyMPIV1(const std::string inputFileName);
	virtual ~SimulationNBodyMPIV1();

protected:
	virtual void initIteration();
	virtual void computeLocalBodiesAcceleration();
	virtual void computeNeighborBodiesAcceleration();

	static inline void computeAccelerationBetweenTwoBodies(const T &G,
	                                                       const T &qiX, const T &qiY, const T &qiZ,
	                                                             T &aiX,       T &aiY,       T &aiZ,
	                                                             T &closNeighi,
	                                                       const T &mj,
	                                                       const T &qjX, const T &qjY, const T &qjZ);

private:
	void init();
};

#include "SimulationNBodyMPIV1.hxx"

#endif /* SIMULATION_N_BODY_MPI_V1_H_ */
