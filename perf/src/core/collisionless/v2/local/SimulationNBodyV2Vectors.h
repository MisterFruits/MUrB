/*
 * Do not remove.
 * Gabriel Hautreux, CINES, gabrielhautreux@gmail.com 
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SIMULATION_N_BODY_V2_VECTORS_H_
#define SIMULATION_N_BODY_V2_VECTORS_H_

#include <string>

#include "../../SimulationNBodyLocal.h"

template <typename T = double>
class SimulationNBodyV2Vectors : public SimulationNBodyLocal<T>
{
public:
	SimulationNBodyV2Vectors(const unsigned long nBodies);
	SimulationNBodyV2Vectors(const std::string inputFileName);
	virtual ~SimulationNBodyV2Vectors();

protected:
	virtual void initIteration();
	virtual void computeLocalBodiesAcceleration();

	inline void computeAccelerationBetweenTwoBodies(const T &iPosX, const T &iPosY, const T &iPosZ,
	                                                      T &iAccsX,      T &iAccsY,      T &iAccsZ,
	                                                      T &iClosNeiDist,          const T &iMasses,
	                                                const T &jPosX, const T &jPosY, const T &jPosZ,
	                                                      T &jAccsX,      T &jAccsY,      T &jAccsZ,
	                                                      T &jClosNeiDist,          const T &jMasses);

private:
	void reAllocateBuffers();
};

#include "SimulationNBodyV2Vectors.hxx"

#endif /* SIMULATION_N_BODY_V2_VECTORS_H_ */
