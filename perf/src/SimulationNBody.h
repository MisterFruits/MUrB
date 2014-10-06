/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SIMULATION_N_BODY_H_
#define SIMULATION_N_BODY_H_

#include <string>

#include "Bodies.h"

template <typename T = double>
class SimulationNBody
{
public:
	const T G = 6.67384e-11;

	Bodies<T>     bodies;
	vector3<T>    accelerations;
	T            *closestNeighborDist;
	T             dt;
	bool          dtConstant;

public:
	SimulationNBody(const unsigned long nBodies);
	SimulationNBody(const std::string inputFileName);

private:
	void allocateBuffers();

public:
	virtual ~SimulationNBody();

	inline void setDtConstant(T dtVal);
	inline void setDtVariable();
	inline T getDt();

	void computeIteration();

private:
	void computeBodiesAcceleration();
	void computeBodiesAccelerationCB();
	inline void computeAccelerationBetweenTwoBodiesNaive(const unsigned long iBody, const unsigned long jBody);
	inline void computeAccelerationBetweenTwoBodies(const unsigned long iBody, const unsigned long jBody);

	void computeBodiesAccelerationV2();
	inline void computeAccelerationBetweenTwoBodiesNaiveV2(const unsigned long iBody, const unsigned long jBody);
	inline void computeAccelerationBetweenTwoBodiesV2(const unsigned long iBody, const unsigned long jBody);

	void findTimeStep();
	inline T computeTimeStep(const unsigned long iBody);
};

#include "SimulationNBody.hxx"

#endif /* SIMULATION_N_BODY_H_ */
