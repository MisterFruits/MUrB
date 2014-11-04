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
protected:
	const T G = 6.67384e-11;

	Bodies<T>  *bodies;
	vector3<T>  accelerations;
	T          *closestNeighborDist;
	bool        dtConstant;
	T           dt;
	T           minDt;
	// stats
	float       flopsPerIte;
	float       allocatedBytes;
	unsigned    nMaxThreads;

protected:
	SimulationNBody(Bodies<T> *bodies);

public:
	virtual ~SimulationNBody();

	inline const Bodies<T>* getBodies() const;
	inline void setDtConstant(T dtVal);
	inline void setDtVariable(T minDt);
	inline const T& getDt() const;
	inline const float& getFlopsPerIte() const;
	inline const float& getAllocatedBytes() const;

	virtual void computeOneIteration() = 0;

protected:
	virtual void initIteration() = 0;
	virtual void findTimeStep()  = 0;
	inline T computeTimeStep(const unsigned long iBody);

private:
	void allocateBuffers();
};

#include "SimulationNBody.hxx"

#endif /* SIMULATION_N_BODY_H_ */
