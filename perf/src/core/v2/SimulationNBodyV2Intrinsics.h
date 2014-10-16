/*
 * Do not remove.
 * Gabriel Hautreux, CINES, gabrielhautreux@gmail.com 
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SIMULATION_N_BODY_V2_INTRINSICS
#define SIMULATION_N_BODY_V2_INTRINSICS

#include <string>

#include "../../utils/myIntrinsicsPlusPlus.h"

#include "../SimulationNBodyLocal.h"

template <typename T = double>
class SimulationNBodyV2Intrinsics : public SimulationNBodyLocal<T>
{
public:
	SimulationNBodyV2Intrinsics(const unsigned long nBodies);
	SimulationNBodyV2Intrinsics(const std::string inputFileName);
	virtual ~SimulationNBodyV2Intrinsics();

protected:
	virtual void initIteration();
	virtual void computeLocalBodiesAcceleration();

	inline void computeAccelerationBetweenTwoBodiesSelf(const T &iPosX, const T &iPosY, const T &iPosZ,
	                                                          T &iAccsX,      T &iAccsY,      T &iAccsZ,
	                                                          T &iClosNeiDist,          const T &iMasses,
	                                                    const T &jPosX, const T &jPosY, const T &jPosZ,
	                                                          T &jAccsX,      T &jAccsY,      T &jAccsZ,
	                                                          T &jClosNeiDist,          const T &jMasses);

	static inline void computeAccelerationBetweenTwoBodies(const mipp::vec &rG,
	                                                       const mipp::vec &rIPosX,
	                                                       const mipp::vec &rIPosY,
	                                                       const mipp::vec &rIPosZ,
	                                                             mipp::vec &rIAccX,
	                                                             mipp::vec &rIAccY,
	                                                             mipp::vec &rIAccZ,
	                                                             mipp::vec &rIClosNeiDist,
	                                                       const mipp::vec &rIMass,
	                                                       const mipp::vec &rJPosX,
	                                                       const mipp::vec &rJPosY,
	                                                       const mipp::vec &rJPosZ,
	                                                             mipp::vec &rJAccX,
	                                                             mipp::vec &rJAccY,
	                                                             mipp::vec &rJAccZ,
	                                                             mipp::vec &rJClosNeiDist,
	                                                       const mipp::vec &rJMass,
	                                                       const bool  dtConstant);
private:
	void _computeLocalBodiesAcceleration();
	void _reAllocateBuffers();
	void reAllocateBuffers();
	void _initIteration();
};

#include "SimulationNBodyV2Intrinsics.hxx"

#endif /* SIMULATION_N_BODY_V2_INTRINSICS */
