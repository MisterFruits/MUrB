/*!
 * \file    SimulationNBodyV2FineTuned.h
 * \brief   Implementation of SimulationNBodyLocal with intrinsic function calls, fine tuned (n²/2 computations).
 * \author  G. Hautreux
 * \date    2014
 *
 * \section LICENSE
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode).
 */
#ifndef SIMULATION_N_BODY_V2_FINETUNED
#define SIMULATION_N_BODY_V2_FINETUNED

#include <string>

#include "../../../../utils/myIntrinsicsPlusPlus.h"

#include "../../SimulationNBodyLocal.h"

/*!
 * \class  SimulationNBodyV2FineTuned
 * \brief  Implementation of SimulationNBodyLocal with intrinsic function calls, fine tuned (n²/2 computations).
 *
 * \tparam T : Type.
 */
template <typename T = double>
class SimulationNBodyV2FineTuned : public SimulationNBodyLocal<T>
{
public:
	SimulationNBodyV2FineTuned(const unsigned long nBodies);
	SimulationNBodyV2FineTuned(const std::string inputFileName);
	virtual ~SimulationNBodyV2FineTuned();

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

#include "SimulationNBodyV2FineTuned.hxx"

#endif /* SIMULATION_N_BODY_V2_FINETUNED */
