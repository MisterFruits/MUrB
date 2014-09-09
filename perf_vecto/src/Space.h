/*
 * Do not remove.
 * MPI/OpenMP training courses
 * Adrien Cassagne, ASA - CINES, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SPACE_H_
#define SPACE_H_

#include <string>

#include "utils/myIntrinsics.h"
#include "utils/typeVector.h"

template <typename T = double>
struct vector3
{
	vec_t<T>* __restrict x;
	vec_t<T>* __restrict y;
	vec_t<T>* __restrict z;
};

template <typename T = double>
class Space
{
public://TODO: public attributes is not a good idea...
	const T G = 6.67384e-11;

	unsigned long nBodies;
	unsigned long nVecs;
	vec_t<T>*     masses;
	vec_t<T>*     radiuses;
	vector3<T>    positions;
	vector3<T>    speeds;
	vector3<T>    accelerations;
	vec_t<T>*     closestNeighborDist;
	T             dt;
	bool          dtConstant;

public:
	Space(const unsigned long nBodies);
	Space(const std::string inputFileName);

	virtual ~Space();

	inline unsigned long getNBodies();

	inline void setDtConstant(T dtVal);
	inline void setDtVariable();
	inline T getDt();

	void computeBodiesAcceleration();
	void iComputeBodiesAcceleration();
	void findTimeStep();
	void updateBodiesPositionAndSpeed();

	bool read(std::istream& stream);
	void write(std::ostream& stream);
	void writeIntoFile(const std::string outputFileName);

private:
	void allocateBuffers();

	void setBody(const unsigned long &iVec, const unsigned short &iBody,
	             const T &mass, const T &radius,
	             const T &posX, const T &posY, const T &posZ,
	             const T &speedX, const T &speedY, const T &speedZ);

	void initBodiesRandomly();
	void initBodiesFromFile(const std::string inputFileName);
	
	void computeAccelerationBetweenTwoVectorOfBodies(const T* __restrict iVecPosX,
	                                                 const T* __restrict iVecPosY,
	                                                 const T* __restrict iVecPosZ,
	                                                       T* __restrict iVecAccsX,
	                                                       T* __restrict iVecAccsY,
	                                                       T* __restrict iVecAccsZ,
	                                                       T* __restrict iClosNeiDist,
	                                                 const T* __restrict jVecMasses,
	                                                 const T* __restrict jVecPosX,
	                                                 const T* __restrict jVecPosY,
	                                                 const T* __restrict jVecPosZ);

	inline void computeAccelerationBetweenTwoBodies(const T &iPosX, const T &iPosY, const T &iPosZ,
	                                                      T &iAccsX,      T &iAccsY,      T &iAccsZ,
	                                                      T &iClosNeiDist,
	                                                const T &jMasses,
	                                                const T &jPosX, const T &jPosY, const T &jPosZ);

	inline void iComputeAccelerationBetweenTwoVectorOfBodies(const T* __restrict iVecPosX,
	                                                         const T* __restrict iVecPosY,
	                                                         const T* __restrict iVecPosZ,
	                                                               T* __restrict iVecAccsX,
	                                                               T* __restrict iVecAccsY,
	                                                               T* __restrict iVecAccsZ,
	                                                               T* __restrict iClosNeiDist,
	                                                         const T* __restrict jVecMasses,
	                                                         const T* __restrict jVecPosX,
	                                                         const T* __restrict jVecPosY,
	                                                         const T* __restrict jVecPosZ);

	static inline void iComputeAccelerationBetweenTwoBodies(const vec &rG,
	                                                        const vec &rIPosX,
	                                                        const vec &rIPosY,
	                                                        const vec &rIPosZ,
	                                                              vec &rIAccX,
	                                                              vec &rIAccY,
	                                                              vec &rIAccZ,
	                                                              vec &rIClosNeiDist,
	                                                              vec &rJMass,
	                                                              vec &rJPosX,
	                                                              vec &rJPosY,
	                                                              vec &rJPosZ);

	inline T computeTimeStep(const unsigned long iVec);
};

template <typename T>
std::ostream& operator<<(std::ostream &o, const Space<T>& s);

#include "Space.hxx"

#endif /* SPACE_H_ */
