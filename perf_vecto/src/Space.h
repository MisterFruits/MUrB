/*
 * Do not remove.
 * MPI/OpenMP training courses
 * Adrien Cassagne, ASA - CINES, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SPACE_H_
#define SPACE_H_

#include <string>

#include "utils/myIntrinsics.h" // needed for intrinsic prototypes describe below (see EXPERIMENTAL)
#include "utils/typeVector.h"

template <typename T = double>
struct vector3
{
	vec_t<T> *x;
	vec_t<T> *y;
	vec_t<T> *z;
};

template <typename T = double>
class Space
{
public://TODO: public attributes is not a good idea...
	const T G = 6.67384e-11;

	unsigned long nBodies;
	unsigned long nVecs;
	vec_t<T>     *masses;
	vec_t<T>     *radiuses;
	vector3<T>    positions;
	vector3<T>    speeds;
	vector3<T>    accelerations;
	vec_t<T>     *closestNeighborDist;
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
	void findTimeStep();
	void updateBodiesPositionAndSpeed();

	bool read(std::istream& stream);
	void write(std::ostream& stream);
	void writeIntoFile(const std::string outputFileName);

private:
	void allocateBuffers();

	void initBodiesRandomly();
	void initBodiesFromFile(const std::string inputFileName);
	
	inline void computeAccelerationBetweenTwoBodies(const unsigned long  iVec,
	                                                const unsigned long  jVec,
	                                                const unsigned short iBody,
	                                                const unsigned short jBody);
	inline void computeAccelerationBetweenTwoBodiesNaive(const unsigned long  iVec,
	                                                     const unsigned long  jVec,
	                                                     const unsigned short iBody,
	                                                     const unsigned short jBody);

	inline T computeTimeStep(const unsigned long iVec);

	// EXPERIMENTAL ===================================================================================================
	/* TODO: this part is commented because the code does not compile if we don't have AVX2 instructions
public:
	void vectorComputeBodiesAcceleration();
	void intrinComputeBodiesAcceleration();

private:
	inline void vectorComputeAccelerationBetweenBodies(const unsigned long iBody, const unsigned long jBody, const int vecDim);
	inline void selfVectorComputeAccelerationBetweenBodies(const unsigned long iBody, const int vecDim);

	inline void intrinComputeAccelerationBetweenBodies(const unsigned long iBody, const unsigned long jBody, const int vecDim, 
								vec px, vec py, vec pz, vec *accx, vec *accy, vec *accz, vec *closest);
	inline void selfIntrinComputeAccelerationBetweenBodies(const unsigned long iBody, const int vecDim,
								vec px, vec py, vec pz, vec *accx, vec *accy, vec *accz, vec *closest);
	*/
	// EXPERIMENTAL ===================================================================================================
};

template <typename T>
std::ostream& operator<<(std::ostream &o, const Space<T>& s);

#include "Space.hxx"

#endif /* SPACE_H_ */
