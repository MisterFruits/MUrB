/*
 * Do not remove.
 * MPI/OpenMP training courses
 * Adrien Cassagne, ASA - CINES, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SPACE_H_
#define SPACE_H_

#include <string>
#include "./utils/myIntrinsics.h"

template <typename T = double>
struct vector3
{
	T *x;
	T *y;
	T *z;
};

template <typename T = double>
class Space
{
private:
	const T G = 6.67384e-11;

	unsigned long nBodies;
	T            *masses;
	vector3<T>    positions;
	vector3<T>    speeds;
	vector3<T>    accelerations;
	T            *closestNeighborLen;
	T             dt;

public:
	Space(const unsigned long nBodies);
	Space(const std::string inputFileName);

	virtual ~Space();

	inline unsigned long getNBodies();

	void computeBodiesAcceleration();
	void vectorComputeBodiesAcceleration();
	void intrinComputeBodiesAcceleration();
	void findTimeStep();
	void updateBodiesPositionAndSpeed();

	bool read(std::istream& stream);
	void write(std::ostream& stream);
	void writeIntoFile(const std::string outputFileName);

private:
	void allocateBuffers();

	void initBodiesRandomly();
	void initBodiesWithFile(const std::string inputFileName);

	
	inline void computeAccelerationBetweenTwoBodies(const unsigned long iBody, const unsigned long jBody);

	inline void vectorComputeAccelerationBetweenBodies(const unsigned long iBody, const unsigned long jBody, const int vecDim);
	inline void selfVectorComputeAccelerationBetweenBodies(const unsigned long iBody, const int vecDim);

	inline void intrinComputeAccelerationBetweenBodies(const unsigned long iBody, const unsigned long jBody, const int vecDim, 
								vec px, vec py, vec pz, vec *accx, vec *accy, vec *accz, vec *closest);
	inline void selfIntrinComputeAccelerationBetweenBodies(const unsigned long iBody, const int vecDim,
								vec px, vec py, vec pz, vec *accx, vec *accy, vec *accz, vec *closest);

	inline T computeTimeStep(const unsigned long iBody);
};

template <typename T>
std::ostream& operator<<(std::ostream &o, const Space<T>& s);

#include "Space.hxx"

#endif /* SPACE_H_ */
