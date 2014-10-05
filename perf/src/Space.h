/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef SPACE_H_
#define SPACE_H_

#include <string>

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
public://TODO: public attributes is not a good idea...
	const T G = 6.67384e-11;

	unsigned long nBodies;
	T            *masses;
	T            *radiuses;
	vector3<T>    positions;
	vector3<T>    speeds;
	vector3<T>    accelerations;
	T            *closestNeighborDist;
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
	void computeBodiesAccelerationV2();
	void computeBodiesAccelerationCB();
	void findTimeStep();
	void updateBodiesPositionAndSpeed();

	bool read(std::istream& stream);
	void write(std::ostream& stream);
	void writeIntoFile(const std::string outputFileName);

private:
	void allocateBuffers();

	void setBody(const unsigned long &iBody,
	             const T &mass, const T &radius,
	             const T &posX, const T &posY, const T &posZ,
	             const T &speedX, const T &speedY, const T &speedZ);

	void initBodiesRandomly();
	void initBodiesFromFile(const std::string inputFileName);
	
	inline void computeAccelerationBetweenTwoBodies(const unsigned long iBody, const unsigned long jBody);
	void computeAccelerationBetweenTwoBodiesNaive(const unsigned long iBody, const unsigned long jBody);

	void computeAccelerationBetweenTwoBodiesNaiveV2(const unsigned long iBody, const unsigned long jBody);

	inline T computeTimeStep(const unsigned long iBody);
};

template <typename T>
std::ostream& operator<<(std::ostream &o, const Space<T>& s);

#include "Space.hxx"

#endif /* SPACE_H_ */
