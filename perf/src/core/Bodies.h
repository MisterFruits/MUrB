/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef BODIES_H_
#define BODIES_H_

#include <string>

template <typename T = double>
struct vector3
{
	T *x;
	T *y;
	T *z;
};

template <typename T = double>
class Bodies
{
private:
	unsigned long n;
	T            *masses;
	T            *radiuses;
	vector3<T>    positions;
	vector3<T>    velocities;
	unsigned long nVecs;
	// stats
	float         allocatedBytes;

public:
	Bodies(const unsigned long n);
	Bodies(const std::string inputFileName);
	Bodies(const Bodies<T>& bodies);

	virtual ~Bodies();

	Bodies<T>& operator=(const Bodies<T>& bodies);

	inline const unsigned long& getN();
	inline const unsigned long& getNVecs();
	inline const T* getMasses();
	inline const T* getRadiuses();
	inline const T* getPositionsX();
	inline const T* getPositionsY();
	inline const T* getPositionsZ();
	inline const T* getVelocitiesX();
	inline const T* getVelocitiesY();
	inline const T* getVelocitiesZ();
	inline const float& getAllocatedBytes();

	void updatePositionsAndVelocities(const vector3<T> &accelerations, T &dt);

	void write(std::ostream& stream);
	void writeIntoFile(const std::string outputFileName);

private:
	inline void setBody(const unsigned long &iBody,
	                    const T &mass, const T &radius,
	                    const T &posX, const T &posY, const T &posZ,
	                    const T &velocityX, const T &velocityY, const T &velocityZ);
	bool read(std::istream& stream);
	void allocateBuffers();
	void initRandomly();
	void initFromFile(const std::string inputFileName);
};

template <typename T>
std::ostream& operator<<(std::ostream &o, const Bodies<T>& s);

#include "Bodies.hxx"

#endif /* BODIES_H_ */