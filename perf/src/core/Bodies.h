/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef BODIES_H_
#define BODIES_H_

#include <string>

template <typename T>
class SimulationNBodyMPI;

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
	friend SimulationNBodyMPI<T>;

protected:
	unsigned long  n;
	T             *masses;
	T             *radiuses;
	vector3<T>     positions;
	vector3<T>     velocities;
	unsigned long  nVecs;
	unsigned short padding;
	// stats
	float          allocatedBytes;

public:
	Bodies();
	Bodies(const unsigned long n, const unsigned long randInit = 0);
	Bodies(const std::string inputFileName);
	Bodies(const Bodies<T>& bodies);

	virtual ~Bodies();

	Bodies<T>& operator=(const Bodies<T>& bodies);
	void hardCopy(const Bodies<T>& bodies);

	inline const unsigned long& getN() const;
	inline const unsigned long& getNVecs() const;
	inline const unsigned short& getPadding() const;
	inline const T* getMasses() const;
	inline const T* getRadiuses() const;
	inline const T* getPositionsX() const;
	inline const T* getPositionsY() const;
	inline const T* getPositionsZ() const;
	inline const T* getVelocitiesX() const;
	inline const T* getVelocitiesY() const;
	inline const T* getVelocitiesZ() const;
	inline const float& getAllocatedBytes() const;

	void updatePositionsAndVelocities(const vector3<T> &accelerations, T &dt);

	bool readFromFile(const std::string inputFileName);
	void write(std::ostream& stream, bool writeN = true) const;
	bool writeIntoFile(const std::string outputFileName) const;
	bool writeIntoFileMPI(const std::string outputFileName, const unsigned long MPINBodies = 0) const;

private:
	void deallocateBuffers();
	inline void setBody(const unsigned long &iBody,
	                    const T &mi, const T &radi,
	                    const T &qiX, const T &qiY, const T &qiZ,
	                    const T &viX, const T &viY, const T &viZ);
	bool read(std::istream& stream);
	void allocateBuffers();
	void initRandomly(const unsigned long randInit = 0);
	bool initFromFile(const std::string inputFileName);
};

template <typename T>
std::ostream& operator<<(std::ostream &o, const Bodies<T>& s);

#include "Bodies.hxx"

#endif /* BODIES_H_ */
