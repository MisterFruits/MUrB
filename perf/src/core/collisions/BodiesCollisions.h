/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef BODIES_COLLISIONS_H_
#define BODIES_COLLISIONS_H_

#include <string>

#include "../Bodies.h"

template <typename T = double>
class BodiesCollisions : public Bodies<T>
{
public:
	BodiesCollisions();
	BodiesCollisions(const unsigned long n, const unsigned long randInit = 0);
	BodiesCollisions(const std::string inputFileName);
	BodiesCollisions(const BodiesCollisions<T>& bodies);

	virtual ~BodiesCollisions();

	void applyCollisions(std::vector<std::vector<unsigned long>> collisions);
	void applyMultiCollisions(std::vector<std::vector<unsigned long>> collisions);
};

#include "BodiesCollisions.hxx"

#endif /* BODIES_COLLISIONS_H_ */
