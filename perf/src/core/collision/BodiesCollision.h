/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef BODIES_COLLISION_H_
#define BODIES_COLLISION_H_

#include <string>

#include "../Bodies.h"

template <typename T = double>
class BodiesCollision : public Bodies<T>
{
public:
	BodiesCollision();
	BodiesCollision(const unsigned long n, const unsigned long randInit = 0);
	BodiesCollision(const std::string inputFileName);
	BodiesCollision(const BodiesCollision<T>& bodies);

	virtual ~BodiesCollision();

	void applyCollisions(std::vector<std::vector<unsigned long>> collisions);
};

#include "BodiesCollision.hxx"

#endif /* BODIES_COLLISION_H_ */
