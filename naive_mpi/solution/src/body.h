/*
 * Do not remove.
 * MPI/OpenMP training courses
 * Adrien Cassagne, ASA - CINES, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef BODY_H_
#define BODY_H_

typedef struct
{
	double mass; // body mass
	double posX; // body position following x axis
	double posY; // body position following y axis
} body;

typedef struct
{
	body   *b;                 // body mass and body position
	double speedX;             // body speed following x axis
	double speedY;             // body speed following y axis
	double accelerationX;      // body acceleration following x axis
	double accelerationY;      // body acceleration following y axis
	double closestNeighborLen; // contains the distance with the closest neighbor
} localBody;

/* Initialize body b with mass and position */
void initBody(body* b, double mass, double posX, double posY);

/* Initialize local body lb with mass, position and speed */
void initLocalBody(localBody* lb, double mass, double posX, double posY, double speedX, double speedY);

/* Compute acceleration of localBody lb with body b, update closestNeighborLen if necessary */
void computeAcceleration (localBody *lb, const body *b);

/* Compute dt for localBody lb, using of closestNeighborLen */
double computeDt (const localBody lb);

/* Update new position and speed of localBody lb with given dt and his initial position, speed and acceleration */
void updatePositionAndSpeed(localBody *lb, const double dt);

#endif /* BODY_H_ */
