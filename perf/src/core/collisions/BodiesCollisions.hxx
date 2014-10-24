/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#include <cmath>
#include <limits>
#include <string>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sys/stat.h>

#include "../../utils/myIntrinsicsPlusPlus.h"

#include "BodiesCollisions.h"

template <typename T>
BodiesCollisions<T>::BodiesCollisions()
	: Bodies<T>()
{
}

template <typename T>
BodiesCollisions<T>::BodiesCollisions(const unsigned long n, const unsigned long randInit)
	: Bodies<T>(n, randInit)
{
}

template <typename T>
BodiesCollisions<T>::BodiesCollisions(const std::string inputFileName)
	: Bodies<T>(inputFileName)
{
}

template <typename T>
BodiesCollisions<T>::BodiesCollisions(const BodiesCollisions<T>& bodies)
	: Bodies<T>(bodies)
{
}

template <typename T>
BodiesCollisions<T>::~BodiesCollisions()
{
}

/* collisions 2D: work only if mi == mj (pool elastics collisions)
template <typename T>
void BodiesCollisions<T>::applyCollisions(std::vector<std::vector<unsigned long>> collisions)
{
	for(unsigned long iBody = 0; iBody < this->n; iBody++)
	{
		if(collisions[iBody].size())
		{
			std::cout << "Collisions for body n°" << iBody << ":" << std::endl;

			for(unsigned long jBody = 0; jBody < collisions[iBody].size(); jBody++)
			{
				std::cout << "  - body n°" << collisions[iBody][jBody] << std::endl;
			}
		}
	}


	for(unsigned long iBody = 0; iBody < this->n; iBody++)
		if(collisions[iBody].size())
		{
			assert(collisions[iBody].size() < 2);
			unsigned long jBody = collisions[iBody][0];

			// compute iBody ------------------------------------------------------------------------------------------
			T diffPosX = this->positions.x[jBody] - this->positions.x[iBody];
			T diffPosY = this->positions.y[jBody] - this->positions.y[iBody];

			T squareDist = (diffPosX * diffPosX) + (diffPosY * diffPosY);
			T dist = std::sqrt(squareDist);

			T normX = diffPosX / dist;
			T normY = diffPosY / dist;

			T perpendicularX = -normY;
			T perpendicularY =  normX;

			const T iVelocityX = this->velocities.x[iBody];
			const T iVelocityY = this->velocities.y[iBody];

			const T iSquareVelocityDist = (iVelocityX * iVelocityX) + (iVelocityY * iVelocityY);
			const T iVelocityDist = std::sqrt(iSquareVelocityDist);

			const T iCosTeta = ((perpendicularX * iVelocityX) + (perpendicularY * iVelocityY)) / iVelocityDist;
			const T iTeta = std::acos(iCosTeta);

			const T tmpIVelX = perpendicularX * (iCosTeta * iVelocityDist);
			const T tmpIVelY = perpendicularY * (iCosTeta * iVelocityDist);

			const T tmpJVelX = normX * (std::sin(iTeta) * iVelocityDist);
			const T tmpJVelY = normY * (std::sin(iTeta) * iVelocityDist);

			// compute jBody ------------------------------------------------------------------------------------------
			normX = -normX;
			normY = -normY;

			perpendicularX = -normY;
			perpendicularY =  normX;

			const T jVelocityX = this->velocities.x[jBody];
			const T jVelocityY = this->velocities.y[jBody];

			const T jSquareVelocityDist = (jVelocityX * jVelocityX) + (jVelocityY * jVelocityY);
			const T jVelocityDist = std::sqrt(jSquareVelocityDist);

			const T jCosTeta = ((perpendicularX * jVelocityX) + (perpendicularY * jVelocityY)) / jVelocityDist;
			const T jTeta = std::acos(jCosTeta);

			this->velocities.x[jBody] = perpendicularX * (jCosTeta * jVelocityDist) + tmpJVelX;
			this->velocities.y[jBody] = perpendicularY * (jCosTeta * jVelocityDist) + tmpJVelY;

			this->velocities.x[iBody] = normX * (std::sin(jTeta) * jVelocityDist) + tmpIVelX;
			this->velocities.y[iBody] = normY * (std::sin(jTeta) * jVelocityDist) + tmpIVelY;

			collisions[jBody].clear();
		}
}
*/

/* collisions 2D
template <typename T>
void BodiesCollisions<T>::applyCollisions(std::vector<std::vector<unsigned long>> collisions)
{
	for(unsigned long iBody = 0; iBody < this->n; iBody++)
	{
		if(collisions[iBody].size())
		{
			std::cout << "Collisions for body n°" << iBody << ":" << std::endl;

			for(unsigned long jBody = 0; jBody < collisions[iBody].size(); jBody++)
			{
				std::cout << "  - body n°" << collisions[iBody][jBody] << std::endl;
			}
		}
	}

	for(unsigned long iBody = 0; iBody < this->n; iBody++)
		if(collisions[iBody].size())
		{
			assert(collisions[iBody].size() < 2);
			unsigned long jBody = collisions[iBody][0];

			// 1. Find unit normal and unit tangent vectors
			T normX = this->positions.x[jBody] - this->positions.x[iBody];
			T normY = this->positions.y[jBody] - this->positions.y[iBody];

			T dij = std::sqrt((normX * normX) + (normY * normY));

			T uNormX = normX / dij;
			T uNormY = normY / dij;

			T uTangX = -uNormY;
			T uTangY =  uNormX;

			// 3 . Projecting the velocity vectors onto the unit normal and unit tangent vectors
			T viNorm = uNormX * this->velocities.x[iBody] + uNormY * this->velocities.y[iBody];
			T viTang = uTangX * this->velocities.x[iBody] + uTangY * this->velocities.y[iBody];

			T vjNorm = uNormX * this->velocities.x[jBody] + uNormY * this->velocities.y[jBody];
			T vjTang = uTangX * this->velocities.x[jBody] + uTangY * this->velocities.y[jBody];

			// 4. Find the new tangential velocities
			T viNewTang = viTang;
			T vjNewTang = vjTang;

			T iMass = this->masses[iBody];
			T jMass = this->masses[jBody];

			// 5. Find the new normal velocities
			T viNewNorm = (viNorm * (iMass - jMass) + (2.0 * jMass * vjNorm)) / (iMass + jMass);
			T vjNewNorm = (vjNorm * (jMass - iMass) + (2.0 * iMass * viNorm)) / (iMass + jMass);

			// 6. Convert the scalar normal and tangential velocities into vectors
			T viNewNormX = viNewNorm * uNormX;
			T viNewNormY = viNewNorm * uNormY;

			T viNewTangX = viNewTang * uTangX;
			T viNewTangY = viNewTang * uTangY;

			T vjNewNormX = vjNewNorm * uNormX;
			T vjNewNormY = vjNewNorm * uNormY;

			T vjNewTangX = vjNewTang * uTangX;
			T vjNewTangY = vjNewTang * uTangY;

			// 7. Find the final velocity vectors by adding the normal and tangential components
			this->velocities.x[iBody] = viNewNormX + viNewTangX;
			this->velocities.y[iBody] = viNewNormY + viNewTangY;

			this->velocities.x[jBody] = vjNewNormX + vjNewTangX;
			this->velocities.y[jBody] = vjNewNormY + vjNewTangY;

			collisions[jBody].clear();
		}
}
*/

/* collisions 3D */
template <typename T>
void BodiesCollisions<T>::applyCollisions(std::vector<std::vector<unsigned long>> collisions)
{
	/*
	for(unsigned long iBody = 0; iBody < this->n; iBody++)
		if(collisions[iBody].size())
		{
			std::cout << "Collisions for body n°" << iBody << ":" << std::endl;

			for(unsigned long jBody = 0; jBody < collisions[iBody].size(); jBody++)
				std::cout << "  - with body n°" << collisions[iBody][jBody] << std::endl;
		}
	*/

	for(unsigned long iBody = 0; iBody < this->n; iBody++)
		if(collisions[iBody].size())
		{
			assert(collisions[iBody].size() < 2);
			unsigned long jBody = collisions[iBody][0];

			// 1. Find unit normal and unit tangent 1 and tangent 2 vectors
			T normX = this->positions.x[jBody] - this->positions.x[iBody];
			T normY = this->positions.y[jBody] - this->positions.y[iBody];
			T normZ = this->positions.z[jBody] - this->positions.z[iBody];

			T dNorm = std::sqrt((normX * normX) + (normY * normY) + (normZ * normZ));

			T uNormX = normX / dNorm;
			T uNormY = normY / dNorm;
			T uNormZ = normZ / dNorm;
			//std::cout << "uNorm = {" << uNormX << ", " << uNormY << ", " << uNormZ << "}" << std::endl;

			// uTan1 and uTan2 vectors define the collision tangent plan (for 3D collisions)
			T tan1X, tan1Y, tan1Z;
			if(uNormX != 0)
			{
				tan1X = -(uNormY / uNormX);
				tan1Y = 1.0;
				tan1Z = 0.0;
			}
			else if(uNormY != 0)
			{
				tan1X = 1.0;
				tan1Y = -(uNormX / uNormY);
				tan1Z = 0.0;
			}
			else if(uNormZ != 0)
			{
				tan1X = 1.0;
				tan1Y = 0.0;
				tan1Z = -(uNormX / uNormZ);
			}
			else
			{
				std::cout << "Error: uNorm = {" << uNormX << ", " << uNormY << ", " << uNormZ << "}" << std::endl;
				exit(-1);
			}

			T dTan1 = std::sqrt((tan1X * tan1X) + (tan1Y * tan1Y) + (tan1Z * tan1Z));

			//std::cout << "tan1 = {" << tan1X << ", " << tan1Y << ", " << tan1Z << "}, "
			//          << "dTan1 = " << dTan1 << std::endl;

			T uTan1X = tan1X / dTan1;
			T uTan1Y = tan1Y / dTan1;
			T uTan1Z = tan1Z / dTan1;

			//std::cout << "uTan1 = {" << uTan1X << ", " << uTan1Y << ", " << uTan1Z << "}" << std::endl;
			//std::cout << "uNorm.uTan1 = " << ((uNormX * uTan1X) + (uNormY * uTan1Y) + (uNormZ * uTan1Z)) << std::endl;

			// uTan2 vector = (uNormX vector) ^ (uTan1 vector); (cross product or vector product)
			T uTan2X = (uNormY * uTan1Z) - (uNormZ * uTan1Y);
			T uTan2Y = (uNormZ * uTan1X) - (uNormX * uTan1Z);
			T uTan2Z = (uNormX * uTan1Y) - (uNormY * uTan1X);
			//std::cout << "uTan2 = {" << uTan2X << ", " << uTan2Y << ", " << uTan2Z << "}" << std::endl;

			//std::cout << "uNorm.uTan2 = " << ((uNormX * uTan2X) + (uNormY * uTan2Y) + (uNormZ * uTan2Z)) << std::endl;
			//std::cout << "uTan1.uTan2 = " << ((uTan1X * uTan2X) + (uTan1Y * uTan2Y) + (uTan1Z * uTan2Z)) << std::endl;

			// 2. Create the initial (before the collision) velocity vectors, iVel and jVel
			T iVelX = this->velocities.x[iBody];
			T iVelY = this->velocities.y[iBody];
			T iVelZ = this->velocities.z[iBody];

			T jVelX = this->velocities.x[jBody];
			T jVelY = this->velocities.y[jBody];
			T jVelZ = this->velocities.z[jBody];

			// 3. Projecting the velocity vectors onto the unit normal and unit tangent vectors
			// (scalar product or dot product)
			T viNorm = (uNormX * iVelX) + (uNormY * iVelY) + (uNormZ * iVelZ);
			T viTan1 = (uTan1X * iVelX) + (uTan1Y * iVelY) + (uTan1Z * iVelZ);
			T viTan2 = (uTan2X * iVelX) + (uTan2Y * iVelY) + (uTan2Z * iVelZ);

			T vjNorm = (uNormX * jVelX) + (uNormY * jVelY) + (uNormZ * jVelZ);
			T vjTan1 = (uTan1X * jVelX) + (uTan1Y * jVelY) + (uTan1Z * jVelZ);
			T vjTan2 = (uTan2X * jVelX) + (uTan2Y * jVelY) + (uTan2Z * jVelZ);

			// 4. Find the new tangential velocities
			T viNewTan1 = viTan1;
			T viNewTan2 = viTan2;

			T vjNewTan1 = vjTan1;
			T vjNewTan2 = vjTan2;

			T iMass = this->masses[iBody];
			T jMass = this->masses[jBody];

			// 5. Find the new normal velocities (apply elastic collision based on momentum and kinetic energy conservation)
			T viNewNorm = (viNorm * (iMass - jMass) + (2.0 * jMass * vjNorm)) / (iMass + jMass);
			T vjNewNorm = (vjNorm * (jMass - iMass) + (2.0 * iMass * viNorm)) / (iMass + jMass);

			// 6. Convert the scalar normal and tangential velocities into vectors
			T viNewNormX = viNewNorm * uNormX;
			T viNewNormY = viNewNorm * uNormY;
			T viNewNormZ = viNewNorm * uNormZ;

			T viNewTan1X = viNewTan1 * uTan1X;
			T viNewTan1Y = viNewTan1 * uTan1Y;
			T viNewTan1Z = viNewTan1 * uTan1Z;

			T viNewTan2X = viNewTan2 * uTan2X;
			T viNewTan2Y = viNewTan2 * uTan2Y;
			T viNewTan2Z = viNewTan2 * uTan2Z;

			T vjNewNormX = vjNewNorm * uNormX;
			T vjNewNormY = vjNewNorm * uNormY;
			T vjNewNormZ = vjNewNorm * uNormZ;

			T vjNewTan1X = vjNewTan1 * uTan1X;
			T vjNewTan1Y = vjNewTan1 * uTan1Y;
			T vjNewTan1Z = vjNewTan1 * uTan1Z;

			T vjNewTan2X = vjNewTan2 * uTan2X;
			T vjNewTan2Y = vjNewTan2 * uTan2Y;
			T vjNewTan2Z = vjNewTan2 * uTan2Z;

			// 7. Find the final velocity vectors by adding the normal and tangential components
			this->velocities.x[iBody] = viNewNormX + viNewTan1X + viNewTan2X;
			this->velocities.y[iBody] = viNewNormY + viNewTan1Y + viNewTan2Y;
			this->velocities.z[iBody] = viNewNormZ + viNewTan1Z + viNewTan2Z;

			this->velocities.x[jBody] = vjNewNormX + vjNewTan1X + vjNewTan2X;
			this->velocities.y[jBody] = vjNewNormY + vjNewTan1Y + vjNewTan2Y;
			this->velocities.z[jBody] = vjNewNormZ + vjNewTan1Z + vjNewTan2Z;

			collisions[jBody].clear();
		}
}

/* multi collisions 3D: not working */
template <typename T>
void BodiesCollisions<T>::applyMultiCollisions(std::vector<std::vector<unsigned long>> collisions)
{
	/*
	for(unsigned long iBody = 0; iBody < this->n; iBody++)
		if(collisions[iBody].size())
		{
			std::cout << "Collisions for body n°" << iBody << ":" << std::endl;

			for(unsigned long jBody = 0; jBody < collisions[iBody].size(); jBody++)
				std::cout << "  - with body n°" << collisions[iBody][jBody] << std::endl;
		}
	*/

	for(unsigned long iBody = 0; iBody < this->n; iBody++)
	{
		for(unsigned long iCollision = 0; iCollision < collisions[iBody].size(); iCollision++)
		{
			unsigned long jBody = collisions[iBody][iCollision];

			if(jBody > iBody)
			{
				// 1. Find unit normal and unit tangent 1 and tangent 2 vectors
				T normX = this->positions.x[jBody] - this->positions.x[iBody];
				T normY = this->positions.y[jBody] - this->positions.y[iBody];
				T normZ = this->positions.z[jBody] - this->positions.z[iBody];

				T dNorm = std::sqrt((normX * normX) + (normY * normY) + (normZ * normZ));

				T uNormX = normX / dNorm;
				T uNormY = normY / dNorm;
				T uNormZ = normZ / dNorm;

				// uTan1 and uTan2 vectors define the collision tangent plan (for 3D collisions)
				T tan1X, tan1Y, tan1Z;
				if(uNormX != 0)
				{
					tan1X = -(uNormY / uNormX);
					tan1Y = 1.0;
					tan1Z = 0.0;
				}
				else if(uNormY != 0)
				{
					tan1X = 1.0;
					tan1Y = -(uNormX / uNormY);
					tan1Z = 0.0;
				}
				else if(uNormZ != 0)
				{
					tan1X = 1.0;
					tan1Y = 0.0;
					tan1Z = -(uNormX / uNormZ);
				}
				else
				{
					std::cout << "Error: uNorm = {" << uNormX << ", " << uNormY << ", " << uNormZ << "}" << std::endl;
					exit(-1);
				}

				T dTan1 = std::sqrt((tan1X * tan1X) + (tan1Y * tan1Y) + (tan1Z * tan1Z));

				T uTan1X = tan1X / dTan1;
				T uTan1Y = tan1Y / dTan1;
				T uTan1Z = tan1Z / dTan1;

				// uTan2 vector = (uNormX vector) ^ (uTan1 vector); (cross product or vector product)
				T uTan2X = (uNormY * uTan1Z) - (uNormZ * uTan1Y);
				T uTan2Y = (uNormZ * uTan1X) - (uNormX * uTan1Z);
				T uTan2Z = (uNormX * uTan1Y) - (uNormY * uTan1X);

				// 2. Create the initial (before the collision) velocity vectors, iVel and jVel
				T iVelX = this->velocities.x[iBody];
				T iVelY = this->velocities.y[iBody];
				T iVelZ = this->velocities.z[iBody];

				T jVelX = this->velocities.x[jBody];
				T jVelY = this->velocities.y[jBody];
				T jVelZ = this->velocities.z[jBody];

				// 3. Projecting the velocity vectors onto the unit normal and unit tangent vectors
				// (scalar product or dot product)
				T viNorm = (uNormX * iVelX) + (uNormY * iVelY) + (uNormZ * iVelZ);
				T viTan1 = (uTan1X * iVelX) + (uTan1Y * iVelY) + (uTan1Z * iVelZ);
				T viTan2 = (uTan2X * iVelX) + (uTan2Y * iVelY) + (uTan2Z * iVelZ);

				T vjNorm = (uNormX * jVelX) + (uNormY * jVelY) + (uNormZ * jVelZ);
				T vjTan1 = (uTan1X * jVelX) + (uTan1Y * jVelY) + (uTan1Z * jVelZ);
				T vjTan2 = (uTan2X * jVelX) + (uTan2Y * jVelY) + (uTan2Z * jVelZ);

				// 4. Find the new tangential velocities
				T viNewTan1 = viTan1;
				T viNewTan2 = viTan2;

				T vjNewTan1 = vjTan1;
				T vjNewTan2 = vjTan2;

				T iMass = this->masses[iBody];
				T jMass = this->masses[jBody];

				// 5. Find the new normal velocities (apply elastic collision based on momentum and kinetic energy conservation)
				T viNewNorm = (viNorm * (iMass - jMass) + (2.0 * jMass * vjNorm)) / (iMass + jMass);
				T vjNewNorm = (vjNorm * (jMass - iMass) + (2.0 * iMass * viNorm)) / (iMass + jMass);

				// 6. Convert the scalar normal and tangential velocities into vectors
				T viNewNormX = viNewNorm * uNormX;
				T viNewNormY = viNewNorm * uNormY;
				T viNewNormZ = viNewNorm * uNormZ;

				T viNewTan1X = viNewTan1 * uTan1X;
				T viNewTan1Y = viNewTan1 * uTan1Y;
				T viNewTan1Z = viNewTan1 * uTan1Z;

				T viNewTan2X = viNewTan2 * uTan2X;
				T viNewTan2Y = viNewTan2 * uTan2Y;
				T viNewTan2Z = viNewTan2 * uTan2Z;

				T vjNewNormX = vjNewNorm * uNormX;
				T vjNewNormY = vjNewNorm * uNormY;
				T vjNewNormZ = vjNewNorm * uNormZ;

				T vjNewTan1X = vjNewTan1 * uTan1X;
				T vjNewTan1Y = vjNewTan1 * uTan1Y;
				T vjNewTan1Z = vjNewTan1 * uTan1Z;

				T vjNewTan2X = vjNewTan2 * uTan2X;
				T vjNewTan2Y = vjNewTan2 * uTan2Y;
				T vjNewTan2Z = vjNewTan2 * uTan2Z;

				// 7. Find the final velocity vectors by adding the normal and tangential components
				this->velocities.x[iBody] = viNewNormX + viNewTan1X + viNewTan2X;
				this->velocities.y[iBody] = viNewNormY + viNewTan1Y + viNewTan2Y;
				this->velocities.z[iBody] = viNewNormZ + viNewTan1Z + viNewTan2Z;

				this->velocities.x[jBody] = vjNewNormX + vjNewTan1X + vjNewTan2X;
				this->velocities.y[jBody] = vjNewNormY + vjNewTan1Y + vjNewTan2Y;
				this->velocities.z[jBody] = vjNewNormZ + vjNewTan1Z + vjNewTan2Z;
			}
		}
	}
}
