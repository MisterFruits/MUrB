/*
 * SpheresVisuNo.h
 *
 *  Created on: 08 sept. 2014
 *      Author: Adrien Cassagne
 */

#ifndef OGL_SPHERES_VISU_NO_H_
#define OGL_SPHERES_VISU_NO_H_

#include "SpheresVisu.h"

template <typename T = double>
class SpheresVisuNo :  public SpheresVisu {
public:
	SpheresVisuNo();

	virtual ~SpheresVisuNo();

	void refreshDisplay();
	bool windowShouldClose();
};

#include "SpheresVisuNo.hxx"

#endif /* OGL_SPHERES_VISU_NO_H_ */
