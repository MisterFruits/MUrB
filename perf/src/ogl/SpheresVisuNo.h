/*!
 * \file    SpheresVisuNo.h
 * \brief   No visu.
 * \author  A. Cassagne
 * \date    2014
 *
 * \section LICENSE
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode).
 *
 * \section DESCRIPTION
 * This is the traditional entry file for the code execution.
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
