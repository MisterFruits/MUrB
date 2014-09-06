/*
 * OGLSpheresVisuGS.h
 *
 *  Created on: 06 sept. 2014
 *      Author: Adrien Cassagne
 */

#ifndef OGL_SPHERES_VISU_GS_H_
#define OGL_SPHERES_VISU_GS_H_

#include <map>
#include <string>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>

#include "OGLSpheresVisu.h"

template <typename T = double>
class OGLSpheresVisuGS :  public OGLSpheresVisu<T> {
public:
	OGLSpheresVisuGS(const std::string winName,
	                 const int winWidth,
	                 const int winHeight,
	                 const T *positionsX,
	                 const T *positionsY,
	                 const T *positionsZ,
	                 const T *radius,
	                 const unsigned long nSpheres);

	virtual ~OGLSpheresVisuGS();

	void refreshDisplay();
};

#include "OGLSpheresVisuGS.hxx"

#endif /* OGL_SPHERES_VISU_GS_H_ */
