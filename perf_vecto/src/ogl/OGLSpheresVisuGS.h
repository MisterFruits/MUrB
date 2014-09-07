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
	                 const vec_t<T> *positionsX,
	                 const vec_t<T> *positionsY,
	                 const vec_t<T> *positionsZ,
	                 const vec_t<T> *radius,
	                 const unsigned long nSpheres);

	virtual ~OGLSpheresVisuGS();

	void refreshDisplay();
};

#include "OGLSpheresVisuGS.hxx"

#endif /* OGL_SPHERES_VISU_GS_H_ */
