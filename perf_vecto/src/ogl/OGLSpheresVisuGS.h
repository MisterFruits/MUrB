/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
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
