/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef OGL_SPHERES_VISU_NO_H_
#define OGL_SPHERES_VISU_NO_H_

#include <map>
#include <string>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>

#include "OGLSpheresVisu.h"

template <typename T = double>
class OGLSpheresVisuNo :  public OGLSpheresVisu<T> {
public:
	OGLSpheresVisuNo();

	virtual ~OGLSpheresVisuNo();

	void refreshDisplay();
};

#include "OGLSpheresVisuNo.hxx"

#endif /* OGL_SPHERES_VISU_NO_H_ */
