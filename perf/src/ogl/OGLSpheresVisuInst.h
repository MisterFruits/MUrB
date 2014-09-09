/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef OGL_SPHERES_VISU_INST_H_
#define OGL_SPHERES_VISU_INST_H_

#include <map>
#include <string>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>

#include "OGLSpheresVisu.h"

template <typename T = double>
class OGLSpheresVisuInst :  public OGLSpheresVisu<T> {
private:
	const float          PI               = 3.1415926;
	const unsigned long  nPointsPerCircle = 22;
	unsigned long        vertexModelSize;
	GLfloat             *vertexModel;
	GLuint               modelBufferRef;

public:
	OGLSpheresVisuInst(const std::string winName,
	                   const int winWidth,
	                   const int winHeight,
	                   const T *positionsX,
	                   const T *positionsY,
	                   const T *positionsZ,
	                   const T *radius,
	                   const unsigned long nSpheres);

	virtual ~OGLSpheresVisuInst();
	void refreshDisplay();
};

#include "OGLSpheresVisuInst.hxx"

#endif /* OGL_SPHERES_VISU_INST_H_ */
