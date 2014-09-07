/*
 * OGLSpheresVisuInst.h
 *
 *  Created on: 06 sept. 2014
 *      Author: Adrien Cassagne
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
	                   const vec_t<T> *positionsX,
	                   const vec_t<T> *positionsY,
	                   const vec_t<T> *positionsZ,
	                   const vec_t<T> *radius,
	                   const unsigned long nSpheres);

	virtual ~OGLSpheresVisuInst();
	void refreshDisplay();
};

#include "OGLSpheresVisuInst.hxx"

#endif /* OGL_SPHERES_VISU_INST_H_ */
