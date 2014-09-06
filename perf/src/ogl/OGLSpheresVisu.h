/*
 * OGLSpheresVisu.h
 *
 *  Created on: 30 août 2014
 *      Author: Adrien Cassagne
 */

#ifndef OGL_SPHERES_VISU_H_
#define OGL_SPHERES_VISU_H_

#include <map>
#include <string>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>

#include "OGLControl.h"

template <typename T = double>
class OGLSpheresVisu {
protected:
	GLFWwindow *window;

	const T *positionsX;
	float   *positionsXBuffer;
	const T *positionsY;
	float   *positionsYBuffer;
	const T *positionsZ;
	float   *positionsZBuffer;
	const T *radius;
	float   *radiusBuffer;

	const unsigned long nSpheres;

	GLuint vertexArrayRef;
	GLuint positionBufferRef[3];
	GLuint radiusBufferRef;
	GLuint mvpRef;
	GLuint shaderProgramRef;

	glm::mat4 mvp;

	OGLControl *control;

protected:
	OGLSpheresVisu(const std::string winName,
	               const int winWidth,
	               const int winHeight,
	               const T *positionsX,
	               const T *positionsY,
	               const T *positionsZ,
	               const T *radius,
	               const unsigned long nSpheres);

public:
	virtual ~OGLSpheresVisu();
	virtual void refreshDisplay() = 0;
	inline bool windowShouldClose();

protected:
	bool compileShaders(const std::vector<GLenum> shadersType, const std::vector<std::string> shadersFiles);
	void updatePositions();
};

#include "OGLSpheresVisu.hxx"

#endif /* OGL_SPHERES_VISU_H_ */
