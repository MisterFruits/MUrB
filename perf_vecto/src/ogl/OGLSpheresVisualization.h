/*
 * OGLSpheresVisualization.h
 *
 *  Created on: 30 août 2014
 *      Author: Adrien Cassagne
 */

#ifndef OGL_SPHERES_VISUALIZATION_H_
#define OGL_SPHERES_VISUALIZATION_H_

#include <map>
#include <string>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>

#include "OGLControl.h"

#include "../utils/typeVector.h"

template <typename T = double>
class OGLSpheresVisualization {
private:
	GLFWwindow *window;

	const vec_t<T> *positionsX;
	float          *positionsXBuffer;
	const vec_t<T> *positionsY;
	float          *positionsYBuffer;
	const vec_t<T> *positionsZ;
	float          *positionsZBuffer;
	const vec_t<T> *radius;
	float          *radiusBuffer;

	const unsigned long nSpheres;

	GLuint vertexArrayRef;
	GLuint positionBufferRef[3];
	GLuint radiusBufferRef;
	GLuint mvpRef;
	GLuint shaderProgramRef;

	glm::mat4 mvp;

	OGLControl *control;

public:
	OGLSpheresVisualization(const std::string winName,
	                        const int winWidth,
	                        const int winHeight,
	                        const vec_t<T> *positionsX,
	                        const vec_t<T> *positionsY,
	                        const vec_t<T> *positionsZ,
	                        const vec_t<T> *radius,
	                        const unsigned long nSpheres);

	virtual ~OGLSpheresVisualization();

	void refreshDisplay();

	inline bool windowShouldClose();

private:
	bool compileShaders(const std::vector<GLenum> shadersType, const std::vector<std::string> shadersFiles);

};

#include "OGLSpheresVisualization.hxx"

#endif /* OGL_SPHERES_VISUALIZATION_H_ */
