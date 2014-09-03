/*
 * OGLTools.h
 *
 *  Created on: 29 août 2014
 *      Author: adrien
 */

#ifndef OGLTOOLS_H_
#define OGLTOOLS_H_

#include <string>
#include <vector>

#include <GLFW/glfw3.h>

class OGLTools {
public:
	static GLFWwindow* initAndMakeWindow(const int winWidth, const int winHeight, const std::string winName);

	static GLuint loadShaderFromFile(const GLenum shaderType, const std::string shaderFilePath);

	static GLuint linkShaders(const std::vector<GLuint> shaders);
};

#endif /* OGLTOOLS_H_ */
