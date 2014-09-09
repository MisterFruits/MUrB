/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
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
