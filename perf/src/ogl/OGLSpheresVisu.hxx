/*
 * OGLSpheresVisu.cpp
 *
 *  Created on: 30 août 2014
 *      Author: Adrien Cassagne
 */

#include <string>
#include <chrono>
#include <thread>
#include <vector>
#include <cassert>
#include <iostream>

#define GLM_FORCE_RADIANS
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>

#include "OGLTools.h"

#include "OGLSpheresVisu.h"

template <typename T>
OGLSpheresVisu<T>::OGLSpheresVisu(const string winName,
                                  const int winWidth,
                                  const int winHeight,
                                  const T *positionsX,
                                  const T *positionsY,
                                  const T *positionsZ,
                                  const T *radius,
                                  const unsigned long nSpheres)
	: window           (NULL),
	  positionsX       (positionsX),
	  positionsXBuffer (NULL),
	  positionsY       (positionsY),
	  positionsYBuffer (NULL),
	  positionsZ       (positionsZ),
	  positionsZBuffer (NULL),
	  radius           (radius),
	  radiusBuffer     (NULL),
	  nSpheres         (nSpheres),
	  vertexArrayRef   ((GLuint) 0),
	  positionBufferRef{(GLuint) 0, (GLuint) 0, (GLuint) 0},
	  radiusBufferRef  ((GLuint) 0),
	  mvpRef           ((GLuint) 0),
	  shaderProgramRef ((GLuint) 0),
	  mvp              (glm::mat4(1.0f)),
	  control          (NULL)
{
	assert(winWidth > 0);
	assert(winHeight > 0);
	assert(positionsX);
	assert(positionsY);
	assert(positionsZ);
	assert(radius);

	if(sizeof(T) == sizeof(float))
	{
#ifndef NBODY_DOUBLE //TODO: delete this define, this is just a patch to compile when using double
		this->positionsXBuffer = const_cast<T*>(positionsX); //TODO: do not use const_cast !
		this->positionsYBuffer = const_cast<T*>(positionsY); //TODO: do not use const_cast !
		this->positionsZBuffer = const_cast<T*>(positionsZ); //TODO: do not use const_cast !
		this->radiusBuffer     = const_cast<T*>(radius);     //TODO: do not use const_cast !
#endif
	}
	else
	{
		this->positionsXBuffer = new float[this->nSpheres];
		this->positionsYBuffer = new float[this->nSpheres];
		this->positionsZBuffer = new float[this->nSpheres];
		this->radiusBuffer     = new float[this->nSpheres];

		for(unsigned long iVertex = 0; iVertex < this->nSpheres; iVertex++)
			this->radiusBuffer[iVertex] = (float) this->radius[iVertex];
	}

	this->window = OGLTools::initAndMakeWindow(winWidth, winHeight, winName.c_str());

	if(this->window)
	{
		glGenVertexArrays(1, &(this->vertexArrayRef));
		glBindVertexArray(this->vertexArrayRef);

		// bind position and radius buffers to OpenGL system
		glGenBuffers(1, &this->positionBufferRef[0]); // can change over iterations, so binding is in refreshDisplay()
		glGenBuffers(1, &this->positionBufferRef[1]); // can change over iterations, so binding is in refreshDisplay()
		glGenBuffers(1, &this->positionBufferRef[2]); // can change over iterations, so binding is in refreshDisplay()

		glGenBuffers(1, &(this->radiusBufferRef));
		glBindBuffer(GL_ARRAY_BUFFER, this->radiusBufferRef);
		glBufferData(GL_ARRAY_BUFFER, this->nSpheres * sizeof(GLfloat), this->radiusBuffer, GL_STATIC_DRAW);

		// set background color to black
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

		// Create a control object in order to use mouse and keyboard (move in space)
		this->control = new OGLControl(this->window);

		// Enable depth test
		glEnable(GL_DEPTH_TEST);

		// Accept fragment if it closer to the camera than the former one
		glDepthFunc(GL_LESS);
	}
}

template <typename T>
OGLSpheresVisu<T>::~OGLSpheresVisu()
{
	if(this->window)
		glfwDestroyWindow(this->window);

	if(this->control)
		delete this->control;

	if(sizeof(T) != sizeof(float))
	{
		if(this->positionsXBuffer)
			delete[] this->positionsXBuffer;
		if(this->positionsYBuffer)
			delete[] this->positionsYBuffer;
		if(this->positionsZBuffer)
			delete[] this->positionsZBuffer;
		if(this->radiusBuffer)
			delete[] this->radiusBuffer;
	}
}

//TODO: use map instead of two vectors ;-)
template <typename T>
bool OGLSpheresVisu<T>::compileShaders(const std::vector<GLenum> shadersType,
                                       const std::vector<std::string> shadersFiles)
{
	assert(shadersType.size() == shadersFiles.size());

	bool isFine = true;
	std::vector<GLuint> shaders;

	// load and compile shader programs
	for(int iShader = 0; iShader < shadersType.size(); iShader++)
	{
		GLuint shader = OGLTools::loadShaderFromFile(shadersType[iShader], shadersFiles[iShader]);
		if(shader == 0)
			isFine = false;
		shaders.push_back(shader);
	}

	// link shader program
	if((unsigned) (this->shaderProgramRef = OGLTools::linkShaders(shaders)) == 0)
		isFine = false;

	// ProjectionMatrix * ViewMatrix * ModelMatrix => MVP pattern (Model = identity here)
	// Get a handle for our "MVP" uniform
	this->mvpRef = glGetUniformLocation(this->shaderProgramRef, "MVP");

	for(int iShader = 0; iShader < shaders.size(); iShader++)
		glDeleteShader(shaders[iShader]);

	return isFine;
}

template <typename T>
void OGLSpheresVisu<T>::updatePositions()
{
	// convert positions in float (if necessary)
	if(sizeof(T) != sizeof(float))
		for(unsigned long iVertex = 0; iVertex < this->nSpheres; iVertex++)
		{
			this->positionsXBuffer[iVertex] = (float) this->positionsX[iVertex];
			this->positionsYBuffer[iVertex] = (float) this->positionsY[iVertex];
			this->positionsZBuffer[iVertex] = (float) this->positionsZ[iVertex];
		}

	// bind position buffers to GPU
	glBindBuffer(GL_ARRAY_BUFFER, this->positionBufferRef[0]);
	glBufferData(GL_ARRAY_BUFFER, this->nSpheres * sizeof(GLfloat), this->positionsXBuffer, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, this->positionBufferRef[1]);
	glBufferData(GL_ARRAY_BUFFER, this->nSpheres * sizeof(GLfloat), this->positionsYBuffer, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, this->positionBufferRef[2]);
	glBufferData(GL_ARRAY_BUFFER, this->nSpheres * sizeof(GLfloat), this->positionsZBuffer, GL_STATIC_DRAW);
}

template <typename T>
bool OGLSpheresVisu<T>::windowShouldClose()
{
	if(this->window)
		return (bool) glfwWindowShouldClose(this->window);
	else
		return false;
}