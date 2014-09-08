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
                                  const vec_t<T> *positionsX,
                                  const vec_t<T> *positionsY,
                                  const vec_t<T> *positionsZ,
                                  const vec_t<T> *radius,
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

	this->positionsXBuffer = new float[this->nSpheres];
	this->positionsYBuffer = new float[this->nSpheres];
	this->positionsZBuffer = new float[this->nSpheres];
	this->radiusBuffer     = new float[this->nSpheres];

	// copying the radius values
	for(unsigned long iVec = 0; iVec < ceil(this->nSpheres * 1.0 / VECTOR_SIZE); iVec++)
		for(unsigned short iSphere = 0; iSphere < VECTOR_SIZE; iSphere++)
		{
			unsigned long realSpherePos = iSphere + iVec * VECTOR_SIZE;
			if(realSpherePos < this->nSpheres)
				this->radiusBuffer[realSpherePos] = (float) this->radius[iVec].vec_data[iSphere];
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
OGLSpheresVisu<T>::OGLSpheresVisu()
	: window           (NULL),
	  positionsX       (NULL),
	  positionsXBuffer (NULL),
	  positionsY       (NULL),
	  positionsYBuffer (NULL),
	  positionsZ       (NULL),
	  positionsZBuffer (NULL),
	  radius           (NULL),
	  radiusBuffer     (NULL),
	  nSpheres         (0),
	  vertexArrayRef   ((GLuint) 0),
	  positionBufferRef{(GLuint) 0, (GLuint) 0, (GLuint) 0},
	  radiusBufferRef  ((GLuint) 0),
	  mvpRef           ((GLuint) 0),
	  shaderProgramRef ((GLuint) 0),
	  mvp              (glm::mat4(1.0f)),
	  control          (NULL)
{
}

template <typename T>
OGLSpheresVisu<T>::~OGLSpheresVisu()
{
	if(this->window)
		glfwDestroyWindow(this->window);

	if(this->control)
		delete this->control;

	if(this->positionsXBuffer)
		delete[] this->positionsXBuffer;
	if(this->positionsYBuffer)
		delete[] this->positionsYBuffer;
	if(this->positionsZBuffer)
		delete[] this->positionsZBuffer;
	if(this->radiusBuffer)
		delete[] this->radiusBuffer;
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
	// convert positions in float buffers
	for(unsigned long iVec = 0; iVec < ceil(this->nSpheres * 1.0 / VECTOR_SIZE); iVec++)
		for(unsigned short iSphere = 0; iSphere < VECTOR_SIZE; iSphere++)
		{
			unsigned long realSpherePos = iSphere + iVec * VECTOR_SIZE;
			if(realSpherePos < this->nSpheres)
			{
				this->positionsXBuffer[realSpherePos] = (float) this->positionsX[iVec].vec_data[iSphere];
				this->positionsYBuffer[realSpherePos] = (float) this->positionsY[iVec].vec_data[iSphere];
				this->positionsZBuffer[realSpherePos] = (float) this->positionsZ[iVec].vec_data[iSphere];
			}
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
