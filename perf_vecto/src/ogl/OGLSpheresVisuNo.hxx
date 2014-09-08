/*
 * OGLSpheresVisuNo.cpp
 *
 *  Created on: 08 sept. 2014
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
#include "OGLSpheresVisuNo.h"

template <typename T>
OGLSpheresVisuNo<T>::OGLSpheresVisuNo()
	: OGLSpheresVisu<T>()
{
}

template <typename T>
OGLSpheresVisuNo<T>::~OGLSpheresVisuNo()
{
}

template <typename T>
void OGLSpheresVisuNo<T>::refreshDisplay()
{
}
