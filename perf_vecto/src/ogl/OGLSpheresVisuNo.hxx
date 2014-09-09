/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
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
