/*
 * SpheresVisuNo.cpp
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

#include "SpheresVisuNo.h"

template <typename T>
SpheresVisuNo<T>::SpheresVisuNo()
	: SpheresVisu()
{
}

template <typename T>
SpheresVisuNo<T>::~SpheresVisuNo()
{
}

template <typename T>
void SpheresVisuNo<T>::refreshDisplay()
{
}

template <typename T>
bool SpheresVisuNo<T>::windowShouldClose()
{
	return false;
}
