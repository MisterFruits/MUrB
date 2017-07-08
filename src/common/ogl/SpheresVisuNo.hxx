/*!
 * \file    SpheresVisuNo.hxx
 * \brief   No visu.
 * \author  A. Cassagne
 * \date    2014
 *
 * \section LICENSE
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode).
 *
 * \section DESCRIPTION
 * This is the traditional entry file for the code execution.
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

template <typename T>
bool SpheresVisuNo<T>::pressedSpaceBar()
{
	return false;
}

template <typename T>
bool SpheresVisuNo<T>::pressedPageUp()
{
	return false;
}

template <typename T>
bool SpheresVisuNo<T>::pressedPageDown()
{
	return false;
}
