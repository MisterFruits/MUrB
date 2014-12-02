/*!
 * \file    fragment130.glsl
 * \brief   Very simple fragment shader.
 * \author  A. Cassagne
 * \date    2014
 *
 * \section LICENSE
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode).
 *
 * \section DESCRIPTION
 * This is the traditional entry file for the code execution.
 */
#version 130
//#extension GL_ARB_explicit_attrib_location : enable

// Ouput data
out vec3 outColor;

void main()
{
	vec3 white = vec3(1, 1, 1);
	outColor = white;
}
