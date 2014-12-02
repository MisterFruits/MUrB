/*!
 * \file    vertex150.glsl
 * \brief   This shader just gives vertices and radiuses to the geometry shader.
 * \author  A. Cassagne
 * \date    2014
 *
 * \section LICENSE
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode).
 *
 * \section DESCRIPTION
 * This is the traditional entry file for the code execution.
 */
#version 150
//#extension GL_ARB_explicit_attrib_location : enable

// Input vertex data, different for all executions of this shader.
in float positionXPerVertex;
in float positionYPerVertex;
in float positionZPerVertex;

// Input colors
in float radiusPerVertex;

// Output data ; will be interpolated for each fragment.
out float gRadius;

void main()
{
	// Output position of the vertex, in clip space : MVP * position
	//gl_Position = MVP * vec4(vertexPosition, 1);
	gl_Position = vec4(positionXPerVertex, positionYPerVertex, positionZPerVertex, 1);

	gRadius = radiusPerVertex;
}
