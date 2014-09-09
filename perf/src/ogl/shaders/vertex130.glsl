/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */
 
#version 130
//#extension GL_ARB_explicit_attrib_location : enable

// Input vertex data, different for all executions of this shader.
in vec3	 modelPerVertex;
in float positionXPerVertex;
in float positionYPerVertex;
in float positionZPerVertex;

// Input colors
in float radiusPerVertex;

// Values that stay constant for the whole mesh.
uniform mat4 MVP;

void main()
{
	float scale = radiusPerVertex * 1.0e-7f;

	 // trick in order to avoid start point at the center of the sphere
	if(modelPerVertex.x == 0.0f && modelPerVertex.y == 0.0f && modelPerVertex.z == 0.0f)
		gl_Position = MVP * vec4(1.0e-8f * positionXPerVertex + (scale * 0),
		                         1.0e-8f * positionYPerVertex + (scale * 1),
		                         1.0e-8f * positionZPerVertex + (scale * 0),
								 1);
	else
		gl_Position = MVP * vec4(1.0e-8f * positionXPerVertex + (scale * modelPerVertex.x),
		                         1.0e-8f * positionYPerVertex + (scale * modelPerVertex.y),
		                         1.0e-8f * positionZPerVertex + (scale * modelPerVertex.z),
		                         1);
}
