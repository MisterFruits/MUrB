/*
 * Do not remove.
 * MPI/OpenMP training courses
 * Adrien Cassagne, ASA - CINES, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifdef NBODY_FLOAT
#define TYPE float
#else
#define TYPE double
#endif

#include <map>
#include <cmath>
#include <string>
#include <chrono>
#include <thread>
#include <vector>
#include <cassert>
#include <iostream>
using namespace std;

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
using namespace glm;

#include "ogl/OGLTools.h"
#include "ogl/OGLControl.h"
#include "ogl/OGLSpheresVisualization.h"

#include "utils/Perf.h"
#include "utils/ArgumentsReader.h"

#include "Space.h"

string        InputFileName;
string        OutputBaseName;
unsigned long NBodies;
unsigned long NIterations;
bool          Verbose = false;

int WinWidth  = 1280;
int WinHeight = 720;

/*
 * read args from command line and set global variables
 * usage: ./nbody -f fileName -i nIterations [-v] [-w]
 * */
bool argsReader1(int argc, char** argv)
{
	map<string, string> reqArgs, faculArgs, docArgs;
	ArgumentsReader argsReader(argc, argv);

	reqArgs  ["f"]     = "inputFileName";
	docArgs  ["f"]     = "bodies file name to read (you can put 'data/in/np1/in.testcase1.0.dat').";
	reqArgs  ["i"]     = "nIterations";
	docArgs  ["i"]     = "number of iterations to compute.";
	faculArgs["n"]     = "nBodies";
	docArgs  ["n"]     = "approximate number of bodies randomly generated (do not use with -f).";
	faculArgs["v"]     = "";
	docArgs  ["v"]     = "activate verbose mode.";
	faculArgs["w"]     = "outputFileName";
	docArgs  ["w"]     = "base name of body file(s) to write (you can put 'data/out/out').";
	faculArgs["h"]     = "";
	docArgs  ["h"]     = "display this help.";
	faculArgs["-help"] = "";
	docArgs  ["-help"] = "display this help.";

	if(argsReader.parseArguments(reqArgs, faculArgs))
	{
		if(argsReader.existArgument("h") || argsReader.existArgument("-help"))
		{
			if(argsReader.parseDocArgs(docArgs))
				argsReader.printUsage();
			else
				cout << "A problem was encountered when parsing arguments documentation... exiting." << endl;
			exit(-1);
		}

		InputFileName = argsReader.getArgument("f");
		NIterations   = stoi(argsReader.getArgument("i"));

		if(argsReader.existArgument("v"))
			Verbose = true;
		if(argsReader.existArgument("w"))
			OutputBaseName = argsReader.getArgument("w");
	}
	else
	{
		if(argsReader.parseDocArgs(docArgs))
			return false;
		else
			cout << "A problem was encountered when parsing arguments documentation... exiting." << endl;
		exit(-1);
	}

	return true;
}

/*
 * read args from command line and set global variables
 * usage: ./nbody -n nBodies -i nIterations [-v] [-w]
 * */
void argsReader2(int argc, char** argv)
{
	map<string, string> reqArgs, faculArgs, docArgs;
	ArgumentsReader argsReader(argc, argv);

	reqArgs  ["n"]     = "nBodies";
	docArgs  ["n"]     = "approximate number of bodies randomly generated.";
	reqArgs  ["i"]     = "nIterations";
	docArgs  ["i"]     = "number of iterations to compute.";
	faculArgs["f"]     = "inputFileName";
	docArgs  ["f"]     = "bodies file name to read (do not use with -n).";
	faculArgs["v"]     = "";
	docArgs  ["v"]     = "activate verbose mode.";
	faculArgs["w"]     = "outputFileName";
	docArgs  ["w"]     = "base name of body file(s) to write (you can put 'data/out/out').";
	faculArgs["h"]     = "";
	docArgs  ["h"]     = "display this help.";
	faculArgs["-help"] = "";
	docArgs  ["-help"] = "display this help.";

	if(argsReader.parseArguments(reqArgs, faculArgs))
	{
		if(argsReader.existArgument("h") || argsReader.existArgument("-help"))
		{
			if(argsReader.parseDocArgs(docArgs))
				argsReader.printUsage();
			else
				cout << "A problem was encountered when parsing arguments documentation... exiting." << endl;
			exit(-1);
		}

		NBodies       = stoi(argsReader.getArgument("n"));
		NIterations   = stoi(argsReader.getArgument("i"));
		InputFileName = "";

		if(argsReader.existArgument("v"))
			Verbose = true;
		if(argsReader.existArgument("w"))
			OutputBaseName = argsReader.getArgument("w");
	}
	else
	{
		if(argsReader.parseDocArgs(docArgs))
			argsReader.printUsage();
		else
			cout << "A problem was encountered when parsing arguments documentation... exiting." << endl;
		exit(-1);
	}
}

int main(int argc, char** argv)
{
	Perf perfIte, perfTotal;

	// read arguments from command line
	if(!argsReader1(argc, argv)) // usage: ./nbody -f fileName -i nIterations [-v] [-w]
		argsReader2(argc, argv); // usage: ./nbody -n nBodies  -i nIterations [-v] [-w]

	// create plan by reading file or generate bodies randomly
	Space<TYPE> *space;
	if(InputFileName.empty())
		space = new Space<TYPE>(NBodies);
	else
		space = new Space<TYPE>(InputFileName);
	NBodies = space->getNBodies();

	// compute MB used for this simulation
	float Mbytes = (11 * sizeof(TYPE) * NBodies) / 1024.f / 1024.f;

	// compute flops per iteration
	unsigned long flopsPerIte = NBodies * ((NBodies * 16) + 16 + 18);

	// display simulation configuration
	cout << "N-body simulation started !" << endl;
	cout << "---------------------------" << endl;
	if(!InputFileName.empty())
		cout << "  -> inputFileName(s) : " << InputFileName  << ".*.dat" << endl;
	else
		cout << "  -> random mode      : on"                             << endl;
	if(!OutputBaseName.empty())
		cout << "  -> outputFileName(s): " << OutputBaseName << ".*.dat" << endl;
	cout <<     "  -> nBodies          : " << NBodies                    << endl;
	cout <<     "  -> nIterations      : " << NIterations                << endl;
	cout <<     "  -> verbose          : " << Verbose                    << endl;
	cout <<     "  -> mem. used        : " << Mbytes         << " MB"    << endl << endl;


	// draw up visualization window
	TYPE *radius = new TYPE[NBodies]; //TODO: think to delete this buffer before exiting code
	for(unsigned long iBody = 0; iBody < NBodies; iBody++)
		radius[iBody] = space->masses[iBody] * 100.0f;

	// initialize visualization of bodies (with spheres in space)
	OGLSpheresVisualization<TYPE> visu("N-body", WinWidth, WinHeight,
	                                   space->positions.x, space->positions.y, space->positions.z, radius,
	                                   NBodies);

	// specify shaders path and compile them
	vector<GLenum> shadersType(3);
	vector<string> shadersFiles(3);
	shadersType[0] = GL_VERTEX_SHADER;   shadersFiles[0] = "src/ogl/shaders/vertex.glsl";
	shadersType[1] = GL_GEOMETRY_SHADER; shadersFiles[1] = "src/ogl/shaders/geometry.glsl";
	shadersType[2] = GL_FRAGMENT_SHADER; shadersFiles[2] = "src/ogl/shaders/fragment.glsl";

	visu.compileShaders(shadersType, shadersFiles);

	cout << endl << "Simulation started..." << endl;

	// write initial bodies into file
	if(!OutputBaseName.empty())
	{
		std::string outputFileName = OutputBaseName + ".0.dat";
		space->writeIntoFile(outputFileName);
	}

	unsigned long iIte;
	for(iIte = 1; iIte <= NIterations && !visu.windowShouldClose(); iIte++)
	{
		// refresh display in OpenGL window
		visu.refreshDisplay();

		perfIte.start();
		//-----------------------------//
		//-- Simulation computations --//
		space->computeBodiesAcceleration();
		space->findTimeStep();
		space->updateBodiesPositionAndSpeed();
		//-- Simulation computations --//
		//-----------------------------//
		perfIte.stop();
		perfTotal += perfIte;

		if(Verbose)
			cout << "Processing step " << iIte << " took " << perfIte.getElapsedTime() << " ms "
				 << "(" << perfIte.getGflops(flopsPerIte) << " Gflop/s)." << endl;

		// write iteration results into file
		if(!OutputBaseName.empty())
		{
			std::string outputFileName = OutputBaseName + "." + to_string(iIte) + ".dat";
			space->writeIntoFile(outputFileName);
		}

		iIte++;
	}
	cout << "Simulation ended." << endl << endl;

	cout << "Entire simulation took " << perfTotal.getElapsedTime() << " ms "
	     << "(" << perfTotal.getGflops(flopsPerIte * iIte) << " Gflop/s)" << endl;

	delete[] radius;
	delete   space;

	return EXIT_SUCCESS;
}
