/*
 * Do not remove.
 * MPI/OpenMP training courses
 * Adrien Cassagne, ASA - CINES, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifdef NBODY_DOUBLE
#define TYPE double
#else
#define TYPE float
#endif

#include <map>
#include <cmath>
#include <string>
#include <vector>
#include <cassert>
#include <iostream>
using namespace std;

#include "ogl/OGLSpheresVisualization.h"

#include "utils/Perf.h"
#include "utils/ArgumentsReader.h"

#include "Space.h"

string        InputFileName;
string        OutputBaseName;
unsigned long NBodies;
unsigned long NIterations;
bool          Verbose = false;
TYPE          Dt      = 3600; //in sec, here 3600 sec = 1 hour

int WinWidth  = 800;
int WinHeight = 600;

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
	faculArgs["-dt"]   = "timeStep";
	docArgs  ["-dt"]   = "select a fixed time step in second.";

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
		if(argsReader.existArgument("-dt"))
			Dt = stof(argsReader.getArgument("-dt"));
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
	faculArgs["-dt"]   = "timeStep";
	docArgs  ["-dt"]   = "select a fixed time step in second.";

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
		if(argsReader.existArgument("-dt"))
			Dt = stof(argsReader.getArgument("-dt"));
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

string strDate(TYPE timestamp)
{
	unsigned int days;
	unsigned int hours;
	unsigned int minutes;
	TYPE rest;

	days = timestamp / (24 * 60 * 60);
	rest = timestamp - (days * 24 * 60 * 60);

	hours = rest / (60 * 60);
	rest = rest - (hours * 60 * 60);

	minutes = rest / 60;
	rest = rest - (minutes * 60);

	return to_string(days)    + "d " +
	       to_string(hours)   + "h " +
	       to_string(minutes) + "m " +
	       to_string(rest)    + "s";
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
	float flopsPerIte = (float) NBodies * (float) (NBodies * 17);

	// display simulation configuration
	cout << "N-body simulation started !" << endl;
	cout << "---------------------------" << endl;
	if(!InputFileName.empty())
		cout << "  -> inputFileName    : " << InputFileName              << endl;
	else
		cout << "  -> random mode      : on"                             << endl;
	if(!OutputBaseName.empty())
		cout << "  -> outputFileName(s): " << OutputBaseName << ".*.dat" << endl;
	cout <<     "  -> nBodies          : " << NBodies                    << endl;
	cout <<     "  -> nIterations      : " << NIterations                << endl;
	cout <<     "  -> verbose          : " << Verbose                    << endl;
#ifdef NBODY_DOUBLE
	cout <<     "  -> precision        : double"                         << endl;
#else
	cout <<     "  -> precision        : simple"                         << endl;
#endif
	cout <<     "  -> mem. used        : " << Mbytes         << " MB"    << endl << endl;

	// initialize visualization of bodies (with spheres in space)
	OGLSpheresVisualization<TYPE> visu("N-body", WinWidth, WinHeight,
	                                   space->positions.x, space->positions.y, space->positions.z, space->radiuses,
	                                   NBodies);

	cout << endl << "Simulation started..." << endl;

	// write initial bodies into file
	if(!OutputBaseName.empty())
	{
		std::string outputFileName = OutputBaseName + ".0.dat";
		space->writeIntoFile(outputFileName);
	}

	// constant timestep (easier for the visualization)
	space->setDtConstant(Dt);

	TYPE physicTime = 0.0;
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

		physicTime += space->getDt();

		if(Verbose)
			cout << "Processing step " << iIte << " took " << perfIte.getElapsedTime() << " ms "
				 << "(" << perfIte.getGflops(flopsPerIte)  << " Gflop/s), "
				 << "physic time: " << strDate(physicTime) << endl;

		// write iteration results into file
		if(!OutputBaseName.empty())
		{
			std::string outputFileName = OutputBaseName + "." + to_string(iIte) + ".dat";
			space->writeIntoFile(outputFileName);
		}
	}
	cout << "Simulation ended." << endl << endl;

	cout << "Entire simulation took " << perfTotal.getElapsedTime() << " ms "
	     << "(" << perfTotal.getGflops(flopsPerIte * iIte) << " Gflop/s)" << endl;

	delete space;

	return EXIT_SUCCESS;
}
