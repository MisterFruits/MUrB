/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifdef NBODY_DOUBLE
using floatType = double;
#else
using floatType = float;
#endif

#include <map>
#include <cmath>
#include <string>
#include <vector>
#include <cassert>
#include <iostream>
using namespace std;

#include "ogl/OGLSpheresVisuInst.h"
#include "ogl/OGLSpheresVisuGS.h"
#include "ogl/OGLSpheresVisuNo.h"

#include "utils/Perf.h"
#include "utils/ArgumentsReader.h"

#include "core/Bodies.h"
#include "core/SimulationNBody.h"
#include "core/v1/SimulationNBodyV1.h"
#include "core/v1/SimulationNBodyV1CB.h"
#include "core/v1/SimulationNBodyV1Vectors.h"
#include "core/v1/SimulationNBodyV1Intrinsics.h"
#include "core/v2/SimulationNBodyV2.h"
#include "core/v2/SimulationNBodyV2CB.h"
#include "core/v2/SimulationNBodyV2Vectors.h"

#ifdef USE_MPI
#include <mpi.h>
#include "core/v1mpi/SimulationNBodyV1MPI.h"
#else
#ifndef NO_MPI
#define NO_MPI
namespace MPI
{
	inline void Init()     {}
	inline void Finalize() {}

	class Comm {
	public:
		                   Comm     (       ) {            };
		virtual            ~Comm    (       ) {            };
		static inline int  Get_rank (       ) { return 0;  };
		static inline int  Get_size (       ) { return 1;  };
		static inline void Abort    (int val) { exit(val); };
	};

	Comm COMM_WORLD;
}
#endif
#endif

/* global variables */
string        RootInputFileName;
string        RootOutputFileName;
unsigned long NBodies;
unsigned long NIterations;
unsigned int  ImplId     = 10;
bool          Verbose    = false;
bool          GSEnable   = false;
bool          VisuEnable = true;
bool          DtVariable = false;
floatType     Dt         = 3600; //in sec, 3600 sec = 1 hour
unsigned int  WinWidth   = 800;
unsigned int  WinHeight  = 600;

/*
 * read args from command line and set global variables
 * usage: ./nbody -n nBodies  -i nIterations [-v] [-w] ...
 * usage: ./nbody -f fileName -i nIterations [-v] [-w] ...
 * */
void argsReader(int argc, char** argv)
{
	map<string, string> reqArgs1, reqArgs2, faculArgs, docArgs;
	ArgumentsReader argsReader(argc, argv);

	reqArgs1 ["n"]     = "nBodies";
	docArgs  ["n"]     = "the number of bodies randomly generated.";
	reqArgs1 ["i"]     = "nIterations";
	docArgs  ["i"]     = "the number of iterations to compute.";

	reqArgs2 ["f"]     = "rootInputFileName";
	docArgs  ["f"]     = "the root file name of the body file(s) to read, do not use with -n "
	                     "(you can put 'data/in/p1/8bodies').";
	reqArgs2 ["i"]     = "nIterations";

	faculArgs["v"]     = "";
	docArgs  ["v"]     = "enable verbose mode.";
	faculArgs["w"]     = "rootOutputFileName";
	docArgs  ["w"]     = "the root file name of the body file(s) to write (you can put 'data/out/bodies').";
	faculArgs["h"]     = "";
	docArgs  ["h"]     = "display this help.";
	faculArgs["-help"] = "";
	docArgs  ["-help"] = "display this help.";
	faculArgs["-dt"]   = "timeStep";
	docArgs  ["-dt"]   = "select a fixed time step in second (default is " + to_string(Dt) + " sec).";
	faculArgs["-gs"]   = "";
	docArgs  ["-gs"]   = "enable geometry shader for visu, "
	                     "this is faster than the standard way but not all GPUs can support it.";
	faculArgs["-ww"]   = "winWidth";
	docArgs  ["-ww"]   = "the width of the window in pixel (default is " + to_string(WinWidth) + ").";
	faculArgs["-wh"]   = "winHeight";
	docArgs  ["-wh"]   = "the height of the window in pixel (default is " + to_string(WinHeight) + ").";
	faculArgs["-nv"]   = "";
	docArgs  ["-nv"]   = "no visualization (disable visu).";
	faculArgs["-vdt"]   = "";
	docArgs  ["-vdt"]   = "enable variable time step.";
	faculArgs["-im"]   = "ImplId";
	docArgs  ["-im"]   = "code implementation id (value should be 10, 11, 12, 13, 20, 21, 22 or 100).";

	if(argsReader.parseArguments(reqArgs1, faculArgs))
	{
		NBodies           = stoi(argsReader.getArgument("n")) / MPI::COMM_WORLD.Get_size();
		NIterations       = stoi(argsReader.getArgument("i"));
		RootInputFileName = "";
	}
	else if(argsReader.parseArguments(reqArgs2, faculArgs))
	{
		RootInputFileName = argsReader.getArgument("f");
		NIterations       = stoi(argsReader.getArgument("i"));
	}
	else
	{
		if(argsReader.parseDocArgs(docArgs))
			argsReader.printUsage();
		else
			cout << "A problem was encountered when parsing arguments documentation... exiting." << endl;
		exit(-1);
	}

	if(argsReader.existArgument("h") || argsReader.existArgument("-help"))
	{
		if(argsReader.parseDocArgs(docArgs))
			argsReader.printUsage();
		else
			cout << "A problem was encountered when parsing arguments documentation... exiting." << endl;
		exit(-1);
	}

	if(argsReader.existArgument("v"))
		Verbose = true;
	if(argsReader.existArgument("w"))
		RootOutputFileName = argsReader.getArgument("w");
	if(argsReader.existArgument("-dt"))
		Dt = stof(argsReader.getArgument("-dt"));
	if(argsReader.existArgument("-gs"))
		GSEnable = true;
	if(argsReader.existArgument("-ww"))
		WinWidth = stoi(argsReader.getArgument("-ww"));
	if(argsReader.existArgument("-wh"))
		WinHeight = stoi(argsReader.getArgument("-wh"));
	if(argsReader.existArgument("-nv"))
		VisuEnable = false;
	if(argsReader.existArgument("-vdt"))
		DtVariable = true;
	if(argsReader.existArgument("-im"))
		ImplId = stoi(argsReader.getArgument("-im"));
}

template <typename T>
string strDate(T timestamp)
{
	unsigned int days;
	unsigned int hours;
	unsigned int minutes;
	T rest;

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

template <typename T>
SimulationNBody<T>* selectImplementationAndAllocateSimulation()
{
	SimulationNBody<floatType> *simu;

	string inputFileName = RootInputFileName + ".p" + to_string(MPI::COMM_WORLD.Get_rank()) + ".dat";

	switch(ImplId)
	{
		case 10:
			if(!MPI::COMM_WORLD.Get_rank())
				cout << "Selected implementation: V1 - O(n²)" << endl << endl;
			if(RootInputFileName.empty())
				simu = new SimulationNBodyV1<T>(NBodies);
			else
				simu = new SimulationNBodyV1<T>(inputFileName);
			break;
		case 11:
			if(!MPI::COMM_WORLD.Get_rank())
				cout << "Selected implementation: V1 + cache blocking - O(n²)" << endl << endl;
			if(RootInputFileName.empty())
				simu = new SimulationNBodyV1CB<T>(NBodies);
			else
				simu = new SimulationNBodyV1CB<T>(inputFileName);
			break;
		case 12:
			if(!MPI::COMM_WORLD.Get_rank())
				cout << "Selected implementation: V1 + vectors - O(n²)" << endl << endl;
			if(RootInputFileName.empty())
				simu = new SimulationNBodyV1Vectors<T>(NBodies);
			else
				simu = new SimulationNBodyV1Vectors<T>(inputFileName);
			break;
		case 13:
			if(!MPI::COMM_WORLD.Get_rank())
				cout << "Selected implementation: V1 + intrinsics - O(n²)" << endl << endl;
			if(RootInputFileName.empty())
				simu = new SimulationNBodyV1Intrinsics<T>(NBodies);
			else
				simu = new SimulationNBodyV1Intrinsics<T>(inputFileName);
			break;
		case 20:
			if(!MPI::COMM_WORLD.Get_rank())
				cout << "Selected implementation: V2 - O(n²/2)" << endl << endl;
			if(RootInputFileName.empty())
				simu = new SimulationNBodyV2<T>(NBodies);
			else
				simu = new SimulationNBodyV2<T>(inputFileName);
			break;
		case 21:
			if(!MPI::COMM_WORLD.Get_rank())
				cout << "Selected implementation: V2 + cache blocking - O(n²/2)" << endl << endl;
			if(RootInputFileName.empty())
				simu = new SimulationNBodyV2CB<T>(NBodies);
			else
				simu = new SimulationNBodyV2CB<T>(inputFileName);
			break;
		case 22:
			if(!MPI::COMM_WORLD.Get_rank())
				cout << "Selected implementation: V2 + vectors - O(n²/2)" << endl << endl;
			if(RootInputFileName.empty())
				simu = new SimulationNBodyV2Vectors<T>(NBodies);
			else
				simu = new SimulationNBodyV2Vectors<T>(inputFileName);
			break;
#ifndef NO_MPI
		case 100:
			if(!MPI::COMM_WORLD.Get_rank())
				cout << "Selected implementation: V1 MPI - O(n²)" << endl << endl;
			if(RootInputFileName.empty())
				simu = new SimulationNBodyV1MPI<T>(NBodies);
			else
				simu = new SimulationNBodyV1MPI<T>(inputFileName);
			break;
#endif
		default:
			if(!MPI::COMM_WORLD.Get_rank())
				cout << "This code implementation does not exist... Exiting." << endl;
			MPI::Finalize();
			exit(-1);
			break;
	}

	if(ImplId < 100 && MPI::COMM_WORLD.Get_size() > 1)
	{
		if(!MPI::COMM_WORLD.Get_rank())
			cout << "Implementation n°" << ImplId << " is not an MPI implementation... Exiting." << endl;
		MPI::Finalize();
		exit(-1);
	}

	return simu;
}

template <typename T>
OGLSpheresVisu<T>* selectImplementationAndAllocateVisu(SimulationNBody<floatType> *simu)
{
	OGLSpheresVisu<T>* visu;

	if(!MPI::COMM_WORLD.Get_rank())
	{
		if(VisuEnable)
		{
			const floatType *positionsX = simu->getBodies().getPositionsX();
			const floatType *positionsY = simu->getBodies().getPositionsY();
			const floatType *positionsZ = simu->getBodies().getPositionsZ();
			const floatType *radiuses   = simu->getBodies().getRadiuses();

			if(GSEnable) // geometry shader = better performances on dedicated GPUs
				visu = new OGLSpheresVisuGS<floatType>("n-body (geometry shader)", WinWidth, WinHeight,
													   positionsX, positionsY, positionsZ,
													   radiuses,
													   NBodies);
			else
				visu = new OGLSpheresVisuInst<floatType>("n-body (instancing)", WinWidth, WinHeight,
														 positionsX, positionsY, positionsZ,
														 radiuses,
														 NBodies);
			cout << endl;
		}
		else
			visu = new OGLSpheresVisuNo<floatType>();
	}
	else
		visu = new OGLSpheresVisuNo<floatType>();

	return visu;
}

int main(int argc, char** argv)
{
	MPI::Init();

	// read arguments from the command line
	// usage: ./nbody -f fileName -i nIterations [-v] [-w] ...
	// usage: ./nbody -n nBodies  -i nIterations [-v] [-w] ...
	argsReader(argc, argv);

	// create the n-body simulation
	SimulationNBody<floatType> *simu = selectImplementationAndAllocateSimulation<floatType>();
	const unsigned long n = simu->getBodies().getN();
	NBodies = n;

	// get MB used for this simulation
	float Mbytes = simu->getAllocatedBytes() / 1024.f / 1024.f;

	// display simulation configuration
	if(!MPI::COMM_WORLD.Get_rank())
	{
		cout << "n-body simulation configuration:" << endl;
		cout << "--------------------------------" << endl;
		if(!RootInputFileName.empty())
			cout << "  -> input file name(s)    : " << RootInputFileName << ".p*.dat" << endl;
		else
			cout << "  -> random mode           : enable" << endl;
		if(!RootOutputFileName.empty())
			cout << "  -> output file name(s)   : " << RootOutputFileName << ".p*.*.dat" << endl;
		cout <<     "  -> total nb. of bodies   : " << NBodies * MPI::COMM_WORLD.Get_size() << endl;
		cout <<     "  -> nb. of bodies per proc: " << NBodies << endl;
		cout <<     "  -> nb. of iterations     : " << NIterations << endl;
		cout <<     "  -> verbose mode          : " << ((Verbose) ? "enable" : "disable") << endl;
		cout <<     "  -> precision             : " << ((sizeof(floatType) == 4) ? "simple" : "double") << endl;
		cout <<     "  -> mem. allocated        : " << Mbytes << " MB" << endl;
		cout <<     "  -> geometry shader       : " << ((GSEnable) ? "enable" : "disable") << endl;
		cout <<     "  -> time step             : " << ((DtVariable) ? "variable" : to_string(Dt) + " sec") << endl;
		cout <<     "  -> nb. of MPI proc       : " << MPI::COMM_WORLD.Get_size() << endl << endl;
	}

	// initialize visualization of bodies (with spheres in space)
	OGLSpheresVisu<floatType> *visu = selectImplementationAndAllocateVisu<floatType>(simu);

	if(!MPI::COMM_WORLD.Get_rank())
		cout << "Simulation started..." << endl;

	// write initial bodies into file
	if(!RootOutputFileName.empty())
	{
		std::string outputFileName = RootOutputFileName + ".p" + to_string(MPI::COMM_WORLD.Get_rank()) + ".0.dat";
		simu->getBodies().writeIntoFile(outputFileName);
	}

	// time step selection
	if(!DtVariable)
		simu->setDtConstant(Dt);

	// loop over the iterations
	Perf perfIte, perfTotal;
	floatType physicTime = 0.0;
	unsigned long iIte;
	for(iIte = 1; iIte <= NIterations && !visu->windowShouldClose(); iIte++)
	{
		// refresh the display in OpenGL window
		visu->refreshDisplay();

		// simulation computations
		perfIte.start();
		simu->computeOneIteration();
		perfIte.stop();
		perfTotal += perfIte;

		// compute the elapsed physic time
		physicTime += simu->getDt();

		// display the status of this iteration
		if(Verbose && !MPI::COMM_WORLD.Get_rank())
			cout << "Processing step " << iIte << " took " << perfIte.getElapsedTime() << " ms "
			     << "(" << perfIte.getGflops(simu->getFlopsPerIte()) << " Gflop/s), "
			     << "physic time: " << strDate(physicTime) << endl;

		// write iteration results into file
		if(!RootOutputFileName.empty())
		{
			std::string outputFileName = RootOutputFileName + ".p" + to_string(MPI::COMM_WORLD.Get_rank()) + "." + to_string(iIte) + ".dat";
			simu->getBodies().writeIntoFile(outputFileName);
		}
	}
	if(!MPI::COMM_WORLD.Get_rank())
		cout << "Simulation ended." << endl << endl;

	if(!MPI::COMM_WORLD.Get_rank())
		cout << "Entire simulation took " << perfTotal.getElapsedTime() << " ms "
		     << "(" << perfTotal.getGflops(simu->getFlopsPerIte() * (iIte -1)) << " Gflop/s)" << endl;

	// free resources
	delete visu;
	delete simu;

	if(NIterations > (iIte +1))
		MPI::COMM_WORLD.Abort(0);

	MPI::Finalize();

	return EXIT_SUCCESS;
}
