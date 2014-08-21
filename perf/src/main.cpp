/*
 * Do not remove.
 * MPI/OpenMP training courses
 * Adrien Cassagne, ASA - CINES, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#include <cmath>
#include <string>
#include <cassert>
#include <iostream>

#include "utils/Perf.h"
#include "utils/ArgumentsReader.h"

#include "Space.h"

using namespace std;

string        InputFileName;
string        OutputBaseName;
unsigned long NBodies;
unsigned long NIterations;
bool          Verbose = false;

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
		NIterations   = atoi(argsReader.getArgument("i").c_str());

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

		NBodies       = atoi(argsReader.getArgument("n").c_str());
		NIterations   = atoi(argsReader.getArgument("i").c_str());
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

	// compute MB used for this simulation
	float Mbytes = (11 * sizeof(TYPE) * NBodies) / 1024.f / 1024.f;

	// display simulation configuration
	cout << "N-body simulation started !" << endl;
	cout << "---------------------------" << endl;
	if(!InputFileName.empty())
		cout << "  -> inputFileName(s) : " << InputFileName  << ".*.dat" << endl;
	else
		cout << "  -> random mode      : on"                             << endl;
	if(!OutputBaseName.empty())
		cout << "  -> outputFileName(s): " << OutputBaseName << ".*.dat" << endl;
	cout <<     "  -> nBodies          : " << space->getNBodies()        << endl;
	cout <<     "  -> nIterations      : " << NIterations                << endl;
	cout <<     "  -> verbose          : " << Verbose                    << endl;
	cout <<     "  -> mem. used        : " << Mbytes         << " MB"    << endl << endl;

	// write initial bodies into file
	if(!OutputBaseName.empty())
	{
		std::string outputFileName = OutputBaseName + ".0.dat";
		space->writeIntoFile(outputFileName);
	}

	cout << "Simulation started..." << endl;
	perfTotal.start();
	for(unsigned long iIte = 1; iIte <= NIterations; iIte++) {
		/*******************************/
		/*** Simulation computations ***/
		perfIte.start();

		space->computeBodiesAcceleration();
		space->findTimeStep();
		space->updateBodiesPositionAndSpeed();

		perfIte.stop();
		/*** Simulation computations ***/
		/*******************************/

		if(Verbose)
			cout << "Processing step " << iIte << " took " << perfIte.getElapsedTime() << " ms." << endl;

		// write iteration results into file
		if(!OutputBaseName.empty())
		{
			std::string outputFileName = OutputBaseName + "." + std::to_string(iIte) + ".dat";
			space->writeIntoFile(outputFileName);
		}
	}
	perfTotal.stop();
	cout << "Simulation ended." << endl << endl;

	unsigned long flops = NBodies * NBodies * 12 * NIterations;
	cout << "Entire simulation took " << perfTotal.getElapsedTime() << " ms "
	     << "(" << perfTotal.getGflops(flops) << " Gflop/s)" << endl;

	delete space;

	return EXIT_SUCCESS;
}
