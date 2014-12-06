/*!
 * \file    ArgumentsReader.cpp
 * \brief   Command line arguments management class.
 * \author  A. Cassagne
 * \date    2014
 *
 * \section LICENSE
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode).
 *
 * \section DESCRIPTION
 * This is the traditional entry file for the code execution.
 */
#include <iostream>
#include <cassert>

#include "ArgumentsReader.h"

using namespace std;

ArgumentsReader::ArgumentsReader(int argc, char** argv)
	: m_argv(argc)
{
	assert(argc > 0);

	this->m_programName = argv[0];

	for(unsigned short i = 0; i < argc; ++i)
		this->m_argv[i] = argv[i];
}

ArgumentsReader::~ArgumentsReader()
{
}

bool ArgumentsReader::parseArguments(map<string, string> requireArgs, map<string, string> facultativeArgs)
{
	//assert(requireArgs.size() > 0); // useless, it is possible to have no require arguments
	unsigned short int nReqArg = 0;

	this->clearArguments();

	this->m_requireArgs = requireArgs;
	this->m_facultativeArgs = facultativeArgs;

	for(unsigned short i = 0; i < this->m_argv.size(); ++i)
	{
		if(this->subParseArguments(this->m_requireArgs, i))
			nReqArg++;
		this->subParseArguments(this->m_facultativeArgs, i);
	}

	return nReqArg >= requireArgs.size();
}

bool ArgumentsReader::subParseArguments(map<string, string> args, unsigned short posArg)
{
	assert(posArg < this->m_argv.size());

	bool isFound = false;

	map<string, string>::iterator it;
	for(it = args.begin(); it != args.end(); ++it)
	{
		string curArg = "-" + it->first;
		if(curArg == this->m_argv[posArg])
		{
			if(it->second != "")
			{
				if(posArg != (this->m_argv.size() -1))
				{
					this->m_args[it->first] = this->m_argv[posArg +1];
					isFound = true;
				}
			}
			else
			{
				this->m_args[it->first] = "";
				isFound = true;
			}
		}
	}

	return isFound;
}

bool ArgumentsReader::existArgument(std::string tag)
{
	return (this->m_args.find(tag) != this->m_args.end());
}

string ArgumentsReader::getArgument(string tag)
{
	return this->m_args[tag];
}

bool ArgumentsReader::parseDocArgs(std::map<std::string, std::string> docArgs)
{
	bool reVal = true;

	if(docArgs.empty())
		reVal = false;

	map<string, string>::iterator it;
	for(it = this->m_requireArgs.begin(); it != this->m_requireArgs.end(); ++it)
		if(!(docArgs.find(it->first) != docArgs.end()))
			reVal = false;
		else
			this->m_docArgs[it->first] = docArgs[it->first];

	for(it = this->m_facultativeArgs.begin(); it != this->m_facultativeArgs.end(); ++it)
		if(!(docArgs.find(it->first) != docArgs.end()))
			reVal = false;
		else
			this->m_docArgs[it->first] = docArgs[it->first];

	return reVal;
}

void ArgumentsReader::printUsage()
{
	cout << "Usage: " << this->m_programName;

	map<string, string>::iterator it;
	for(it = this->m_requireArgs.begin(); it != this->m_requireArgs.end(); ++it)
		if(it->second != "")
			cout << " -" << it->first << " " << it->second;
		else
			cout << " -" << it->first;

	for(it = this->m_facultativeArgs.begin(); it != this->m_facultativeArgs.end(); ++it)
		if(it->second != "")
			cout << " [-" << it->first << " " << it->second << "]";
		else
			cout << " [-" << it->first << "]";

	cout << endl;

	if(!this->m_docArgs.empty())
	{
		cout << endl;
		for(it = this->m_requireArgs.begin(); it != this->m_requireArgs.end(); ++it)
			if(this->m_docArgs.find(it->first) != this->m_docArgs.end())
				cout << "\t-" << it->first << "\t\t" << this->m_docArgs.find(it->first)->second << endl;

		for(it = this->m_facultativeArgs.begin(); it != this->m_facultativeArgs.end(); ++it)
			if(this->m_docArgs.find(it->first) != this->m_docArgs.end())
				cout << "\t-" << it->first << "\t\t" << this->m_docArgs.find(it->first)->second << endl;
	}
}

void ArgumentsReader::clearArguments()
{
	this->m_requireArgs.clear();
	this->m_facultativeArgs.clear();
	this->m_args.clear();
	this->m_docArgs.clear();
}
