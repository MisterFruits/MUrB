/*
 * ArgumentsReader.h
 *
 *  Created on: Apr 10, 2013
 *      Author: Adrien CASSAGNE
 */

#ifndef ARGUMENTS_READER_H_
#define ARGUMENTS_READER_H_

#include <string>
#include <vector>
#include <map>

class ArgumentsReader {
private:
	/* Simple copie des données de "char** argv". */
	std::vector<std::string>           m_argv;

	/* La liste des arguments obligatoires. */
	std::map<std::string, std::string> m_requireArgs;

	/* La liste des arguments facultatifs. */
	std::map<std::string, std::string> m_facultativeArgs;

	/* La liste des arguments et des valeurs de ces derniers (après parsing). */
	std::map<std::string, std::string> m_args;

	/* La documentation des arguments si l'utilisateur l'a renseignée */
	std::map<std::string, std::string> m_docArgs;

	/* Le nom de l'executable du programme. */
	std::string                        m_programName;

public:
	/*
	 * Le contructeur prend les fameux "int argc" et "char** argv" de la fonction main :-).
	 */
	ArgumentsReader(int argc, char** argv);

	/*
	 * Le destructeur ne fait rien...
	 */
	virtual ~ArgumentsReader();

	/*
	 * Parse "m_argv" selon une liste des arguments requis et facultatifs.
	 * Retourne vrai si tous les arguments requis sont bien présent.
	 */
	bool parseArguments(std::map<std::string, std::string> requireArgs,
	                    std::map<std::string, std::string> facultativeArgs);

	/*
	 * Retourne vrai si l'argument existe (à utiliser après parseArguments).
	 */
	bool existArgument(std::string tag);

	/*
	 * Retourne la valeur d'un argument avec son tag (à utiliser après parseArguments).
	 */
	std::string getArgument(std::string tag);


	/*
	 * Définie la documentation pour les arguments traités par le programme. Retourne faux si docArgs ne contient rien
	 * ou si un des arguments de docArgs ne correspond pas à m_args (à utiliser après parseArguments).
	 */
	bool parseDocArgs(std::map<std::string, std::string> docArgs);

	/*
	 * Affiche une aide pour l'utilisation de la commande.
	 */
	void printUsage();

private:
	/*
	 * Retourne vrai si l'argument "m_argv[posArg]" est dans args.
	 */
	bool subParseArguments(std::map<std::string, std::string> args,
	                       unsigned short posArg);

	/*
	 * Clear m_requireArgs, m_facultativeArgs, m_args and m_docArgs
	 */
	void clearArguments();

};

#endif /* ARGUMENTS_READER_H_ */
