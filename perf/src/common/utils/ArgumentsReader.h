/*!
 * \file    ArgumentsReader.h
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
#ifndef ARGUMENTS_READER_H_
#define ARGUMENTS_READER_H_

#include <string>
#include <vector>
#include <map>

/*!
 * \class  ArgumentsReader
 * \brief  This class allow us to easily manage arguments from the command line.
 */
class ArgumentsReader {
private:
	std::vector<std::string>           m_argv;            /*!< Simple copie des données de "char** argv". */
	std::map<std::string, std::string> m_requireArgs;     /*!< La liste des arguments obligatoires. */
	std::map<std::string, std::string> m_facultativeArgs; /*!< La liste des arguments facultatifs. */
	std::map<std::string, std::string> m_args;            /*!< La liste des arguments et des valeurs de ces derniers (après parsing). */
	std::map<std::string, std::string> m_docArgs;         /*!< La documentation des arguments si l'utilisateur l'a renseignée. */
	std::string                        m_programName;     /*!< Le nom de l'executable du programme. */

public:
    /*!
     *  \brief Constructeur.
     *
     *  Le contructeur prend les fameux "int argc" et "char** argv" de la fonction main :-).
     *
     *  \param argc : Le nombre d'arguments.
     *  \param argv : Le tableau des arguments
     */
	ArgumentsReader(int argc, char** argv);

    /*!
     *  \brief Destructeur.
     *
     *  Le destructeur ne fait rien...
     */
	virtual ~ArgumentsReader();

    /*!
     *  \brief Parse "m_argv".
     *
     *  Parse "m_argv" selon une liste des arguments requis et facultatifs.
     *
     *  \param requireArgs     : Dictionnaire des arguments attendus obligatoires.
     *  \param facultativeArgs : Dictionnaire des arguments attendus facultatifs.
     *
     *  \return Vrai si tous les arguments requis sont bien présents.
     */
	bool parseArguments(std::map<std::string, std::string> requireArgs,
	                    std::map<std::string, std::string> facultativeArgs);

    /*!
     *  \brief Cherche si un agument existe.
     *
     *  \param tag : Tag de l'argument recherché.
     *
     *  \return Vrai si l'argument existe (à utiliser après parseArguments).
     */
	bool existArgument(std::string tag);

    /*!
     *  \brief Retourne la valeur d'un argument.
     *
     *  \param tag : Tag de l'argument recherché.
     *
     *  \return La valeur d'un argument avec son tag (à utiliser après parseArguments).
     */
	std::string getArgument(std::string tag);

    /*!
     *  \brief Définie la documentation pour les arguments traités par le programme.
     *
     *  \param docArgs : Dictionnaire des arguments à documenter.
     *
     *  \return Faux si docArgs ne contient rien ou si un des arguments de docArgs ne correspond pas à m_args
     *  (à utiliser après parseArguments).
     */
	bool parseDocArgs(std::map<std::string, std::string> docArgs);

    /*!
     *  \brief Affiche une aide pour l'utilisation de la commande.
     */
	void printUsage();

private:
    /*!
     *  \brief Retourne vrai si l'argument "m_argv[posArg]" est dans args.
     *
     *  \param args   : Dictionnaire d'arguments.
     *  \param posArg : La position de l'argument recherché dans m_argv[posArg].
     *
     *  \return Vrai si l'argument "m_argv[posArg]" est dans args.
     */
	bool subParseArguments(std::map<std::string, std::string> args,
	                       unsigned short posArg);

    /*!
     *  \brief Clear m_requireArgs, m_facultativeArgs, m_args and m_docArgs.
     */
	void clearArguments();

};

#endif /* ARGUMENTS_READER_H_ */
