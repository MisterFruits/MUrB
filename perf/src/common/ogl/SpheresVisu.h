/*!
 * \file    SpheresVisu.h
 * \brief   Abstract class for spheres visualization.
 * \author  A. Cassagne
 * \date    2014
 *
 * \section LICENSE
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode).
 *
 * \section DESCRIPTION
 * This is the traditional entry file for the code execution.
 */
#ifndef SPHERES_VISU_H_
#define SPHERES_VISU_H_

class SpheresVisu {
protected:
	SpheresVisu(){}
public:
	virtual      ~SpheresVisu(){}
	virtual void refreshDisplay()    = 0;
	virtual bool windowShouldClose() = 0;
};

#endif /* SPHERES_VISU_H_ */
