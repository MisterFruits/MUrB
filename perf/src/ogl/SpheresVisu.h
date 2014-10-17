/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
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
