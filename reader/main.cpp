/*
 * Do not remove.
 * MPI/OpenMP training courses
 * Adrien CASSAGNE, ASA - CINES, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#include <iostream>
#include <unistd.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <sstream>
#include <cassert>
#include <fstream>
#include <SDL/SDL.h>

using namespace std;

struct Point
{
	double mass;
	double x;
	double y;
};

int counter = 0; // Pas beau

/*
 * This is a 32-bit pixel function created with help from this
 * website: http://www.libsdl.org/intro.en/usingvideo.html
 *
 * You will need to make changes if you want it to work with
 * 8-, 16- or 24-bit surfaces.  Consult the above website for
 * more information.
 */
void set_pixel(SDL_Surface *surface, int x, int y, Uint32 pixel)
{
    if(x >= 0 && y >= 0 && x < surface->w && y < surface->h) {
    	Uint8 *target_pixel = (Uint8 *)surface->pixels + y * surface->pitch + x * 4;
    	*(Uint32 *)target_pixel = pixel;
    }
}

/*
 * This is an implementation of the Midpoint Circle Algorithm
 * found on Wikipedia at the following link:
 *
 *   http://en.wikipedia.org/wiki/Midpoint_circle_algorithm
 *
 * The algorithm elegantly draws a circle quickly, using a
 * set_pixel function for clarity.
 */
void draw_circle(SDL_Surface *surface, int n_cx, int n_cy, int radius, Uint32 pixel)
{
    // if the first pixel in the screen is represented by (0,0) (which is in sdl)
    // remember that the beginning of the circle is not in the middle of the pixel
    // but to the left-top from it:

    double error = (double)-radius;
    double x = (double)radius -0.5;
    double y = (double)0.5;
    double cx = n_cx - 0.5;
    double cy = n_cy - 0.5;

    while (x >= y)
    {
        set_pixel(surface, (int)(cx + x), (int)(cy + y), pixel);
        set_pixel(surface, (int)(cx + y), (int)(cy + x), pixel);

        if (x != 0)
        {
            set_pixel(surface, (int)(cx - x), (int)(cy + y), pixel);
            set_pixel(surface, (int)(cx + y), (int)(cy - x), pixel);
        }

        if (y != 0)
        {
            set_pixel(surface, (int)(cx + x), (int)(cy - y), pixel);
            set_pixel(surface, (int)(cx - y), (int)(cy + x), pixel);
        }

        if (x != 0 && y != 0)
        {
            set_pixel(surface, (int)(cx - x), (int)(cy - y), pixel);
            set_pixel(surface, (int)(cx - y), (int)(cy - x), pixel);
        }

        error += y;
        ++y;
        error += y;

        if (error >= 0)
        {
            --x;
            error -= x;
            error -= x;
        }
    }
}

void draw_mass(SDL_Surface *surface, int x, int y, double mass)
{
	int startColor = 100;
	int incr = 10;

	SDL_Color color;
	color.r = startColor;
	color.g = startColor;
	color.b = startColor;

	double falseMass = mass;
	int i = startColor;
	while( (falseMass > 1) && (i < (255 - startColor)) )
	{
		color.r += incr;
		color.g += incr;
		color.b += incr;

		falseMass /= 10;
		i += incr;
	}

	falseMass *= 20;

	draw_circle(surface, x, y, ceil(falseMass), SDL_MapRGB(surface->format, color.r, color.g, color.b));
}

bool readVectorFile(vector<Point>& vect, string& prefixName)
{
	vect.clear();

	stringstream ss;
	ss << prefixName << "." << counter << ".dat";
	string filename = ss.str();

	ifstream vectorFile;
	vectorFile.open(filename.c_str(), std::ios::in);

	if(!vectorFile.is_open())
		return false;

	while(vectorFile.good())
	{
		Point p;
		double trash;

		vectorFile >> p.mass;
		vectorFile >> p.x;
		vectorFile >> p.y;
		vectorFile >> trash;
		vectorFile >> trash;

		p.mass = p.mass/1000.0;

		vect.push_back(p);
	}

	vectorFile.close();

	counter++;
	return true;
}

void display(vector<Point>& vect, SDL_Surface *screen, unsigned step)
{
	SDL_FillRect(screen, NULL, SDL_MapRGB(screen->format, 0, 0, 0));
#ifdef DEBUG
	cout << "---------------------------------------------------" << endl;
	cout << "Step " << step << " (" << vect.size() << " masses): "<< endl;
#endif
	for(unsigned long i=0; i < vect.size(); ++i) {
		draw_mass(screen, vect[i].x, vect[i].y, vect[i].mass);
#ifdef DEBUG
		cout << vect[i].mass << " " << vect[i].x << " " << vect[i].y << endl;
#endif
	}

	SDL_Flip(screen);
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "usage : reader width height prefix_filename" << endl;
		exit(0);
	}

	int xPixelSize = atoi(argv[1]);
	int yPixelSize = atoi(argv[2]);
	string prefixName = string (argv[3]);

	// INIT SDL //////////////////////////////////////
	if (SDL_Init(SDL_INIT_VIDEO) == -1) // Démarrage de la SDL. Si erreur :
	{
		cerr << "Erreur d'initialisation de la SDL : " << SDL_GetError() << endl; // Écriture de l'erreur
		exit(EXIT_FAILURE); // On quitte le programme
	}
	SDL_Surface *screen = NULL;
	screen = SDL_SetVideoMode(xPixelSize, yPixelSize, 32, SDL_HWSURFACE | SDL_DOUBLEBUF);
	SDL_WM_SetCaption("Body Reader", NULL);
	//////////////////////////////////////////////////

	vector<Point> vect;
	int halt = 0;
	SDL_Event event;
	int i=0;
	bool anim = false;
	unsigned int initDelay = 120;
	unsigned int curDelay = initDelay;

	while(readVectorFile(vect, prefixName) && !halt)
	{
		display(vect, screen, i++);
		
		while (!halt)
		{
			if(!anim)
				SDL_WaitEvent(&event);
			else
			{
				SDL_Delay(curDelay);
				SDL_PollEvent(&event);
			}

			if(event.type == SDL_QUIT)
				halt = 1;
			else if(event.type == SDL_KEYDOWN)
			{
				if(event.key.keysym.sym == SDLK_ESCAPE)
				{
					halt = 1;
				}
				else if(event.key.keysym.sym == SDLK_RIGHT)
				{
					anim = true;
					curDelay /= 2;
				}
				else if(event.key.keysym.sym == SDLK_LEFT)
				{
					anim = false;
					curDelay = initDelay;
				}
				else if(event.key.keysym.sym == SDLK_UP)
					break;
			}
			if(anim)
				break;
		}
	}

	SDL_Quit();
	return EXIT_SUCCESS;
}
