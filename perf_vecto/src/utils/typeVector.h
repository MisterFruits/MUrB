/*
 * Do not remove.
 * Optimization training courses 2014 (CINES)
 * Adrien Cassagne, adrien.cassagne@cines.fr
 * This file is under CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
 */

#ifndef TYPE_VECTOR_H_
#define TYPE_VECTOR_H_

#include "myIntrinsics.h" // needed for intrinsic prototypes describe below (see EXPERIMENTAL)

template <typename T = double>
struct alignas(REQUIRED_ALIGNEMENT) vec_t
{
  T vec_data[VECTOR_SIZE] __attribute__ ((aligned (REQUIRED_ALIGNEMENT)));
};

#endif /* TYPE_VECTOR_H_ */
