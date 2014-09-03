#ifndef TYPE_VECTOR_H_
#define TYPE_VECTOR_H_

#include "myIntrinsics.h" // needed for intrinsic prototypes describe below (see EXPERIMENTAL)

template <typename T = double>
struct alignas(REQUIRED_ALIGNEMENT) vec_t
{
  T vec_data[VECTOR_SIZE];
};

#endif /* TYPE_VECTOR_H_ */
