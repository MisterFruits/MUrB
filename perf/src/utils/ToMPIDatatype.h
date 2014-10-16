#ifndef TO_MPI_DATATYPE_H_
#define TO_MPI_DATATYPE_H_

#include <mpi.h>

template<typename T>
class ToMPIDatatype {
public:
	static inline MPI::Datatype value();
};

template<>
inline MPI::Datatype ToMPIDatatype<double>::value() {
	return MPI_DOUBLE;
}

template<>
inline MPI::Datatype ToMPIDatatype<float>::value() {
	return MPI_FLOAT;
}

template<>
inline MPI::Datatype ToMPIDatatype<long>::value() {
	return MPI_LONG;
}

template<>
inline MPI::Datatype ToMPIDatatype<int>::value() {
	return MPI_INTEGER;
}

#endif /* TO_MPI_DATATYPE_H_ */
