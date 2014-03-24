// create MPI data structure
MPI_Datatype bodyMPI;
body b;
int blockLengths[3] = {1, 1, 1};
MPI_Aint displacements[3];
displacements[0] = (char*) &b.mass - (char*) &b; // char* because sizeof(char) = 1
displacements[1] = (char*) &b.posX - (char*) &b; // char* because sizeof(char) = 1
displacements[2] = (char*) &b.posY - (char*) &b; // char* because sizeof(char) = 1
MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
MPI_Type_create_struct(3, blockLengths, displacements, types, &bodyMPI);

// commit MPI data structure
MPI_Type_commit(&bodyMPI);

// computations and use of MPI structure

// release structure
MPI_Type_free(&bodyMPI);
