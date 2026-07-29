#ifndef FIM2D_IO_H
#define FIM2D_IO_H

#include "FIM2D_mpi.h"

#include <mpi.h>

int FIMIO_Init(MPI_Comm comm);

int FIMIO_Finalize(void);

int FIMIO_Can_restart(void);

int FIMIO_Checkpoint(char* prefix, MPIDomain* domain, EikonalGrid* g_processus, int cycle, MPI_Info info);

#endif /* FIM2D_IO_H */