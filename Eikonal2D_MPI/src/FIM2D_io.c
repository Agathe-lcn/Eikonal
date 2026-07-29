#include "../include/FIM2D_io.h"

#include <stdio.h>
#include <string.h>

static int FIMIO_Type_create_rowblk(EikonalGrid* g_processus, MPIDomain* domain, MPI_Datatype* newtype){
    int err;
    int n_overlap = domain->n_overlap;
    int m = domain->m;

    MPI_Datatype vectype;
    MPI_Aint disp;

    err = MPI_Type_vector(n_overlap, m, m, MPI_DOUBLE, &vectype);

    if (err != MPI_SUCCESS) 
        return err;

    int len = 1;

    MPI_Get_address(g_processus->T, &disp);
    err = MPI_Type_create_hindexed(1, &len, &disp, vectype, newtype);

    err = MPI_Type_free(&vectype);

    return err;
}

static int FIMIO_Type_create_hdr_rowblk(EikonalGrid* g_processus, MPIDomain* domain, MPI_Datatype *newtype){
    int err;
    int n = domain->n;
    int m = domain->m;
    double h = domain->h;
    int n_overlap = domain->n_overlap;

    int lens[4] = { 1, 1, 1, 1 };
    MPI_Aint disps[4];
    MPI_Datatype types[4];
    MPI_Datatype rowblk;

    FIMIO_Type_create_rowblk(g_processus, domain, &rowblk);

    MPI_Get_address(&n, &disps[0]);
    MPI_Get_address(&m, &disps[1]);
    MPI_Get_address(&h, &disps[2]);
    disps[3] = (MPI_Aint) MPI_BOTTOM ;

    types[0] = MPI_INT;
    types[1] = MPI_INT;
    types[2] = MPI_DOUBLE;
    types[3] = rowblk;

    err = MPI_Type_create_struct(4, lens, disps, types, newtype);

    err = MPI_Type_free(&rowblk);

    return err;
}

static MPI_Comm fimio_comm = MPI_COMM_NULL;

int FIMIO_Init(MPI_Comm comm){
    int err;
    err = MPI_Comm_dup(comm, &fimio_comm);
    return err;
}

int FIMIO_Finalize(void){
    int err = MPI_SUCCESS;
    err = MPI_Comm_free(&fimio_comm);
    return err;
}

int FIMIO_Can_restart(void){
    return 1;
}

int FIMIO_Checkpoint(char* prefix, MPIDomain* domain, EikonalGrid* g_processus, int cycle, MPI_Info info){
    int err = MPI_SUCCESS;
    int amode = MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_UNIQUE_OPEN;

    MPI_File fh;
    MPI_Datatype type;
    MPI_Offset myfileoffset;

    char filename[256];

    snprintf(filename, 255, "%s-%d.chkpt", prefix, cycle);

    err = MPI_File_open(fimio_comm, filename, amode, info, &fh);

    if (err != MPI_SUCCESS){
        fprintf(stderr, "Error opening %s.\n", filename);
        return err;
    }

    if (domain->rank == 0){
        FIMIO_Type_create_hdr_rowblk(g_processus, domain, &type);
        myfileoffset = 0;
    }
    else{
        FIMIO_Type_create_rowblk(g_processus, domain, &type);
        myfileoffset = 3 * (MPI_Offset)sizeof(int) + (MPI_Offset)domain->i_owned_start * domain->m * (MPI_Offset)sizeof(double);
    }

    err = MPI_Type_commit(&type);

    err = MPI_File_write_at_all(fh, myfileoffset, MPI_BOTTOM, 1, type, MPI_STATUS_IGNORE);

    err = MPI_Type_free(&type);

    err = MPI_File_close(&fh);

    return err;
}

