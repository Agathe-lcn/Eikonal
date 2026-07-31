
#include "../include/FIM2D_io.h"

#include <stdio.h>
#include <string.h>

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

    // Construire le nom du fichier
    if (cycle < 0)
        snprintf(filename, 255, "%s.chkpt", prefix);
    else
        snprintf(filename, 255, "%s-%d.chkpt", prefix, cycle);

    // Le processus 0 écrit l'en-tête 
    if (domain->rank == 0){
        FILE* file = fopen(filename, "wb");
        if (!file){
            fprintf(stderr, "Erreur: Impossible d'ouvrir le fichier %s.\n", filename);
            return MPI_ERR_IO;
        }

        // Ecriture de l'en-tête: n, m et h
        fwrite(&domain->n, sizeof(int), 1, file);
        fwrite(&domain->m, sizeof(int), 1, file);
        fwrite(&domain->h, sizeof(double), 1, file);

        fclose(file);
    }

    // Barrière pour s'assurer que l'en-tête est bien écrite
    MPI_Barrier(fimio_comm);

    // Chaque processus écrit ses données
    err = MPI_File_open(fimio_comm, filename, amode, info, &fh);

    if (err != MPI_SUCCESS){
        fprintf(stderr, "Error opening %s.\n", filename);
        return err;
    }

    // Calcul de l'offset
    myfileoffset = 16 + (MPI_Offset)domain->i_owned_start * domain->m *sizeof(double);

    err = MPI_File_write_at_all(fh, myfileoffset, g_processus->T, domain->n_overlap * domain->m, MPI_DOUBLE, MPI_STATUS_IGNORE);

    err = MPI_File_close(&fh);

    return err;
}