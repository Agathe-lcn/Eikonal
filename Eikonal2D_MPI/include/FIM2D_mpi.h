#ifndef FIM2D_MPI_H
#define FIM2D_MPI_H

#include <mpi.h>
#include "../../Eikonal2D/include/SolveEikonal2D.h"
#include "../../Eikonal2D/include/FIM2D.h"

typedef struct{
    MPI_Comm comm;
    int rank;
    int nproc;
    int proc_coord; // Numéro de la bande du processus
    int overlap;    // Taille du recouvrement

    // Domaine global
    int n;
    int m;
    double h;

    // Domaine sans recouvrement
    int i_owned_start;
    int i_owned_end;
    int n_owned;

    // Domaine avec recouvrement
    int i_local_start;
    int i_local_end;
    int n_local;

    int top_ghost;  // Nombre de lignes de recouvrement en haut
    int bottom_ghost;   // Nombre de lignes de recouvrement en bas
    int up_rank;    // Rang du voisin du dessus
    int down_rank;  // Rang du voisin du dessous
} MPIDomain;

// Crée la topologie 1D en bandes
MPIDomain* topology_create(int n, int m, double h, int overlap);

#endif /* FIM2D_MPI_H */