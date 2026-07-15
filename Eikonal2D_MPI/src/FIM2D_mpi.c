#include "../include/FIM2D_mpi.h"

#include <stdlib.h>
#include <stdio.h>

#define TAG_BAS 100
#define TAG_HAUT 1000

MPIDomain* topology_create(int n, int m, double h, int overlap){
    MPIDomain* domain = (MPIDomain*)malloc(sizeof(MPIDomain));
    if (!domain)
        return NULL;

    // Overlap entre 1 et 3
    if (overlap < 1)
        overlap = 1;
    if (overlap > 3)
        overlap = 3;
    
    // Topologie cartésienne 1D (sous forme de bandes)
    int ndims=1;
    int nproc_per_dim[1] = {0};
    int periods_per_dim[1] = {0};
    int rank, nproc;

    MPI_Comm_rank(MPI_COMM_WORLD, &(rank));
    MPI_Comm_size(MPI_COMM_WORLD, &(nproc));

    int err = 0;

    // Création de la topologie cartésienne adaptée au nombre de processus
    nproc_per_dim[1] = nproc;

    MPI_Comm GRID_COMM;

    err = MPI_Cart_create(MPI_COMM_WORLD, ndims, nproc_per_dim, periods_per_dim, 0, &GRID_COMM);

    if (rank == 0){
        printf("\n dims: %d, nbproc: %d \n", ndims, nproc_per_dim[0]);
        fflush(stdout);
    }

    int proc_coords;
    // Récupérer les coordonnées du rang dans la topologie
    err = MPI_Cart_coords(GRID_COMM, rank, ndims, &proc_coords);

    // Calcul des rangs des voisins avec MPI_Cart_shift
    int up_rank;
    int down_rank;
    MPI_Cart_shift(GRID_COMM, 0, 1, &up_rank, &down_rank);

    // Distribution équilibrée des lignes entre les processus
    int Q = n / nproc;
    int R = n % nproc;
    int i_start;
    int n_owned;
    
    if (proc_coords < R){
        n_owned = Q + 1;
        i_start = proc_coords * (Q+1);
    }
    else{
        n_owned = Q;
        i_start = R * (Q+1) + (proc_coords - R) * Q;
    }

    domain->comm = GRID_COMM;
    domain->rank = rank;
    domain->nproc = nproc;
    domain->proc_coord = proc_coords;
    domain->overlap = overlap;
    domain->n = n;
    domain->m = m;
    domain->h = h;
    domain->i_owned_start = i_start;
    domain->i_owned_end = i_start + n_owned - 1;
    domain->n_owned = n_owned;

    if (up_rank != MPI_PROC_NULL)
        domain->top_ghost = overlap;
    else
        domain->top_ghost = 0;

    if (down_rank != MPI_PROC_NULL)
        domain->bottom_ghost = overlap;
    else
        domain->bottom_ghost = 0;


    // Vérification au cas où un sous-domaine est plus petit que le recouvrement
    if (domain->top_ghost > domain->i_owned_start)
        domain->top_ghost = domain->i_owned_start;
    if (domain->bottom_ghost > n)
        domain->bottom_ghost = n;

    domain->i_local_start = domain->i_owned_start - domain->top_ghost;
    domain->i_local_end = domain->i_owned_end + domain->bottom_ghost;
    domain->n_local = domain->i_local_end - domain->i_local_start + 1;
    domain->up_rank = up_rank;
    domain->down_rank = down_rank;

    return domain;
}


void topology_free(MPIDomain* domain){
    if (!domain)
        return;

    free(domain);
}