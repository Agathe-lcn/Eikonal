#include "../include/FIM2D_mpi.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv){
    MPI_Init(&argc, &argv);

    int n = 20;      // nombre de lignes total
    int m = 10;       // nombre de colonnes
    double h = 0.05;
    int overlap = 2;

    MPIDomain* domain = topology_create(n, m, h, overlap);

    if (domain){
        printf("[rank %d] Domaine cree avec succes\n", domain->rank);
        free(domain);
    } else {
        printf("Erreur : allocation du domaine echouee\n");
    }

    MPI_Finalize();
    return 0;
}