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

    EikonalGrid g;
    g.n = domain->n_overlap;
    g.m = domain->m;
    g.h = domain->h;
    g.T = (double*)malloc((size_t)g.n * g.m * sizeof(double));
    g.F = NULL;

    // On met le rang comme valeur de T (non réaliste mais sert d'exemple)
    for (int k = 0; k < g.n * g.m; k++)
        g.T[k] = (double)domain->rank;

    printf("[rank %d] T avant = %.1f partout\n", domain->rank, (double)domain->rank);
    fflush(stdout);

    int* changed_cells = (int*)malloc((size_t)domain->n_overlap * m * sizeof(int));
    int changed = exchange_overlap(domain, &g, changed_cells);

    printf("[rank %d] changed = %d\n", domain->rank, changed);
    for (int i = 0; i < g.n; i++)
        printf("[rank %d] ligne %d : T = %.1f\n", domain->rank, i, g.T[i * m]);
    fflush(stdout);

    free(g.T);
    free(changed_cells);
    topology_free(domain);


    MPI_Finalize();
    return 0;
}