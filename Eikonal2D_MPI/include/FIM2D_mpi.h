#ifndef FIM2D_MPI_H
#define FIM2D_MPI_H

#include <mpi.h>
#include <stdbool.h>

#include "../../Eikonal2D/include/SolveEikonal2D.h"
#include "../../Eikonal2D/include/FIM2D.h"

typedef struct{
    MPI_Comm comm;
    MPI_Comm exch_comm; // Communicateur pour les échanges de recouvrement
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
    int i_start_overlap;
    int i_end_overlap;
    int n_overlap;

    int top_ghost;  // Nombre de lignes de recouvrement en haut
    int bottom_ghost;   // Nombre de lignes de recouvrement en bas
    int up_rank;    // Rang du voisin du dessus
    int down_rank;  // Rang du voisin du dessous
} MPIDomain;



typedef struct{
    int n;
    int m;
    double h;

    // Sources globales
    int nsources;
    int* src_i;
    int* src_j;

    // Sources locales
    int nsources_overlap;
    int* src_i_overlap;
    int* src_j_overlap;

    int nwalls;
    int* wall_c1;
    int* wall_c2;
    int* wall_r1;
    int* wall_r2;

    bool valid;
}Config2;





// Crée la topologie 1D en bandes
MPIDomain* topology_create(int n, int m, double h, int overlap);

// Libère la topologie MPI
void topology_free(MPIDomain* domain);

// Echange la bande de recouvrement avec les voisins en haut et en bas, et applique le minimum sur T
// Retourne 1 si au moins une valeur a été améliorée pour le processus et 0 sinon
int exchange_overlap(MPIDomain* domain, EikonalGrid* g_processus, int* changed_cells);

// Propagation de l'onde à partir des mailles start sur autant de pixels que la valeur du recouvrement
// (Par exemple, si on a un recouvrement de 3 pixels alors chaque sous domain MPI propage l'onde sur 3 pixels)
// depth[k] permet de connaitre la distance entre la maille k et la source
// frontier[k] vaut 1 si la maille k a été parcourue lors du dernier tour de la FIM (elle deviendra donc une maille de départ lors du prochain appel à local_propagate), 0 sinon
void local_propagate(EikonalGrid* g_processus, const int* start, int overlap, double epsilon, int* depth, int* frontier);


// FIM avec utilisation du MPI
void fim_solve_mpi(MPIDomain* domain, EikonalGrid* g_processus, Config2* cfg, int* start, double espilon, int nb_cycles);

#endif /* FIM2D_MPI_H */