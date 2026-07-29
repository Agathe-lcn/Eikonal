#include "../include/FIM2D_mpi.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define TAG 0

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

    // Test
    //printf("[rank %d]: n = %d, m = %d, h = %f, overlap = %d, nproc = %d \n", rank, n, m, h, overlap, nproc);
    //fflush(stdout);

    int err = 0;

    // Création de la topologie cartésienne adaptée au nombre de processus
    nproc_per_dim[0] = nproc;

    MPI_Comm GRID_COMM;

    err = MPI_Cart_create(MPI_COMM_WORLD, ndims, nproc_per_dim, periods_per_dim, 0, &GRID_COMM);

    // Test
    //printf("[rank %d]: MPI_Carte_create ->err = %d \n", rank, err);
    //fflush(stdout);

    int proc_coords;
    // Récupérer les coordonnées du rang dans la topologie
    err = MPI_Cart_coords(GRID_COMM, rank, ndims, &proc_coords);

    // Test 
    //printf("[rank %d] proc_coords=%d (err=%d)\n", rank, proc_coords, err);
    //fflush(stdout);

    // Calcul des rangs des voisins avec MPI_Cart_shift
    int up_rank;
    int down_rank;
    MPI_Cart_shift(GRID_COMM, 0, 1, &up_rank, &down_rank);

    // Test
    //printf("[rank %d] up_rank=%d down_rank=%d\n", rank, up_rank, down_rank);
    //fflush(stdout);

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

    // Test
    //printf("[rank %d] Q=%d R=%d n_owned=%d i_start=%d\n", rank, Q, R, n_owned, i_start);
    //fflush(stdout);

    domain->comm = GRID_COMM;
    MPI_Comm_dup(GRID_COMM, &domain->exch_comm);
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

    domain->i_start_overlap = domain->i_owned_start - domain->top_ghost;
    domain->i_end_overlap = domain->i_owned_end + domain->bottom_ghost;
    domain->n_overlap = domain->i_end_overlap - domain->i_start_overlap + 1;
    domain->up_rank = up_rank;
    domain->down_rank = down_rank;

    // Test
    //printf("[rank %d] top_ghost=%d bottom_ghost=%d i_owned=[%d,%d] i_overlap=[%d,%d] n_overlap=%d\n", rank, domain->top_ghost, domain->bottom_ghost, domain->i_owned_start, domain->i_owned_end, domain->i_start_overlap, domain->i_end_overlap, domain->n_overlap);
    //fflush(stdout);

    return domain;
}


void topology_free(MPIDomain* domain){
    if (!domain)
        return;

    MPI_Comm_free(&domain->comm);
    MPI_Comm_free(&domain->exch_comm);
    free(domain);
}


int exchange_overlap(MPIDomain* domain, EikonalGrid* g_processus, int* changed_cells){
    int m = domain->m;
    int overlap = domain->overlap;
    int band = 2*overlap;   // Band est la largeur de la zone en commun entre un processus et son voisin
    int any_change = 0;

    printf("[rank %d] exchange_overlap : m=%d overlap=%d band=%d n_overlap=%d\n", domain->rank, m, overlap, band, domain->n_overlap);
    fflush(stdout);

    // On remet tout le tableau à zéro
    memset(changed_cells, 0, domain->n_overlap*m*sizeof(int));

    MPI_Request reqs[4];

    double* send_up = NULL;
    double* send_down = NULL;
    double* recv_up = NULL;
    double* recv_down = NULL;


    // Pour l'échange avec le voisin du dessus, on vérifie que le processus n'est pas tout en haut, que le recouvrement est strictement positif et que son sous-domaine possède suffisamment de lignes pour l'échange de données
    int do_up = domain->up_rank != MPI_PROC_NULL && overlap > 0 && domain->n_overlap >= band;
    if (do_up){
        send_up = (double*)malloc(band * m * sizeof(double));
        recv_up = (double*)malloc(band * m * sizeof(double));
        memcpy(send_up, g_processus->T, band * m * sizeof(double));

        printf("[rank %d] Sendrecv avec up_rank=%d (%d doubles)\n", domain->rank, domain->up_rank, band*m);
        fflush(stdout);

        MPI_Sendrecv(send_up, band*m, MPI_DOUBLE, domain->up_rank, TAG, recv_up, band*m, MPI_DOUBLE, domain->up_rank, TAG, domain->exch_comm, MPI_STATUS_IGNORE);
    }


    // Pour l'échange avec le voisin du dessous, on vérifie que le processus n'est pas tout en bas, que le recouvrement est strictement positif et que son sous domaine possède suffisamment de lignes pour l'échange de données
    int do_down = domain->down_rank != MPI_PROC_NULL && overlap > 0 && domain->n_overlap >= band;
    if (do_down){
        int beginning = (domain->n_overlap - band)*m;
        send_down = (double*)malloc(band * m * sizeof(double));
        recv_down = (double*)malloc(band * m * sizeof(double));
        memcpy(send_down, g_processus->T + beginning, band * m * sizeof(double));

        printf("[rank %d] Sendrecv avec down_rank=%d (%d doubles)\n", domain->rank, domain->down_rank, band*m);
        fflush(stdout);

        MPI_Sendrecv(send_down, band*m, MPI_DOUBLE, domain->down_rank, TAG, recv_down, band*m, MPI_DOUBLE, domain->down_rank, TAG, domain->exch_comm, MPI_STATUS_IGNORE);
    }


    // Test
    printf("[rank %d] échanges Send/Recv termines\n", domain->rank);
    fflush(stdout);



    // Comparaison des valeurs de T entre le processus courant et son voisin du haut pour garder le minimum à chaque cellule
    if (recv_up){
        for (int k = 0; k < band*m; k++){
            double old_T = g_processus->T[k];
            double new_T = recv_up[k];
            if (new_T < old_T - EIKONAL_EPS){
                g_processus->T[k] = new_T;
                any_change = 1;
                changed_cells[k] = 1;

                // Test
                printf("[rank %d] up: cellule k=%d mise à jour %.4f -> %.4f\n", domain->rank, k, old_T, new_T);
            }
        }

        free(send_up);
        free(recv_up);
    }


    // Comparaison des valeurs de T entre le processus courant et son voisin du bas pour garder le minimum à chaque cellule
    if (recv_down){
        int beginning = (domain->n_overlap - band) * m;
        for (int k = beginning; k < beginning + band*m; k++){
            double old_T = g_processus->T[k];
            double new_T = recv_down[k - beginning];
            if (new_T < old_T - EIKONAL_EPS){
                g_processus->T[k] = new_T;
                any_change = 1;
                changed_cells[k] = 1;

                // Test
                printf("[rank %d] down: cellule k=%d mise à jour %.4f -> %.4f\n", domain->rank, k, old_T, new_T);
            }
        }

        free(send_down);
        free(recv_down);
    }

    // Test
    printf("[rank %d] exchange_overlap terminé, any_change=%d\n", domain->rank, any_change);
    fflush(stdout);

    return any_change;
}


int local_propagate(EikonalGrid* g_processus, const int* start, int overlap, double epsilon, int* depth, int* frontier){
    int n = g_processus->n;
    int m = g_processus->m;
    int ncell = n*m;

    // Initialisation de depth et frontier
    for (int k=0; k < ncell; k++){
        depth[k] = -1;
        frontier[k] = 0;
    }

    // Création de la narrow band
    NodeList* narrow = list_create(ncell);

    // On met les cellules de départ à une profondeur 0
    for (int k=0; k < ncell; k++){
        if (!start[k])
            continue;

        depth[k] = 0;
        if (!list_contains(narrow,k))
            list_push_back(narrow,k);
    }

    while (!list_is_empty(narrow)){
        int index = list_pop_front(narrow);
        int i = index / m;
        int j = index % m;

        double T_old = g_processus->T[index];
        g_processus->T[index] = eikonal_solve_local(g_processus, i, j);
        double diff = fabs(g_processus->T[index] - T_old);

        if (diff <= epsilon){   // Convergence
            // On parcourt les voisins qui sont dans la zone de propagation voulue
            int depth_neighbor = depth[index] + 1;
            if (depth_neighbor <= overlap){
                int neighbors[4][2] = {{i-1, j}, {i+1, j}, {i, j-1}, {i, j+1}};
                for (int k=0; k < 4; k++){
                    int ni = neighbors[k][0];
                    int nj = neighbors[k][1];
                    if (ni >= 0 && ni < n && nj >= 0 && nj < m){
                        int index_neighbor = ni * m + nj;
                        
                        // On met à jour la valeur de la profondeur du voisin
                        if (depth[index_neighbor] < 0)
                            depth[index_neighbor] = index_neighbor;

                        double T_neighbor_new = eikonal_solve_local(g_processus, ni, nj);
                        if (T_neighbor_new < g_processus->T[index_neighbor]){
                            g_processus->T[index_neighbor] = T_neighbor_new;
                            if (!list_contains(narrow, index_neighbor))
                                list_push_back(narrow, index_neighbor);
                        }
                    }
                }
            }
        }
        else
            list_push_front(narrow, index);
    }

    // Les mailles qui ont été atteintes à la profondeur maximale (c'est-à-dire qui sont à la "frontière") sont ajoutées au tableau frontier pour devenir les points de départ de la propagation suivante
    for (int k=0; k < ncell; k++){
        if (depth[k] == overlap)
            frontier[k] = 1;
    }

    // Nettoyage
    list_free(narrow);
}



void fim_solve_mpi(MPIDomain* domain, EikonalGrid* g_processus, Config2* cfg_processus, double epsilon, int nb_cycles){
    int n = g_processus->n;
    int m = g_processus->m;
    int ncell = n * m;
 
    int* start = (int*)malloc(ncell * sizeof(int));
    int* frontier = (int*)malloc(ncell * sizeof(int));
    int* depth = (int*)malloc(ncell * sizeof(int));
    int* changed_cells = (int*)malloc(domain->n_overlap * domain->m * sizeof(int));

    int cycle = 0;
    while(true){
        // On fait la propagation sur 'overlap' cellules de distance
        local_propagate(g_processus, start, domain->overlap, epsilon, depth, frontier);

        // Communication entre les processus
        int changed = exchange_overlap(domain, g_processus, changed_cells);

        // Les mailles de départ du prochain cycle sont celles sur lequelles on s'est arrêté au cycle précédent et les mailles qui ont été modifiées pendant la communication
        memset(start, 0, ncell * sizeof(int));
        int continue_local = 0;
        for (int k = 0; k < ncell; k++){
            if (frontier[k]){
                start[k] = 1;
                continue_local = 1;
            }
            if (changed_cells[k]){
                start[k] = 1;
                continue_local = 1;
            }
        }

        cycle++;
        if (0 < nb_cycles <= cycle)
            break;

        // Si tous les processus n'ont plus de travail alors on arrête tout
        int continue_global = 0;
        MPI_Allreduce(&continue_local, &continue_global, 1, MPI_INT, MPI_MAX, domain->comm);
        if (!continue_global)
            break;
    }

    free(start);
    free(frontier);
    free(depth);
    free(changed_cells);
}