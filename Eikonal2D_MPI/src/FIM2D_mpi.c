#include "../include/FIM2D_mpi.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define TAG 0
#define TEST_MODE 0     // 1 pour désactiver les communications et 0 pour les activer
#define DEBUG_COMM_PHASES 2   // arrêter après la 2e communication MPI

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


static void write_exchange_trace(MPIDomain* domain, int cycle, int communication_count, const char* direction, const char* stage,
                                 int local_row, int local_col, double value, double old_T, double new_T, int updated){
    char filename[256];
    snprintf(filename, sizeof(filename), "mpi_exchange_rank%d.txt", domain->rank);

    FILE* file = fopen(filename, "a");
    if (!file)
        return;

    int global_row = domain->i_start_overlap + local_row;
    fprintf(file,
            "cycle=%d comm=%d rank=%d direction=%s stage=%s local=(%d,%d) global=(%d,%d) value=%.6f old_T=%.6f new_T=%.6f updated=%d\n",
            cycle, communication_count, domain->rank, direction, stage,
            local_row, local_col, global_row, local_col, value, old_T, new_T, updated);
    fclose(file);
}

int exchange_overlap(MPIDomain* domain, EikonalGrid* g_processus, int* changed_cells, int cycle, int communication_count){
    // Si on veut tester sans les communications, alors on ne fait rien
    #if TEST_MODE
        return 0;
    #endif

    int m = domain->m;
    int overlap = domain->overlap;
    int any_change = 0;

    memset(changed_cells, 0, domain->n_overlap*m*sizeof(int));

    double* send_up = NULL;
    double* send_down = NULL;
    double* recv_up = NULL;
    double* recv_down = NULL;

    int do_up = domain->up_rank != MPI_PROC_NULL && domain->top_ghost > 0;
    if (do_up){
        send_up = (double*)malloc(overlap * m * sizeof(double));
        recv_up = (double*)malloc(overlap * m * sizeof(double));
        memcpy(send_up, g_processus->T + domain->top_ghost*m, overlap * m * sizeof(double));

        for (int offset = 0; offset < overlap; offset++){
            int local_row = domain->top_ghost + offset;
            for (int col = 0; col < m; col++){
                int idx = offset * m + col;
                if (cycle == 0 || cycle == 1)
                    write_exchange_trace(domain, cycle, communication_count, "up", "send",
                                     local_row, col, send_up[idx], 0.0, 0.0, 0);
            }
        }

        MPI_Sendrecv(send_up, overlap*m, MPI_DOUBLE, domain->up_rank, TAG, recv_up, overlap*m, MPI_DOUBLE, domain->up_rank, TAG, domain->exch_comm, MPI_STATUS_IGNORE);
    }

    int do_down = domain->down_rank != MPI_PROC_NULL && domain->bottom_ghost > 0;
    if (do_down){
        send_down = (double*)malloc(overlap * m * sizeof(double));
        recv_down = (double*)malloc(overlap * m * sizeof(double));
        int send_offset = (domain->n_overlap - domain->bottom_ghost - overlap) * m;
        memcpy(send_down, g_processus->T + send_offset, overlap * m * sizeof(double));

        for (int offset = 0; offset < overlap; offset++){
            int local_row = domain->n_overlap - domain->bottom_ghost - overlap + offset;
            for (int col = 0; col < m; col++){
                int idx = offset * m + col;
                if (cycle == 0 || cycle == 1)
                    write_exchange_trace(domain, cycle, communication_count, "down", "send",
                                     local_row, col, send_down[idx], 0.0, 0.0, 0);
            }
        }

        MPI_Sendrecv(send_down, overlap*m, MPI_DOUBLE, domain->down_rank, TAG, recv_down, overlap*m, MPI_DOUBLE, domain->down_rank, TAG, domain->exch_comm, MPI_STATUS_IGNORE);
    }

    if (recv_up){
        for (int k = 0; k < overlap*m; k++){
            int local_row = k / m;
            int local_col = k % m;
            int global_row = domain->i_start_overlap + local_row;
            double old_T = g_processus->T[k];
            double new_T = recv_up[k];
            if (cycle == 0 || cycle == 1)
                write_exchange_trace(domain, cycle, communication_count, "up", "recv",
                                 local_row, local_col, new_T, old_T, new_T, 0);
            if (new_T < old_T - EIKONAL_EPS){
                g_processus->T[k] = new_T;
                any_change = 1;
                changed_cells[k] = 1;
                if (cycle == 0 || cycle == 1)
                    write_exchange_trace(domain, cycle, communication_count, "up", "apply",
                                     local_row, local_col, new_T, old_T, new_T, 1);
            } else {
                if (cycle == 0 || cycle == 1)
                    write_exchange_trace(domain, cycle, communication_count, "up", "keep",
                                     local_row, local_col, new_T, old_T, new_T, 0);
            }
        }

        free(send_up);
        free(recv_up);
    }

    if (recv_down){
        int beginning = (domain->n_overlap - overlap) * m;
        for (int k = beginning; k < beginning + overlap*m; k++){
            int local_row = (k / m);
            int local_col = k % m;
            double old_T = g_processus->T[k];
            double new_T = recv_down[k - beginning];
            if (cycle == 0 || cycle == 1)
                write_exchange_trace(domain, cycle, communication_count, "down", "recv",
                                 local_row, local_col, new_T, old_T, new_T, 0);
            if (new_T < old_T - EIKONAL_EPS){
                g_processus->T[k] = new_T;
                any_change = 1;
                changed_cells[k] = 1;
                if (cycle == 0 || cycle == 1)
                    write_exchange_trace(domain, cycle, communication_count, "down", "apply",
                                     local_row, local_col, new_T, old_T, new_T, 1);
            } else {
                if (cycle == 0 || cycle == 1)
                    write_exchange_trace(domain, cycle, communication_count, "down", "keep",
                                     local_row, local_col, new_T, old_T, new_T, 0);
            }
        }

        free(send_down);
        free(recv_down);
    }

    return any_change;
}


int local_propagate(EikonalGrid* g_processus, const int* start, int overlap, double epsilon, int* depth, int* frontier){
    int n = g_processus->n;
    int m = g_processus->m;
    int ncell = n*m;
    int local_changed = 0;

    // TEST
    /*printf("[Processus] === DEBUT local_propagate ===\n");
    printf("[Processus] overlap=%d, epsilon=%f\n", overlap, epsilon);
    int nb_start = 0;
    for (int k = 0; k < ncell; k++) {
        if (start[k]) nb_start++;
    }
    printf("[Processus] start contient %d cellules\n", nb_start);*/

    // Initialisation de depth et frontier
    for (int k=0; k < ncell; k++){
        depth[k] = -1;
        frontier[k] = 0;
    }

    // Création de la narrow band
    NodeList* narrow = list_create(ncell);

    // On met les cellules de départ à une profondeur 0

    // TEST
    //int nb_in_narrow = 0;

    for (int k=0; k < ncell; k++){
        if (start[k]){
            depth[k] = 0;
            list_push_back(narrow,k);
        }

            // TEST
            //nb_in_narrow++;
            //printf("[Processus] Cellule %d ajoutée à narrow (source)\n", k);
    }

    // TEST
    /*printf("[Processus] %d cellules dans narrow\n", nb_in_narrow);
    if (list_is_empty(narrow)) {
        printf("[Processus] ERREUR: narrow est vide ! Aucune source trouvée.\n");
        list_free(narrow);
        return;
    }
    int iterations = 0;*/

    while (!list_is_empty(narrow)){
        // TEST
        //iterations++;

        int index = list_pop_front(narrow);
        int i = index / m;
        int j = index % m;

        // TEST
        //printf("[Processus] Iteration %d: traitement de la cellule %d (%d,%d)\n", iterations, index, i, j);

        if (eikonal_grid_is_obstacle(g_processus, i, j))
            continue;

        double T_old = g_processus->T[index];

        // Si c'est une source alors on ne la recalcule pas
        if (T_old != 0.0){
            g_processus->T[index] = eikonal_solve_local(g_processus, i, j);
        }

        double diff = fabs(g_processus->T[index] - T_old);
        if (diff > epsilon)
            local_changed = 1;

        //TEST
        //printf("[Processus] T_old=%f, T_new=%f, diff=%f\n", T_old, g_processus->T[index], diff);

        if (diff <= epsilon){   // Convergence
            // On parcourt les voisins qui sont dans la zone de propagation voulue
            int current_depth = depth[index];

            // Si on a atteint la profondeur max de l'overlap, alors on ne propage plus aux voisins
            if (current_depth >= overlap){
                frontier[index] = 1;
                continue;
            }

            // TEST
            //printf("[Processus] Convergence! depth_neighbor=%d, overlap=%d\n", depth_neighbor, overlap);
            int neighbors[4][2] = {{i-1, j}, {i+1, j}, {i, j-1}, {i, j+1}};
            for (int k=0; k < 4; k++){
                int ni = neighbors[k][0];
                int nj = neighbors[k][1];
                if (ni >= 0 && ni < n && nj >= 0 && nj < m){
                    if (eikonal_grid_is_obstacle(g_processus, ni, nj))
                        continue;

                    int index_neighbor = ni * m + nj;
                    double T_neighbor_new = eikonal_solve_local(g_processus, ni, nj);
                    
                    if (T_neighbor_new < g_processus->T[index_neighbor] - epsilon){
                        g_processus->T[index_neighbor] = T_neighbor_new;
                        local_changed = 1;

                        // TEST
                        //printf("[Processus] Voisin (%d,%d) mis à jour: %f\n", ni, nj, T_neighbor_new);

                        // On met à jour la profondeur du voisin
                        depth[index_neighbor] = current_depth + 1;

                        if (!list_contains(narrow, index_neighbor)){
                            list_push_back(narrow, index_neighbor);

                            // TEST
                            //printf("[Processus] Voisin ajouté à narrow\n");
                        }
                    }
                }
            }
        }
        else{
            list_push_front(narrow, index);

            // TEST
            //printf("[Processus] Pas de convergence, remis en tête de narrow\n");
        }
    }

    // Les mailles qui ont été atteintes à la profondeur maximale (c'est-à-dire qui sont à la "frontière") sont ajoutées au tableau frontier pour devenir les points de départ de la propagation suivante
    
    //TEST
    //int nb_frontier = 0;

    for (int k=0; k < ncell; k++){
        if (depth[k] == overlap){
            frontier[k] = 1;

            // TEST
            //nb_frontier++;
        }
    }

    // TEST
    //printf("[Processus] %d cellules en frontière (profondeur=%d)\n", nb_frontier, overlap);
    //printf("[Processus] === FIN local_propagate (%d itérations) ===\n", iterations);

    // Nettoyage
    list_free(narrow);
    return local_changed;
}


static void write_debug_trace(MPIDomain* domain, int cycle, const char* phase, int local_changed, int changed, int continue_local, int continue_global, int communication_count){
    char filename[256];
    snprintf(filename, sizeof(filename), "mpi_debug_rank%d.txt", domain->rank);

    FILE* file = fopen(filename, "a");
    if (!file)
        return;

    fprintf(file,
            "cycle=%d phase=%s comm=%d rank=%d local_changed=%d changed=%d continue_local=%d continue_global=%d\n",
            cycle, phase, communication_count, domain->rank, local_changed, changed, continue_local, continue_global);
    fclose(file);
}

static void write_cycle_summary(MPIDomain* domain, int cycle, int communication_count, int local_changed, int changed, int continue_local, int continue_global){
    char filename[256];
    snprintf(filename, sizeof(filename), "mpi_cycle_summary.txt");

    FILE* file = fopen(filename, "a");
    if (!file)
        return;

    fprintf(file,
            "cycle=%d comm=%d rank=%d local_changed=%d changed=%d continue_local=%d continue_global=%d\n",
            cycle, communication_count, domain->rank, local_changed, changed, continue_local, continue_global);
    fclose(file);
}

static void dump_boundary_state(MPIDomain* domain, EikonalGrid* g_processus, int cycle, int communication_count, const int* start, const int* frontier, const int* changed_cells){
    char filename[256];
    snprintf(filename, sizeof(filename), "mpi_boundary_rank%d.txt", domain->rank);

    FILE* file = fopen(filename, "a");
    if (!file)
        return;

    fprintf(file, "=== cycle %d comm %d rank %d ===\n", cycle, communication_count, domain->rank);
    fprintf(file, "owned=[%d,%d] overlap=[%d,%d] n_overlap=%d top_ghost=%d bottom_ghost=%d\n",
            domain->i_owned_start, domain->i_owned_end,
            domain->i_start_overlap, domain->i_end_overlap,
            domain->n_overlap, domain->top_ghost, domain->bottom_ghost);

    int m = domain->m;

    fprintf(file, "top boundary:\n");
    for (int i = 0; i < domain->top_ghost; i++){
        int local_row = i;
        int global_row = domain->i_start_overlap + local_row;
        for (int j = 0; j < m; j++){
            int idx = local_row * m + j;
            fprintf(file, "  local=(%d,%d) global=(%d,%d) T=%.6f start=%d frontier=%d changed=%d\n",
                    local_row, j, global_row, j, g_processus->T[idx], start[idx], frontier[idx], changed_cells[idx]);
        }
    }

    fprintf(file, "bottom boundary:\n");
    for (int i = 0; i < domain->bottom_ghost; i++){
        int local_row = g_processus->n - 1 - i;
        int global_row = domain->i_start_overlap + local_row;
        for (int j = 0; j < m; j++){
            int idx = local_row * m + j;
            fprintf(file, "  local=(%d,%d) global=(%d,%d) T=%.6f start=%d frontier=%d changed=%d\n",
                    local_row, j, global_row, j, g_processus->T[idx], start[idx], frontier[idx], changed_cells[idx]);
        }
    }

    fclose(file);
}

static void dump_active_cells(MPIDomain* domain, EikonalGrid* g_processus, int cycle, int communication_count, const int* start, const int* frontier, const int* changed_cells){
    char filename[256];
    snprintf(filename, sizeof(filename), "mpi_active_rank%d.txt", domain->rank);

    FILE* file = fopen(filename, "a");
    if (!file)
        return;

    fprintf(file, "=== cycle %d comm %d rank %d ===\n", cycle, communication_count, domain->rank);

    int count = 0;
    for (int k = 0; k < g_processus->n * g_processus->m; k++){
        if (start[k] || frontier[k] || changed_cells[k]){
            int local_row = k / domain->m;
            int local_col = k % domain->m;
            int global_row = domain->i_start_overlap + local_row;
            int global_col = local_col;
            if (count < 160){
                fprintf(file, "  cell[%d] local=(%d,%d) global=(%d,%d) start=%d frontier=%d changed=%d T=%.6f\n",
                        k, local_row, local_col, global_row, global_col, start[k], frontier[k], changed_cells[k], g_processus->T[k]);
            }
            count++;
        }
    }

    fprintf(file, "active_count=%d\n", count);
    fclose(file);
}

static void dump_partial_grid(MPIDomain* domain, EikonalGrid* g_processus, int cycle, int communication_count){
    char filename[256];
    snprintf(filename, sizeof(filename), "mpi_grid_rank%d.txt", domain->rank);

    FILE* file = fopen(filename, "a");
    if (!file)
        return;

    fprintf(file, "=== cycle %d comm %d rank %d ===\n", cycle, communication_count, domain->rank);
    for (int i = 0; i < g_processus->n; i++){
        fprintf(file, "row %d: ", i);
        for (int j = 0; j < g_processus->m; j++){
            fprintf(file, "%.6f ", g_processus->T[i * g_processus->m + j]);
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
    fclose(file);
}

static void dump_halo_exchange(MPIDomain* domain, EikonalGrid* g_processus, int cycle, int communication_count, const int* changed_cells){
    char filename[256];
    snprintf(filename, sizeof(filename), "mpi_halo_rank%d.txt", domain->rank);

    FILE* file = fopen(filename, "a");
    if (!file)
        return;

    fprintf(file, "=== cycle %d comm %d rank %d ===\n", cycle, communication_count, domain->rank);
    fprintf(file, "neighbors: up=%d down=%d\n", domain->up_rank, domain->down_rank);

    int m = domain->m;
    fprintf(file, "top halo values:\n");
    for (int i = 0; i < domain->top_ghost; i++){
        int idx = i * m;
        fprintf(file, "  [%d] T=%.6f changed=%d\n", idx, g_processus->T[idx], changed_cells[idx]);
    }

    fprintf(file, "bottom halo values:\n");
    for (int i = 0; i < domain->bottom_ghost; i++){
        int idx = (g_processus->n - 1 - i) * m;
        fprintf(file, "  [%d] T=%.6f changed=%d\n", idx, g_processus->T[idx], changed_cells[idx]);
    }

    fclose(file);
}

void fim_solve_mpi(MPIDomain* domain, EikonalGrid* g_processus, Config2* cfg_processus, int* start, double epsilon, int nb_cycles){
    int n = g_processus->n;
    int m = g_processus->m;
    int ncell = n * m;

    // TEST
    /*printf("========== DEBUT fim_solve_mpi ==========\n");
    printf("[Processus %d] n=%d, m=%d, ncell=%d\n", domain->rank, n, m, ncell);
    int nb_start_init = 0;
    for (int k = 0; k < ncell; k++) {
        if (start[k]) nb_start_init++;
    }
    printf("[Processus %d] start contient %d cellules marquées au début\n", domain->rank, nb_start_init);*/
 
    int* frontier = (int*)malloc(ncell * sizeof(int));
    int* depth = (int*)malloc(ncell * sizeof(int));
    int* changed_cells = (int*)calloc(domain->n_overlap * domain->m, sizeof(int));

    int cycle = 0;
    int communication_count = 0;
    while(true){

        // TEST
        //printf("[Processus %d] Cycle %d - avant local_propagate\n", domain->rank, cycle);

        // On fait la propagation sur 'overlap' cellules de distance
        int local_changed = local_propagate(g_processus, start, domain->overlap, epsilon, depth, frontier);

        // On initialise continue_local avant de l'utiliser dans write_cycle_summary
        int continue_local = local_changed;

        if (cycle == 0 || cycle == 1){
            printf("[rank %d] cycle %d - before communication (local_changed=%d)\n", domain->rank, cycle, local_changed);
            fflush(stdout);
            write_cycle_summary(domain, cycle, communication_count, local_changed, 0, continue_local, 0);
            dump_boundary_state(domain, g_processus, cycle, communication_count, start, frontier, changed_cells);
            dump_active_cells(domain, g_processus, cycle, communication_count, start, frontier, changed_cells);
        }

        // Communication entre les processus
        int changed = exchange_overlap(domain, g_processus, changed_cells, cycle, communication_count);

        memset(start, 0, ncell * sizeof(int));

        communication_count++;
        if (cycle == 0 || cycle == 1){
            printf("[rank %d] cycle %d - after communication %d (changed=%d)\n", domain->rank, cycle, communication_count, changed);
            fflush(stdout);
            write_cycle_summary(domain, cycle, communication_count, local_changed, changed, continue_local, 0);
            dump_boundary_state(domain, g_processus, cycle, communication_count, start, frontier, changed_cells);
            dump_active_cells(domain, g_processus, cycle, communication_count, start, frontier, changed_cells);
            dump_halo_exchange(domain, g_processus, cycle, communication_count, changed_cells);
            dump_partial_grid(domain, g_processus, cycle, communication_count);
        }

        // Les mailles de départ du prochain cycle sont celles sur lesquelles on s'est arrêté au cycle précédent et les mailles qui ont été modifiées pendant la communication
        //memset(start, 0, ncell * sizeof(int));
        continue_local = local_changed || changed;
        for (int k = 0; k < ncell; k++){
            if (frontier[k] || changed_cells[k]){
                start[k] = 1;
            }
        }

        cycle++;

        // Dans le cas où on teste sans communications, on continue tant qu'il y a du travail
        #if TEST_MODE
            if (!continue_local){
                // TEST
                printf("[Processus %d] Cycle %d - plus de travail, arrêt\n", domain->rank, cycle);

                break;
            }
            if (nb_cycles > 0 && cycle >= nb_cycles){
                // TEST
                printf("[Processus %d] Cycle %d - NB_CYCLES ATTEINT, arrêt\n", domain->rank, cycle);

                break;
            }
            continue;
        #endif

        if (nb_cycles > 0 && cycle >= nb_cycles)
            break;

        // TEST
        //printf("[Processus %d] Arrivé au cycle %d, continue_local = %d\n", domain->rank, cycle, continue_local);
        //fflush(stdout);

        // Si tous les processus n'ont plus de travail alors on arrête tout
        int continue_global = 0;
        MPI_Allreduce(&continue_local, &continue_global, 1, MPI_INT, MPI_MAX, domain->comm);

        if (cycle == 0 || cycle == 1){
            printf("[rank %d] cycle %d - before global decision: local=%d global=%d\n", domain->rank, cycle, continue_local, continue_global);
            fflush(stdout);
            write_debug_trace(domain, cycle, "after_comm", local_changed, changed, continue_local, continue_global, communication_count);
        }

        // TEST
        //printf("[Processus %d] Sortie de Allreduce, continue_global = %d\n", domain->rank, cycle, continue_global);
        //fflush(stdout);
        if (!continue_global){
            if (cycle == 0 || cycle == 1){
                printf("[rank %d] cycle %d - stop because continue_global=0\n", domain->rank, cycle);
                fflush(stdout);
                write_debug_trace(domain, cycle, "stop_global", local_changed, changed, continue_local, continue_global, communication_count);
            }
            break;
        }
/*
        if (communication_count >= DEBUG_COMM_PHASES){
            printf("[rank %d] cycle %d - stop after communication phase %d\n", domain->rank, cycle, communication_count);
            fflush(stdout);
            write_debug_trace(domain, cycle, "stop_debug_phase", local_changed, changed, continue_local, continue_global, communication_count);
            break;
        }*/

        if (cycle >= 5000) {
            printf("[Processus %d] Arrêt forcé de test à 5000 cycles.\n", domain->rank);
            break;
        }
    }

    //TEST
    //printf("[Processus %d] FIN fim_solve_mpi après %d cycles\n", domain->rank, cycle);
    //printf("========== FIN fim_solve_mpi ==========\n");

    free(frontier);
    free(depth);
    free(changed_cells);
}