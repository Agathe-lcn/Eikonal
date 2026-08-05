#include "../include/FIM2D_mpi.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define TAG 0
#define TEST_MODE 0     // 1 pour désactiver les communications et 0 pour les activer

static int is_owned_local_row(const MPIDomain* domain, int i_local){
    return i_local >= domain->top_ghost && i_local < domain->top_ghost + domain->n_owned;
}

static void report_topology_error(int rank, const char* message){
    if (rank == 0){
        fprintf(stderr, "%s\n", message);
        fflush(stderr);
    }
}

static int source_depth_is_better(int new_depth, int old_depth){
    return new_depth >= 0 && (old_depth < 0 || new_depth < old_depth);
}

typedef struct {
    int relative_index;
    double value;
    int depth;
} OptimHaloCell;

static int optimized_stop_check_period(const MPIDomain* domain){
    if (!domain->optim_com_mpi)
        return 1;
    if (domain->overlap < 1)
        return 1;
    return domain->overlap;
}

static int halo_value_is_send_worthy(double current_t, int current_depth, double cached_t, int cached_depth){
    if (current_depth < 0)
        return 0;
    if (current_t >= EIKONAL_INF)
        return 0;
    if (cached_t >= EIKONAL_INF)
        return 1;
    if (current_t < cached_t - EIKONAL_EPS)
        return 1;
    if (fabs(current_t - cached_t) <= EIKONAL_EPS && source_depth_is_better(current_depth, cached_depth))
        return 1;
    return 0;
}

static int init_last_sent_cache(double** t_cache, int** depth_cache, int band_size){
    *t_cache = (double*)malloc((size_t)band_size * sizeof(double));
    *depth_cache = (int*)malloc((size_t)band_size * sizeof(int));
    if (!*t_cache || !*depth_cache){
        free(*t_cache);
        free(*depth_cache);
        *t_cache = NULL;
        *depth_cache = NULL;
        return 0;
    }

    for (int k = 0; k < band_size; k++){
        (*t_cache)[k] = EIKONAL_INF;
        (*depth_cache)[k] = -1;
    }

    return 1;
}

static void record_halo_sendrecv(MPIDomain* domain, unsigned long long send_bytes, unsigned long long recv_bytes){
    domain->halo_sendrecv_calls++;
    if (send_bytes > 0){
        domain->halo_messages_sent++;
        domain->halo_bytes_sent += send_bytes;
    }
    if (recv_bytes > 0){
        domain->halo_messages_received++;
        domain->halo_bytes_received += recv_bytes;
    }
}

static int build_sparse_halo_payload(const MPIDomain* domain, const EikonalGrid* g_processus, const int* source_depth,
    int local_row_start, int row_count, double* cached_t, int* cached_depth,
    OptimHaloCell* out_cells){
    int band_size = row_count * domain->m;
    int send_count = 0;

    for (int rel = 0; rel < band_size; rel++){
        int local_index = local_row_start * domain->m + rel;
        double current_t = g_processus->T[local_index];
        int current_depth = source_depth[local_index];

        if (!halo_value_is_send_worthy(current_t, current_depth, cached_t[rel], cached_depth[rel]))
            continue;

        out_cells[send_count].relative_index = rel;
        out_cells[send_count].value = current_t;
        out_cells[send_count].depth = current_depth;
        cached_t[rel] = current_t;
        cached_depth[rel] = current_depth;
        send_count++;
    }

    return send_count;
}

static int apply_sparse_halo_payload(MPIDomain* domain, EikonalGrid* g_processus, int* changed_cells, int* source_depth,
    int base_offset, int row_count, const OptimHaloCell* cells, int recv_count){
    int any_change = 0;
    int band_size = row_count * domain->m;

    for (int entry = 0; entry < recv_count; entry++){
        int rel = cells[entry].relative_index;
        int target;
        double old_t;
        double new_t;
        int new_depth;
        int depth_better;

        if (rel < 0 || rel >= band_size)
            continue;

        target = base_offset + rel;
        old_t = g_processus->T[target];
        new_t = cells[entry].value;
        new_depth = cells[entry].depth;
        depth_better = source_depth_is_better(new_depth, source_depth[target]);

        if (new_t < old_t - EIKONAL_EPS){
            g_processus->T[target] = new_t;
            source_depth[target] = new_depth;
            changed_cells[target] = 1;
            any_change = 1;
            domain->halo_cells_updated++;
            continue;
        }

        if (fabs(new_t - old_t) <= EIKONAL_EPS && depth_better){
            source_depth[target] = new_depth;
            changed_cells[target] = 1;
            any_change = 1;
            domain->halo_cells_updated++;
        }
    }

    return any_change;
}

static int activate_owned_from_changed_cell(const MPIDomain* domain, int* start, int* source_depth, int max_depth, int index){
    int i = index / domain->m;
    int j = index % domain->m;
    int source_cell_depth = source_depth[index];

    if (is_owned_local_row(domain, i)){
        if (max_depth >= 0 && (source_cell_depth < 0 || source_cell_depth > max_depth))
            return 0;
        start[index] = 1;
        return 1;
    }

    if (i < domain->top_ghost){
        int target = domain->top_ghost * domain->m + j;
        int target_depth = source_cell_depth + (domain->top_ghost - i);
        if (source_cell_depth < 0)
            return 0;
        if (max_depth >= 0 && target_depth > max_depth)
            return 0;
        if (source_depth[target] < 0 || source_depth[target] > target_depth)
            source_depth[target] = target_depth;
        start[target] = 1;
        return 1;
    }

    if (i >= domain->top_ghost + domain->n_owned){
        int last_owned = domain->top_ghost + domain->n_owned - 1;
        int target = last_owned * domain->m + j;
        int target_depth = source_cell_depth + (i - last_owned);
        if (source_cell_depth < 0)
            return 0;
        if (max_depth >= 0 && target_depth > max_depth)
            return 0;
        if (source_depth[target] < 0 || source_depth[target] > target_depth)
            source_depth[target] = target_depth;
        start[target] = 1;
        return 1;
    }

    return 0;
}

MPIDomain* topology_create(int n, int m, double h, int overlap){
    // Topologie cartésienne 1D (sous forme de bandes)
    int ndims=1;
    int nproc_per_dim[1] = {0};
    int periods_per_dim[1] = {0};
    int rank, nproc;

    MPI_Comm_rank(MPI_COMM_WORLD, &(rank));
    MPI_Comm_size(MPI_COMM_WORLD, &(nproc));

    if (nproc < 1){
        if (rank == 0){
            fprintf(stderr, "Erreur: nombre de processus MPI invalide (%d).\n", nproc);
            fflush(stderr);
        }
        return NULL;
    }
    if (n < 1){
        if (rank == 0){
            fprintf(stderr, "Erreur: n=%d invalide, le nombre de lignes doit etre >= 1.\n", n);
            fflush(stderr);
        }
        return NULL;
    }
    if (m < 1){
        if (rank == 0){
            fprintf(stderr, "Erreur: m=%d invalide, le nombre de colonnes doit etre >= 1.\n", m);
            fflush(stderr);
        }
        return NULL;
    }
    if (h <= 0.0){
        if (rank == 0){
            fprintf(stderr, "Erreur: h=%g invalide, le pas de grille doit etre > 0.\n", h);
            fflush(stderr);
        }
        return NULL;
    }
    if (nproc > n){
        if (rank == 0){
            fprintf(stderr, "Erreur: nombre de processus MPI trop grand (%d) pour n=%d, certains rangs auraient 0 ligne owned.\n", nproc, n);
            fflush(stderr);
        }
        return NULL;
    }

    int min_owned = n / nproc;
    if (overlap < 1){
        report_topology_error(rank, "Erreur: overlap doit etre un entier >= 1.");
        return NULL;
    }
    if (overlap > min_owned){
        if (rank == 0){
            fprintf(stderr, "Erreur: overlap=%d invalide, il doit etre <= n_owned minimal=%d pour ce decoupage MPI.\n", overlap, min_owned);
            fflush(stderr);
        }
        return NULL;
    }

    MPIDomain* domain = (MPIDomain*)malloc(sizeof(MPIDomain));
    if (!domain)
        return NULL;

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
    int Q = min_owned;
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

    if (n_owned < 1){
        if (rank == 0){
            fprintf(stderr, "Erreur: decoupage MPI invalide, le rang %d aurait %d ligne owned.\n", rank, n_owned);
            fflush(stderr);
        }
        free(domain);
        return NULL;
    }
    if (i_start < 0 || i_start + n_owned > n){
        if (rank == 0){
            fprintf(stderr, "Erreur: decoupage MPI incoherent pour le rang %d (i_start=%d, n_owned=%d, n=%d).\n", rank, i_start, n_owned, n);
            fflush(stderr);
        }
        free(domain);
        return NULL;
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
    domain->max_depth = -1;
    domain->optim_com_mpi = 0;
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
    domain->solver_cycles = 0;
    domain->halo_exchange_rounds = 0;
    domain->halo_sendrecv_calls = 0;
    domain->halo_messages_sent = 0;
    domain->halo_messages_received = 0;
    domain->halo_bytes_sent = 0;
    domain->halo_bytes_received = 0;
    domain->halo_cells_updated = 0;
    domain->allreduce_calls = 0;
    domain->allreduce_payload_bytes = 0;
    domain->allreduce_skipped_cycles = 0;
    domain->last_sent_up_t = NULL;
    domain->last_sent_down_t = NULL;
    domain->last_sent_up_depth = NULL;
    domain->last_sent_down_depth = NULL;

    if (domain->up_rank != MPI_PROC_NULL && domain->top_ghost > 0){
        if (!init_last_sent_cache(&domain->last_sent_up_t, &domain->last_sent_up_depth, overlap * m)){
            if (rank == 0)
                fprintf(stderr, "Erreur: impossible d'allouer le cache d'envoi halo vers le haut.\n");
            MPI_Comm_free(&domain->comm);
            MPI_Comm_free(&domain->exch_comm);
            free(domain);
            return NULL;
        }
    }

    if (domain->down_rank != MPI_PROC_NULL && domain->bottom_ghost > 0){
        if (!init_last_sent_cache(&domain->last_sent_down_t, &domain->last_sent_down_depth, overlap * m)){
            if (rank == 0)
                fprintf(stderr, "Erreur: impossible d'allouer le cache d'envoi halo vers le bas.\n");
            free(domain->last_sent_up_t);
            free(domain->last_sent_up_depth);
            MPI_Comm_free(&domain->comm);
            MPI_Comm_free(&domain->exch_comm);
            free(domain);
            return NULL;
        }
    }

    if (domain->i_start_overlap < 0 || domain->i_end_overlap >= n || domain->n_overlap < domain->n_owned){
        if (rank == 0){
            fprintf(stderr, "Erreur: recouvrement MPI incoherent (rank=%d, owned=[%d,%d], overlap=[%d,%d], n_overlap=%d, n=%d).\n",
                rank,
                domain->i_owned_start,
                domain->i_owned_end,
                domain->i_start_overlap,
                domain->i_end_overlap,
                domain->n_overlap,
                n);
            fflush(stderr);
        }
        free(domain);
        return NULL;
    }

    // Test
    //printf("[rank %d] top_ghost=%d bottom_ghost=%d i_owned=[%d,%d] i_overlap=[%d,%d] n_overlap=%d\n", rank, domain->top_ghost, domain->bottom_ghost, domain->i_owned_start, domain->i_owned_end, domain->i_start_overlap, domain->i_end_overlap, domain->n_overlap);
    //fflush(stdout);

    return domain;
}


void topology_free(MPIDomain* domain){
    if (!domain)
        return;

    free(domain->last_sent_up_t);
    free(domain->last_sent_down_t);
    free(domain->last_sent_up_depth);
    free(domain->last_sent_down_depth);
    MPI_Comm_free(&domain->comm);
    MPI_Comm_free(&domain->exch_comm);
    free(domain);
}


static int exchange_overlap_full(MPIDomain* domain, EikonalGrid* g_processus, int* changed_cells, int* source_depth){
    // Si on veut tester sans les communications, alors on ne fait rien
    #if TEST_MODE
        return 0;
    #endif

    int m = domain->m;
    int overlap = domain->overlap;
    int any_change = 0;

    domain->halo_exchange_rounds++;

    // TEST
    //printf("[rank %d] exchange_overlap : m=%d overlap=%d buffer_size=%d n_overlap=%d\n", domain->rank, m, overlap, overlap * m, domain->n_overlap);
    //fflush(stdout);

    // On remet tout le tableau à zéro
    memset(changed_cells, 0, domain->n_overlap*m*sizeof(int));

    MPI_Request reqs[4];

    double* send_up = NULL;
    double* send_down = NULL;
    double* recv_up = NULL;
    double* recv_down = NULL;
    int* send_up_depth = NULL;
    int* send_down_depth = NULL;
    int* recv_up_depth = NULL;
    int* recv_down_depth = NULL;


    // Pour l'échange avec le voisin du dessus, on vérifie que le processus n'est pas tout en haut
    int do_up = domain->up_rank != MPI_PROC_NULL && domain->top_ghost > 0;
    if (do_up){
        send_up = (double*)malloc(overlap * m * sizeof(double));
        recv_up = (double*)malloc(overlap * m * sizeof(double));
        send_up_depth = (int*)malloc(overlap * m * sizeof(int));
        recv_up_depth = (int*)malloc(overlap * m * sizeof(int));
        memcpy(send_up, g_processus->T + domain->top_ghost*m, overlap * m * sizeof(double));
        memcpy(send_up_depth, source_depth + domain->top_ghost*m, overlap * m * sizeof(int));
        domain->halo_sendrecv_calls += 2;
        domain->halo_messages_sent += 2;
        domain->halo_messages_received += 2;
        domain->halo_bytes_sent += (unsigned long long)(overlap * m * (sizeof(double) + sizeof(int)));
        domain->halo_bytes_received += (unsigned long long)(overlap * m * (sizeof(double) + sizeof(int)));

        // TEST
        //printf("[rank %d] Sendrecv avec up_rank=%d (%d doubles)\n", domain->rank, domain->up_rank, overlap*m);
        //fflush(stdout);

        MPI_Sendrecv(send_up, overlap*m, MPI_DOUBLE, domain->up_rank, TAG, recv_up, overlap*m, MPI_DOUBLE, domain->up_rank, TAG, domain->exch_comm, MPI_STATUS_IGNORE);
        MPI_Sendrecv(send_up_depth, overlap*m, MPI_INT, domain->up_rank, TAG + 1, recv_up_depth, overlap*m, MPI_INT, domain->up_rank, TAG + 1, domain->exch_comm, MPI_STATUS_IGNORE);
    }


    // Pour l'échange avec le voisin du dessous, on vérifie que le processus n'est pas tout en bas
    int do_down = domain->down_rank != MPI_PROC_NULL && domain->bottom_ghost > 0;
    if (do_down){
        send_down = (double*)malloc(overlap * m * sizeof(double));
        recv_down = (double*)malloc(overlap * m * sizeof(double));
        send_down_depth = (int*)malloc(overlap * m * sizeof(int));
        recv_down_depth = (int*)malloc(overlap * m * sizeof(int));
        int send_offset = (domain->n_overlap - domain->bottom_ghost - overlap) * m;
        memcpy(send_down, g_processus->T + send_offset, overlap * m * sizeof(double));
        memcpy(send_down_depth, source_depth + send_offset, overlap * m * sizeof(int));
        domain->halo_sendrecv_calls += 2;
        domain->halo_messages_sent += 2;
        domain->halo_messages_received += 2;
        domain->halo_bytes_sent += (unsigned long long)(overlap * m * (sizeof(double) + sizeof(int)));
        domain->halo_bytes_received += (unsigned long long)(overlap * m * (sizeof(double) + sizeof(int)));

        // TEST
        //printf("[rank %d] Sendrecv avec down_rank=%d (%d doubles)\n", domain->rank, domain->down_rank, overlap*m);
        //fflush(stdout);

        MPI_Sendrecv(send_down, overlap*m, MPI_DOUBLE, domain->down_rank, TAG, recv_down, overlap*m, MPI_DOUBLE, domain->down_rank, TAG, domain->exch_comm, MPI_STATUS_IGNORE);
        MPI_Sendrecv(send_down_depth, overlap*m, MPI_INT, domain->down_rank, TAG + 1, recv_down_depth, overlap*m, MPI_INT, domain->down_rank, TAG + 1, domain->exch_comm, MPI_STATUS_IGNORE);
    }


    // Test
    //printf("[rank %d] échanges Send/Recv termines\n", domain->rank);
    //fflush(stdout);



    // Comparaison des valeurs de T entre le processus courant et son voisin du haut pour garder le minimum à chaque cellule
    if (recv_up){
        for (int k = 0; k < overlap*m; k++){
            double old_T = g_processus->T[k];
            double new_T = recv_up[k];
            if (new_T < old_T - EIKONAL_EPS){
                g_processus->T[k] = new_T;
                source_depth[k] = recv_up_depth[k];
                any_change = 1;
                changed_cells[k] = 1;
                domain->halo_cells_updated++;

                // Test
                //printf("[rank %d] up: cellule k=%d mise à jour %.4f -> %.4f\n", domain->rank, k, old_T, new_T);
            }
        }

        free(send_up);
        free(recv_up);
        free(send_up_depth);
        free(recv_up_depth);
    }


    // Comparaison des valeurs de T entre le processus courant et son voisin du bas pour garder le minimum à chaque cellule
    if (recv_down){
        int beginning = (domain->n_overlap - overlap) * m;
        for (int k = beginning; k < beginning + overlap*m; k++){
            double old_T = g_processus->T[k];
            double new_T = recv_down[k - beginning];
            if (new_T < old_T - EIKONAL_EPS){
                g_processus->T[k] = new_T;
                source_depth[k] = recv_down_depth[k - beginning];
                any_change = 1;
                changed_cells[k] = 1;
                domain->halo_cells_updated++;

                // Test
                //printf("[rank %d] down: cellule k=%d mise à jour %.4f -> %.4f\n", domain->rank, k, old_T, new_T);
            }
        }

        free(send_down);
        free(recv_down);
        free(send_down_depth);
        free(recv_down_depth);
    }

    // Test
    //printf("[rank %d] exchange_overlap terminé, any_change=%d\n", domain->rank, any_change);
    //fflush(stdout);

    return any_change;
}

static int exchange_overlap_optimized(MPIDomain* domain, EikonalGrid* g_processus, int* changed_cells, int* source_depth){
    int any_change = 0;
    int m = domain->m;
    int overlap = domain->overlap;
    int band_size = overlap * m;

    memset(changed_cells, 0, (size_t)domain->n_overlap * m * sizeof(int));
    domain->halo_exchange_rounds++;

    if (domain->up_rank != MPI_PROC_NULL && domain->top_ghost > 0){
        int send_count = 0;
        int recv_count = 0;
        OptimHaloCell* send_cells = (OptimHaloCell*)malloc((size_t)band_size * sizeof(OptimHaloCell));

        if (!send_cells){
            free(send_cells);
            printf("Erreur [rang %d]: allocation impossible pour l'echange halo optimise vers le haut.\n", domain->rank);
            return 0;
        }

        send_count = build_sparse_halo_payload(domain, g_processus, source_depth,
            domain->top_ghost, overlap,
            domain->last_sent_up_t, domain->last_sent_up_depth,
            send_cells);

        MPI_Sendrecv(&send_count, 1, MPI_INT, domain->up_rank, TAG + 10,
            &recv_count, 1, MPI_INT, domain->up_rank, TAG + 10,
            domain->exch_comm, MPI_STATUS_IGNORE);
        record_halo_sendrecv(domain, sizeof(int), sizeof(int));

        if (send_count > 0 || recv_count > 0){
            int recv_alloc = recv_count > 0 ? recv_count : 1;
            OptimHaloCell* recv_cells = (OptimHaloCell*)malloc((size_t)recv_alloc * sizeof(OptimHaloCell));

            if (!recv_cells){
                free(send_cells);
                free(recv_cells);
                printf("Erreur [rang %d]: allocation impossible pour la reception halo optimisee depuis le haut.\n", domain->rank);
                return 0;
            }

            MPI_Sendrecv(send_cells, send_count * (int)sizeof(OptimHaloCell), MPI_BYTE, domain->up_rank, TAG + 11,
                recv_cells, recv_count * (int)sizeof(OptimHaloCell), MPI_BYTE, domain->up_rank, TAG + 11,
                domain->exch_comm, MPI_STATUS_IGNORE);
            record_halo_sendrecv(domain,
                (unsigned long long)send_count * sizeof(OptimHaloCell),
                (unsigned long long)recv_count * sizeof(OptimHaloCell));

            any_change |= apply_sparse_halo_payload(domain, g_processus, changed_cells, source_depth,
                0, overlap, recv_cells, recv_count);

            free(recv_cells);
        }

        free(send_cells);
    }

    if (domain->down_rank != MPI_PROC_NULL && domain->bottom_ghost > 0){
        int send_count = 0;
        int recv_count = 0;
        int send_offset = (domain->n_overlap - domain->bottom_ghost - overlap) * m;
        int recv_base = (domain->n_overlap - overlap) * m;
        int local_row_start = send_offset / m;
        OptimHaloCell* send_cells = (OptimHaloCell*)malloc((size_t)band_size * sizeof(OptimHaloCell));

        if (!send_cells){
            free(send_cells);
            printf("Erreur [rang %d]: allocation impossible pour l'echange halo optimise vers le bas.\n", domain->rank);
            return any_change;
        }

        send_count = build_sparse_halo_payload(domain, g_processus, source_depth,
            local_row_start, overlap,
            domain->last_sent_down_t, domain->last_sent_down_depth,
            send_cells);

        MPI_Sendrecv(&send_count, 1, MPI_INT, domain->down_rank, TAG + 10,
            &recv_count, 1, MPI_INT, domain->down_rank, TAG + 10,
            domain->exch_comm, MPI_STATUS_IGNORE);
        record_halo_sendrecv(domain, sizeof(int), sizeof(int));

        if (send_count > 0 || recv_count > 0){
            int recv_alloc = recv_count > 0 ? recv_count : 1;
            OptimHaloCell* recv_cells = (OptimHaloCell*)malloc((size_t)recv_alloc * sizeof(OptimHaloCell));

            if (!recv_cells){
                free(send_cells);
                free(recv_cells);
                printf("Erreur [rang %d]: allocation impossible pour la reception halo optimisee depuis le bas.\n", domain->rank);
                return any_change;
            }

            MPI_Sendrecv(send_cells, send_count * (int)sizeof(OptimHaloCell), MPI_BYTE, domain->down_rank, TAG + 11,
                recv_cells, recv_count * (int)sizeof(OptimHaloCell), MPI_BYTE, domain->down_rank, TAG + 11,
                domain->exch_comm, MPI_STATUS_IGNORE);
            record_halo_sendrecv(domain,
                (unsigned long long)send_count * sizeof(OptimHaloCell),
                (unsigned long long)recv_count * sizeof(OptimHaloCell));

            any_change |= apply_sparse_halo_payload(domain, g_processus, changed_cells, source_depth,
                recv_base, overlap, recv_cells, recv_count);

            free(recv_cells);
        }

        free(send_cells);
    }

    return any_change;
}

int exchange_overlap(MPIDomain* domain, EikonalGrid* g_processus, int* changed_cells, int* source_depth){
    if (domain->optim_com_mpi)
        return exchange_overlap_optimized(domain, g_processus, changed_cells, source_depth);

    return exchange_overlap_full(domain, g_processus, changed_cells, source_depth);
}


void local_propagate(MPIDomain* domain, EikonalGrid* g_processus, const int* start, int overlap, double epsilon, int* depth, int* frontier, int* source_depth, int max_depth){
    int n = g_processus->n;
    int m = g_processus->m;
    int ncell = n*m;

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
        int is_owned = is_owned_local_row(domain, i);
        int current_source_depth = source_depth[index];

        // TEST
        //printf("[Processus] Iteration %d: traitement de la cellule %d (%d,%d)\n", iterations, index, i, j);

        if (eikonal_grid_is_obstacle(g_processus, i, j))
            continue;
        if (max_depth >= 0 && (current_source_depth < 0 || current_source_depth > max_depth))
            continue;

        double T_old = g_processus->T[index];

        // Les lignes fantômes servent de conditions de bord importées: elles ne sont pas recalculées localement.
        if (is_owned && T_old != 0.0){
            g_processus->T[index] = eikonal_solve_local(g_processus, i, j);
        }

        double diff = is_owned ? fabs(g_processus->T[index] - T_old) : 0.0;

        //TEST
        //printf("[Processus] T_old=%f, T_new=%f, diff=%f\n", T_old, g_processus->T[index], diff);

        if (diff <= epsilon){   // Convergence
            // On parcourt les voisins qui sont dans la zone de propagation voulue
            int current_depth = depth[index];

            if (max_depth >= 0 && current_source_depth >= max_depth)
                continue;

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
                    if (!is_owned_local_row(domain, ni))
                        continue;

                    int next_source_depth = current_source_depth + 1;
                    if (max_depth >= 0 && next_source_depth > max_depth)
                        continue;

                    if (eikonal_grid_is_obstacle(g_processus, ni, nj))
                        continue;

                    int index_neighbor = ni * m + nj;
                    double T_neighbor_new = eikonal_solve_local(g_processus, ni, nj);
                    
                    if (T_neighbor_new < g_processus->T[index_neighbor] - epsilon){
                        g_processus->T[index_neighbor] = T_neighbor_new;
                        if (source_depth[index_neighbor] < 0 || source_depth[index_neighbor] > next_source_depth)
                            source_depth[index_neighbor] = next_source_depth;

                        // TEST
                        //printf("[Processus] Voisin (%d,%d) mis à jour: %f\n", ni, nj, T_neighbor_new);

                        // On met à jour la profondeur du voisin
                        if (depth[index_neighbor] < 0 || depth[index_neighbor] > current_depth +1)
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

    // Nettoyage
    list_free(narrow);
}


void fim_solve_mpi(MPIDomain* domain, EikonalGrid* g_processus, Config2* cfg_processus, int* start, double epsilon, int nb_cycles){
    int n = g_processus->n;
    int m = g_processus->m;
    int ncell = n * m;
    int max_depth = cfg_processus->max_depth;

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
    int* source_depth = (int*)malloc(ncell * sizeof(int));

    if (!frontier || !depth || !changed_cells || !source_depth){
        printf("Erreur [rang %d]: impossible d'allouer les tableaux de travail MPI.\n", domain->rank);
        free(frontier);
        free(depth);
        free(changed_cells);
        free(source_depth);
        return;
    }

    for (int k = 0; k < ncell; k++){
        source_depth[k] = -1;
        if (start[k] && g_processus->T[k] == 0.0)
            source_depth[k] = 0;
    }

    int cycle = 0;
    int stop_check_period = optimized_stop_check_period(domain);
    while(true){

        // TEST
        //printf("[Processus %d] Cycle %d - avant local_propagate\n", domain->rank, cycle);

        // On fait la propagation sur 'overlap' cellules de distance
        local_propagate(domain, g_processus, start, domain->overlap, epsilon, depth, frontier, source_depth, max_depth);

        // Communication entre les processus
        int changed = exchange_overlap(domain, g_processus, changed_cells, source_depth);
        // Les mailles de départ du prochain cycle sont celles sur lequelles on s'est arrêté au cycle précédent et les mailles qui ont été modifiées pendant la communication
        memset(start, 0, ncell * sizeof(int));
        int continue_local = 0;
        int frontier_count = 0;
        int changed_count = 0;
        for (int k = 0; k < ncell; k++){
            if (frontier[k]){
                start[k] = 1;
                continue_local = 1;
                frontier_count++;
            }
            if (changed_cells[k]){
                if (activate_owned_from_changed_cell(domain, start, source_depth, max_depth, k)){
                    continue_local = 1;
                    changed_count++;
                }
            }
        }

        cycle++;
        domain->solver_cycles = (unsigned long long)cycle;

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

        // En mode optimise on espace les Allreduce de critere d'arret avec une periode commune a tous les rangs.
        // Tous les processus restent ainsi alignes sur les memes cycles collectifs.
        if (domain->optim_com_mpi && (cycle % stop_check_period) != 0){
            domain->allreduce_skipped_cycles++;
            if (cycle >= 5000) {
                printf("[Processus %d] Arrêt forcé de test à 5000 cycles.\n", domain->rank);
                break;
            }
            continue;
        }

        // Si tous les processus n'ont plus de travail alors on arrête tout
        int continue_global = 0;
        domain->allreduce_calls++;
        domain->allreduce_payload_bytes += (unsigned long long)sizeof(int);
        MPI_Allreduce(&continue_local, &continue_global, 1, MPI_INT, MPI_MAX, domain->comm);

        // TEST
        //printf("[Processus %d] Sortie de Allreduce, continue_global = %d\n", domain->rank, cycle, continue_global);
        //fflush(stdout);

        if (!continue_global)
            break;

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
    free(source_depth);
}