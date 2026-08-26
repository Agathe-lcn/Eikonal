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
    int max_depth;
    int optim_com_mpi;

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

    // Statistiques de communication du solveur MPI
    unsigned long long solver_cycles;
    unsigned long long halo_exchange_rounds;
    unsigned long long halo_sendrecv_calls;
    unsigned long long halo_messages_sent;
    unsigned long long halo_messages_received;
    unsigned long long halo_bytes_sent;
    unsigned long long halo_bytes_received;
    unsigned long long halo_cells_updated;
    unsigned long long allreduce_calls;
    unsigned long long allreduce_payload_bytes;
    unsigned long long allreduce_skipped_cycles;

    // Profiling: temps cumulés (secondes) mesurés dans fim_solve_mpi
    double time_local_propagate;
    double time_exchange_overlap;
    double time_start_construction;
    double time_allreduce;
    double time_solve_total;

    // Cache du dernier halo envoye pour le mode OPTIM_COM_MPI.
    double* last_sent_up_t;
    double* last_sent_down_t;
    int* last_sent_up_depth;
    int* last_sent_down_depth;
} MPIDomain;

typedef struct{
    int n;
    int m;
    double h;
    int overlap;
    int max_depth;
    int optim_com_mpi;

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

typedef struct{
    int* flags;      // Tableau booleen pour deduplication
    int* indices;    // Liste compacte des indices marques a 1
    int count;       // Nombre d'entrees dans indices
    int capacity;    // Capacite allouee de indices
    int flag_size;   // Taille du tableau flags
} FlagList;

FlagList* flaglist_create(int flag_size, int capacity);
void flaglist_free(FlagList* fl);
void flaglist_set(FlagList* fl, int index);
void flaglist_clear(FlagList* fl);

// Crée la topologie 1D en bandes
MPIDomain* topology_create(int n, int m, double h, int overlap);

// Libère la topologie MPI
void topology_free(MPIDomain* domain);

// Echange la bande de recouvrement avec les voisins en haut et en bas, et applique le minimum sur T
// Retourne 1 si au moins une valeur a été améliorée pour le processus et 0 sinon
int exchange_overlap(MPIDomain* domain, EikonalGrid* g_processus, FlagList* changed_cells, int* source_depth);

// Propagation de l'onde à partir des mailles start sur autant de pixels que la valeur du recouvrement
// (Par exemple, si on a un recouvrement de 3 pixels alors chaque sous domain MPI propage l'onde sur 3 pixels)
// depth[k] permet de connaitre la distance entre la maille k et la source
// frontier[k] vaut 1 si la maille k a été parcourue lors du dernier tour de la FIM (elle deviendra donc une maille de départ lors du prochain appel à local_propagate), 0 sinon
void local_propagate(MPIDomain* domain, EikonalGrid* g_processus, FlagList* start, int overlap, double epsilon, int* depth, int* depth_gen, int cur_gen, FlagList* frontier, NodeList* narrow, int* source_depth, int max_depth);


// FIM avec utilisation du MPI
void fim_solve_mpi(MPIDomain* domain, EikonalGrid* g_processus, Config2* cfg, FlagList* start, double espilon, int nb_cycles);



int parse_named_int_exact(const char* line, const char* key, int* out_value);

int parse_named_bool_flag_exact(const char* line, const char* key, int* out_value);

Config2 read_config_mpi(const char* filename);

void free_config_mpi(Config2* cfg);

void add_walls_local(EikonalGrid* g, Config2 cfg, MPIDomain* domain);

int removing_sources_in_walls_local(EikonalGrid* g_processus, Config2* cfg, MPIDomain* domain);

void save_sources_mpi(Config2 cfg, int rank);

void save_local_result_mpi(const MPIDomain* domain, const EikonalGrid* g_processus);

void save_mpi_communication_report(const MPIDomain* domain);

void save_mpi_profiling_report(const MPIDomain* domain);

void initialize_grid_with_sources(EikonalGrid* g_processus, Config2* cfg_processus, MPIDomain* domain, FlagList* start);

#endif /* FIM2D_MPI_H */