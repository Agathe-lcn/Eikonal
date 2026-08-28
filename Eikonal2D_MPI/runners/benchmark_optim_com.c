#include "../include/FIM2D_mpi.h"
#include "../include/FIM2D_io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <limits.h>

#define EPSILON 1e-12
#define MAX_LINE 1024

// Valeurs de OPTIM_COM_MPI à tester
static const int OPTIM_VALUES[] = {0, 1};
static const int NB_OPTIM_VALUES = (int)(sizeof(OPTIM_VALUES) / sizeof(OPTIM_VALUES[0]));

// Nombre de répétitions par valeur d'OPTIM_COM_MPI (pour moyenner le temps de calcul)
#define NB_RUNS 3

static long read_global_halo_exchange_rounds(const char* filename){
    FILE* f = fopen(filename, "r");
    if (!f){
        printf("Erreur: impossible d'ouvrir %s pour lire le nombre de communications\n", filename);
        return -1;
    }

    char line[MAX_LINE];
    long rounds = -1;

    while (fgets(line, sizeof(line), f)){
        if (strstr(line, "global_halo_exchange_rounds")){
            char* eq = strchr(line, '=');
            if (eq){
                rounds = strtol(eq + 1, NULL, 10);
            }
            break;
        }
    }

    fclose(f);
    return rounds;
}

static long long read_global_halo_payload_bytes_total(const char* filename){
    FILE* f = fopen(filename, "r");
    if (!f){
        printf("Erreur: impossible d'ouvrir %s pour lire le volume de données échangées\n", filename);
        return -1;
    }

    char line[MAX_LINE];
    long long bytes_total = -1;

    while (fgets(line, sizeof(line), f)){
        if (strstr(line, "global_halo_payload_bytes_total")){
            char* eq = strchr(line, '=');
            if (eq){
                bytes_total = strtoll(eq + 1, NULL, 10);
            }
            break;
        }
    }

    fclose(f);
    return bytes_total;
}

int main(int argc, char** argv){
    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    FIMIO_Init(MPI_COMM_WORLD);

    Config2 cfg_processus;

    if (rank == 0){
        const char* config_file = "config.txt";
        if (argc == 2){
            config_file = argv[1];
        }
        else if (argc > 2){
            printf("Erreur: arguments invalides. Usage: benchmark_optim_com [config.txt]\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        cfg_processus = read_config_mpi(config_file);

        if (!cfg_processus.valid){
            free_config_mpi(&cfg_processus);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (cfg_processus.n == 0 || cfg_processus.m == 0){
            printf("Erreur: configuration invalide (n ou m pas défini)\n");
            free_config_mpi(&cfg_processus);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (cfg_processus.h <= 0){
            printf("Erreur: configuration invalide (h pas défini)\n");
            free_config_mpi(&cfg_processus);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (cfg_processus.nsources == 0){
            printf("Erreur: aucune source spécifiée\n");
            free_config_mpi(&cfg_processus);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (cfg_processus.max_depth < -1){
            printf("Erreur: configuration invalide (max_depth doit etre >= 0, ou -1 pour desactiver la limitation)\n");
            free_config_mpi(&cfg_processus);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (cfg_processus.overlap <= 0){
            printf("Erreur: configuration invalide (overlap doit etre defini et > 0)\n");
            free_config_mpi(&cfg_processus);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    MPI_Bcast(&cfg_processus.n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cfg_processus.m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cfg_processus.h, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cfg_processus.max_depth, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cfg_processus.overlap, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cfg_processus.nsources, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cfg_processus.nwalls, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0){
        cfg_processus.src_i = (int*)malloc(cfg_processus.nsources * sizeof(int));
        cfg_processus.src_j = (int*)malloc(cfg_processus.nsources * sizeof(int));
        cfg_processus.wall_c1 = (int*)malloc(cfg_processus.nwalls * sizeof(int));
        cfg_processus.wall_c2 = (int*)malloc(cfg_processus.nwalls * sizeof(int));
        cfg_processus.wall_r1 = (int*)malloc(cfg_processus.nwalls * sizeof(int));
        cfg_processus.wall_r2 = (int*)malloc(cfg_processus.nwalls * sizeof(int));
        cfg_processus.src_i_overlap = (int*)malloc(cfg_processus.nsources * sizeof(int));
        cfg_processus.src_j_overlap = (int*)malloc(cfg_processus.nsources * sizeof(int));
    }

    MPI_Bcast(cfg_processus.src_i, cfg_processus.nsources, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cfg_processus.src_j, cfg_processus.nsources, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cfg_processus.wall_c1, cfg_processus.nwalls, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cfg_processus.wall_c2, cfg_processus.nwalls, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cfg_processus.wall_r1, cfg_processus.nwalls, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cfg_processus.wall_r2, cfg_processus.nwalls, MPI_INT, 0, MPI_COMM_WORLD);

    FILE* txt = NULL;
    if (rank == 0){
        txt = fopen("benchmark_optim_com_results.txt", "w");
        if (!txt){
            printf("Erreur: impossible de créer %s\n", "benchmark_optim_com_results.txt");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        // En-tête du fichier au format tabulaire
        fprintf(txt, "optim_com_mpi\trun\ttemps_sec\tglobal_halo_exchange_rounds\tglobal_halo_payload_bytes_total\n");
    }

    for (int idx = 0; idx < NB_OPTIM_VALUES; idx++){
        int optim_com_mpi = OPTIM_VALUES[idx];

        for (int run = 1; run <= NB_RUNS; run++){
            MPIDomain* domain = topology_create(cfg_processus.n, cfg_processus.m, cfg_processus.h, cfg_processus.overlap);
            if (!domain){
                if (rank == 0)
                    printf("Erreur: topology_create a échoué pour optim_com_mpi=%d\n", optim_com_mpi);
                free_config_mpi(&cfg_processus);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            domain->max_depth = cfg_processus.max_depth;
            domain->optim_com_mpi = optim_com_mpi;

            EikonalGrid* g_processus = eikonal_grid_create(domain->n_overlap, domain->m, domain->h);
            if (!g_processus){
                printf("Erreur [rang %d]: impossible de créer la grille (optim_com_mpi=%d)\n", rank, optim_com_mpi);
                free_config_mpi(&cfg_processus);
                topology_free(domain);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            eikonal_grid_set_speed_constant(g_processus, 1.0);

            add_walls_local(g_processus, cfg_processus, domain);

            int ns_local = removing_sources_in_walls_local(g_processus, &cfg_processus, domain);
            cfg_processus.nsources_overlap = ns_local;

            if (cfg_processus.nsources == 0){
                if (rank == 0)
                    printf("Erreur: toutes les sources sont dans un mur (optim_com_mpi=%d), rien à propager.\n", optim_com_mpi);
                eikonal_grid_free(g_processus);
                topology_free(domain);
                continue;
            }

            int ncell = domain->n_overlap * domain->m;
            int* start = (int*)calloc(ncell, sizeof(int));
            initialize_grid_with_sources(g_processus, &cfg_processus, domain, start);

            MPI_Barrier(MPI_COMM_WORLD);
            double t_start = MPI_Wtime();

            fim_solve_mpi(domain, g_processus, &cfg_processus, start, EPSILON, -1);

            MPI_Barrier(MPI_COMM_WORLD);
            double t_end = MPI_Wtime();
            double elapsed = t_end - t_start;

            save_mpi_communication_report(domain);
            MPI_Barrier(MPI_COMM_WORLD);

            long rounds = -1;
            long long bytes_total = -1;
            if (rank == 0){
                rounds = read_global_halo_exchange_rounds("mpi_communication_report.txt");
                bytes_total = read_global_halo_payload_bytes_total("mpi_communication_report.txt");
                fprintf(txt, "%d\t%d\t%.6f\t%ld\t%lld\n", optim_com_mpi, run, elapsed, rounds, bytes_total);
                fflush(txt);
            }

            free(start);
            eikonal_grid_free(g_processus);
            topology_free(domain);
        }
    }

    free_config_mpi(&cfg_processus);

    FIMIO_Finalize();
    MPI_Finalize();

    return 0;
}