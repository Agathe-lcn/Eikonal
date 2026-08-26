#include "../include/FIM2D_mpi.h"
#include "../include/FIM2D_io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <limits.h>

#define EPSILON 1e-5
#define TEST_MODE 0     // 1 pour désactiver les communications et 0 pour les activer

int main(int argc, char** argv){
    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    FIMIO_Init(MPI_COMM_WORLD);

    Config2 cfg_processus;
    int overlap = 0;

    if (rank == 0){

        // Nom du fichier de configuration
        const char* config_file = "config.txt";
        if (argc == 2){
            config_file = argv[1];
        }
        else if (argc > 2){
            printf("Erreur: arguments invalides. Usage: generate_grid_mpi [config.txt]\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Lecture de la configuration
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

        if (cfg_processus.overlap < 1){
            printf("Erreur: configuration invalide (overlap doit etre un entier >= 1)\n");
            free_config_mpi(&cfg_processus);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (cfg_processus.max_depth < -1){
            printf("Erreur: configuration invalide (max_depth doit etre >= 0, ou -1 pour desactiver la limitation)\n");
            free_config_mpi(&cfg_processus);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        overlap = cfg_processus.overlap;
    }

    // Distribution de overlap et de certains champs de cfg_processus à tous les processus
    MPI_Bcast(&overlap, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cfg_processus.n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cfg_processus.m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cfg_processus.h, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cfg_processus.max_depth, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cfg_processus.optim_com_mpi, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cfg_processus.nsources, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cfg_processus.nwalls, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Allocation et diffusion des tableaux de sources
    if (rank != 0){
        cfg_processus.src_i = (int*)malloc(cfg_processus.nsources * sizeof(int));
        cfg_processus.src_j = (int*)malloc(cfg_processus.nsources * sizeof(int));
        cfg_processus.wall_c1 = (int*)malloc(cfg_processus.nwalls * sizeof(int));
        cfg_processus.wall_c2 = (int*)malloc(cfg_processus.nwalls * sizeof(int));
        cfg_processus.wall_r1 = (int*)malloc(cfg_processus.nwalls * sizeof(int));
        cfg_processus.wall_r2 = (int*)malloc(cfg_processus.nwalls * sizeof(int));
        // Allocation des tableaux locaux
        cfg_processus.src_i_overlap = (int*)malloc(cfg_processus.nsources * sizeof(int));
        cfg_processus.src_j_overlap = (int*)malloc(cfg_processus.nsources * sizeof(int));
    }

    // Distribution des autres champs de cfg_processus (sauf les champs locaux)
    MPI_Bcast(cfg_processus.src_i, cfg_processus.nsources, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cfg_processus.src_j, cfg_processus.nsources, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cfg_processus.wall_c1, cfg_processus.nwalls, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cfg_processus.wall_c2, cfg_processus.nwalls, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cfg_processus.wall_r1, cfg_processus.nwalls, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cfg_processus.wall_r2, cfg_processus.nwalls, MPI_INT, 0, MPI_COMM_WORLD);

    // Définition de la topologie
    MPIDomain* domain = topology_create(cfg_processus.n, cfg_processus.m, cfg_processus.h, overlap);
    if (!domain){
        free_config_mpi(&cfg_processus);
        FIMIO_Finalize();
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    domain->max_depth = cfg_processus.max_depth;
    domain->optim_com_mpi = cfg_processus.optim_com_mpi;

    // Création de la grille
    EikonalGrid* g_processus = eikonal_grid_create(domain->n_overlap, domain->m, domain->h);
    if (!g_processus){
        printf("Erreur [rang %d]: impossible de créer la grille\n", rank);
        free_config_mpi(&cfg_processus);
        topology_free(domain);
        MPI_Finalize();
        return 1;
    }

    // Vitesse constante égale à 1
    eikonal_grid_set_speed_constant(g_processus, 1.0);

    // On ajoute les murs
    add_walls_local(g_processus, cfg_processus, domain);

    // On supprime les sources dans le mur
    int ns_local = removing_sources_in_walls_local(g_processus, &cfg_processus, domain);
    cfg_processus.nsources_overlap = ns_local;

    if (cfg_processus.nsources == 0){
        printf("Erreur: toutes les sources sont dans un mur, il n'y a rien à propager.\n");
        eikonal_grid_free(g_processus);
        free_config_mpi(&cfg_processus);
        topology_free(domain);
        MPI_Finalize();
        return 1;
    }

    // Initialisation des sources
    int ncell = domain-> n_overlap * domain->m;
    FlagList* start = flaglist_create(ncell, ncell);
    initialize_grid_with_sources(g_processus, &cfg_processus, domain, start);

    // Stockage des informations sur les sources
    save_sources_mpi(cfg_processus, rank);

    // Exécution de la FIM
    MPI_Barrier(MPI_COMM_WORLD);
    double t_fim_start = MPI_Wtime();
    fim_solve_mpi(domain, g_processus, &cfg_processus, start, EPSILON, -1);
    double t_fim_end = MPI_Wtime();

    double t_fim_local = t_fim_end - t_fim_start;
    double t_fim_max;
    MPI_Reduce(&t_fim_local, &t_fim_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0){
        printf("Temps FIM: %.6f secondes\n", t_fim_max);
    }

    // Profiling: sauvegarde du rapport de timing détaillé
    save_mpi_profiling_report(domain);

    // Sauvegarde du résultat local de chaque processus pour la visualisation par sous-domaine.
    double t_save_start = MPI_Wtime();
    save_local_result_mpi(domain, g_processus);
    save_mpi_communication_report(domain);
    double t_save_end = MPI_Wtime();

    // Dans le cas où on veut tester dans les communications, chaque processus sauvegarde ses résultats dans un fichier différent
    #if TEST_MODE
        char filename[256];
        snprintf(filename, sizeof(filename), "result_test_rank%d.txt", rank);

        FILE* file = fopen(filename, "w");
        if (file){
            int n = g_processus->n;
            int m = g_processus->m;
            for (int i = 0; i < n; i++){
                for (int j = 0; j < m; j++){
                    fprintf(file, "%.6f ", g_processus->T[i * m + j]);
                }
                fprintf(file, "\n");
            }
            fclose(file);
        }
        else {
            printf("[Processus %d] Erreur: impossible de sauvegarder %s\n", rank, filename);
        }
    #endif

    // Sauvegarde des résultats
    MPI_Info info = MPI_INFO_NULL;
    double t_chkpt_start = MPI_Wtime();
    int err = FIMIO_Checkpoint("matrix_fim.txt", domain, g_processus, -1, info);
    double t_chkpt_end = MPI_Wtime();
    if (err != MPI_SUCCESS && rank == 0)
        printf("Erreur lors de l'enregistrement des résultats de la FIM");

    // Affichage du profiling global (rang 0)
    double t_save_local = t_save_end - t_save_start;
    double t_chkpt_local = t_chkpt_end - t_chkpt_start;
    double t_save_max, t_chkpt_max;
    MPI_Reduce(&t_save_local, &t_save_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_chkpt_local, &t_chkpt_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0){
        printf("Temps sauvegarde (resultats + rapport com): %.6f secondes\n", t_save_max);
        printf("Temps checkpoint MPI-IO: %.6f secondes\n", t_chkpt_max);
        printf("Temps total (FIM + sauvegarde + checkpoint): %.6f secondes\n",
               t_fim_max + t_save_max + t_chkpt_max);
    }

    // Dans le cas où on veut tester dans les communications, chaque processus sauvegarde ses résultats dans un fichier différent
    #if TEST_MODE
        char filename[256];
        snprintf(filename, sizeof(filename), "result_test_rank%d.txt", rank);

        FILE* file = fopen(filename, "w");
        if (file){
            int n = g_processus->n;
            int m = g_processus->m;
            for (int i = 0; i < n; i++){
                for (int j = 0; j < m; j++){
                    fprintf(file, "%.6f ", g_processus->T[i * m + j]);
                }
                fprintf(file, "\n");
            }
            fclose(file);
        }
        else {
            printf("[Processus %d] Erreur: impossible de sauvegarder %s\n", rank, filename);
        }
    #endif

    // Nettoyage
    flaglist_free(start);
    eikonal_grid_free(g_processus);
    free_config_mpi(&cfg_processus);

    FIMIO_Finalize();

    MPI_Finalize();

    return 0;
}
