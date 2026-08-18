#include "../include/config.h"

#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

// Execute une FIM et mesure le temps
double run_fim(const char* config_file, int max_depth){
    // Lecture de la configuration
    Config cfg = read_config(config_file);

    if (cfg.n == 0 || cfg.m == 0){
        printf("Erreur: configuration invalide (n ou m pas défini)\n");
        free_config(&cfg);
        return -1.0;
    }

    if (cfg.h <= 0){
        printf("Erreur: configuration invalide (h pas défini)\n");
        free_config(&cfg);
        return -1.0;
    }

    if (cfg.nsources == 0){
        printf("Erreur: aucune source spécifiée\n");
        free_config(&cfg);
        return -1.0;
    }

    // Création de la grille
    EikonalGrid* g = eikonal_grid_create(cfg.n, cfg.m, cfg.h);
    if (!g){
        printf("Erreur: impossible de créer la grille\n");
        free_config(&cfg);
        return -1.0;
    }

    // Vitesse constante égale à 1
    eikonal_grid_set_speed_constant(g, 1.0);

    // On ajoute les murs
    add_walls(g, cfg);

    // On supprime les sources dans le mur
    removing_sources_in_walls(g, &cfg);

    if (cfg.nsources == 0){
        printf("Erreur: toutes les sources sont dans un mur, il n'y a rien à propager.\n");
        eikonal_grid_free(g);
        free_config(&cfg);
        return -1.0;
    }

    // Stockage des informations sur les sources
    save_sources(cfg);

    save_speed(g, "speed.txt");

    clock_t start = clock();
    fim_solve(g, cfg.src_i, cfg.src_j, cfg.nsources, EPSILON, max_depth);
    clock_t end = clock();

    double cpu_time = ((double)(end - start)) / CLOCKS_PER_SEC;

    // Nettoyage
    eikonal_grid_free(g);
    free_config(&cfg);

    return cpu_time;
}

// Exécute le benchmark pour un fichier de configuration donné (avec sa taille) et un max_depth donné
void run_benchmark(const char* config_file, const char* dataset_name, int size, int max_depth, int num_runs){
    double* times = (double*)malloc(num_runs * sizeof(double));
    int valid_runs = 0;

    // On exécute plusieurs fois la FIM pour pouvoir ensuite faire une moyenne des résultats
    for (int run = 0; run < num_runs; run++){
        double time = run_fim(config_file, max_depth);

        if (time >= 0.0)
            times[valid_runs++] = time;
    }

    // On calcule la moyenne des résultats
    if (valid_runs > 0){
        double sum = 0.0;

        for (int i = 0; i < valid_runs; i++)
            sum += times[i];

        double avg = sum / valid_runs;

        // On enregistre le temps moyen dans un fichier
        FILE* results = fopen("results_fim.txt", "a");
        if (results){
            fprintf(results, "%s\t%d\t%d\t%.6f\t%d\n", dataset_name, size, max_depth, avg, valid_runs);
            fclose(results);
        }
    }
    else
        printf("Erreur: Aucun run valide pour %s\n", config_file);

    free(times);
}

int main(int argc, char** argv){
    // Paramètres du benchmark
    // Tailles disponibles dans configs/datasets/circle et configs/datasets/line
    int sizes[] = {100, 200, 500, 1000, 2000, 5000, 10000};
    int len_sizes = sizeof(sizes) / sizeof(sizes[0]);

    // Les deux scénarios (jeux de données) à comparer
    const char* datasets[] = {"circle", "line"};
    int len_datasets = sizeof(datasets) / sizeof(datasets[0]);

    int num_runs = 10;
    int max_depths[] = {-1, 10};
    int len_max_depths = sizeof(max_depths) / sizeof(max_depths[0]);

    const char* configs_dir = "../configs/datasets";

    // Initialiser le fichier de résultats
    FILE* results = fopen("results_fim.txt", "w");
    if(!results){
        printf("Erreur: Impossible de créer results_fim.txt\n");
        return 1;
    }
    fprintf(results, "dataset\tsize\tmax_depth\tavg_time_seconds\tnum_runs\n");
    fclose(results);

    // On lance le benchmark pour chaque dataset (circle / line), chaque taille et chaque max_depth
    for (int d = 0; d < len_datasets; d++){
        for (int s = 0; s < len_sizes; s++){
            char config_file[512];
            snprintf(config_file, sizeof(config_file), "%s/%s/%s_%d.txt", configs_dir, datasets[d], datasets[d], sizes[s]);

            for (int k = 0; k < len_max_depths; k++){
                printf("Benchmark: dataset=%s size=%d max_depth=%d\n", datasets[d], sizes[s], max_depths[k]);
                run_benchmark(config_file, datasets[d], sizes[s], max_depths[k], num_runs);
            }
        }
    }

    return 0;
}