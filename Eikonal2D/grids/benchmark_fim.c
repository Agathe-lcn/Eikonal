#include "../include/config.h"

#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

// Fonction pour créer un fichier de configuration avec une source au centre
int create_config_file(const char* filename, int n, int m) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        printf("Erreur: Impossible de créer %s\n", filename);
        return 0;
    }

    double h = 1.0 / n;
    
    // Écrire la configuration
    fprintf(file, "n = %d\n", n);
    fprintf(file, "m = %d\n", m);
    fprintf(file, "h = %.10f\n", h);
    fprintf(file, "\n");
    fprintf(file, "sources:\n");
    
    // Source au centre
    int centre_i = n / 2;
    int centre_j = m / 2;
    fprintf(file, "%d %d\n", centre_i, centre_j);
    
    fprintf(file, "\n");
    fprintf(file, "walls:\n");
    // Aucun mur
    
    fclose(file);
    return 1;
}

// Execute une FIM et mesure le temps
double run_fim(const char* config_file){
    // Lecture de la configuration
    Config cfg = read_config(config_file);

    if (cfg.n == 0 || cfg.m == 0){
        printf("Erreur: configuration invalide (n ou m pas défini)\n");
        free_config(&cfg);
        return 1.0;
    }

    if (cfg.h <= 0){
        printf("Erreur: configuration invalide (h pas défini)\n");
        free_config(&cfg);
        return 1.0;
    }

    if (cfg.nsources == 0){
        printf("Erreur: aucune source spécifiée\n");
        free_config(&cfg);
        return 1.0;
    }

    // Création de la grille
    EikonalGrid* g = eikonal_grid_create(cfg.n, cfg.m, cfg.h);
    if (!g){
        printf("Erreur: impossible de créer la grille\n");
        free_config(&cfg);
        return 1.0;
    }

    // Vitesse constante égale à 1
    eikonal_grid_set_speed_constant(g, 1.0);

    // Définir une vitesse qui varie au lieu d'une vitesse constante
    //set_variable_speed(g, cfg.n, cfg.m, cfg.h);

    // On ajoute les murs
    add_walls(g, cfg);

    // On supprime les sources dans le mur
    removing_sources_in_walls(g, &cfg);

    if (cfg.nsources == 0){
        printf("Erreur: toutes les sources sont dans un mur, il n'y a rien à propager.\n");
        eikonal_grid_free(g);
        free_config(&cfg);
        return 1.0;
    }

    // Stockage des informations sur les sources
    save_sources(cfg);

    save_speed(g, "speed.txt");

    clock_t start = clock();
    fim_solve(g, cfg.src_i, cfg.src_j, cfg.nsources, EPSILON, cfg.max_depth);
    clock_t end = clock();

    double cpu_time = ((double)(end - start)) / CLOCKS_PER_SEC;

    // Nettoyage
    eikonal_grid_free(g);
    free_config(&cfg);

    return cpu_time;
}

// Exécute le benchmark pour une config donnée
void run_benchmark(int n, int m, int max_depth, int num_runs){
    char tempo_config[256];
    snprintf(tempo_config, sizeof(tempo_config), "tempo_config_%dx%d.txt", n, m);

    // Création du fichier de configuration temporaire
    if (!create_config_file(tempo_config, n, m)){
        printf("  Erreur: Impossible de créer le fichier de configuration\n");
        return;
    }

    double* times = (double*)malloc(num_runs * sizeof(double));
    int valid_runs = 0;

    // On exécute plusieurs fois la FIM pour pouvoir ensuite faire une moyenne des résultats
    for (int run = 0; run < num_runs; run++){
        double time = run_fim(tempo_config);

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
            fprintf(results, "%d\t%d\t%d\t%.6f\t%d\n", n, m, max_depth, avg, valid_runs);
            fclose(results);
        }
    }
    else
        printf("Erreur: Aucun run valide\n");

    free(times);
    remove(tempo_config);
}



int main(int argc, char** argv){
    // Paramètres du benchmark
    int grid_sizes[] = {50, 200, 500, 1000, 2500, 5000};
    int len_grid_sizes = sizeof(grid_sizes) / sizeof(grid_sizes[0]);
    int num_runs = 10;
    double max_depth[] = {-1.0, 30.0};

    // Initialiser le fichier de résultats
    FILE* results = fopen("results_fim.txt", "w");
    if(!results){
        printf("Erreur: Impossible de créer results_fim.txt\n");
        return 1;
    }
    fprintf(results, "n\tm\tmax_depth\tavg_time_seconds\tnum_runs\n");
    fclose(results);

    // On lance le benchmark
    for (int k = 0; k<2; k++){
        for (int s = 0; s<len_grid_sizes; s++)
            run_benchmark(grid_sizes[s], grid_sizes[s], max_depth[k], num_runs);
    }

    return 0;
}