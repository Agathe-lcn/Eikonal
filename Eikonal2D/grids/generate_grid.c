#include "../include/FIM2D.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EPSILON 1e-12
#define MAX_LINE 1024

typedef struct{
    int n;
    int m;
    double h;

    int nsources;
    int* src_i;
    int* src_j;

    int nwalls;
    int* wall_c1;   // Colonne de départ
    int* wall_c2;   // Colonne de fin
    int* wall_r1;   // Ligne de départ
    int* wall_r2;   // Ligne de fin
} Config;


// Lecture du fichier de configuration
Config read_config(const char* filename){
    Config cfg = {};
    FILE* file = fopen(filename, "r");
    if (!file){
        printf("Erreur: Impossible d'ouvrir %s\n", filename);
        return cfg;
    }

    char line[MAX_LINE];
    int section = 0;    // 0: aucune, 1: sources, 2: murs
    int max_sources = 200;
    int max_walls = 200;

    // Allocations
    cfg.src_i = (int*)malloc(max_sources * sizeof(int));
    cfg.src_j = (int*)malloc(max_sources * sizeof(int));
    cfg.wall_c1 = (int*)malloc(max_walls * sizeof(int));
    cfg.wall_c2 = (int*)malloc(max_walls * sizeof(int));
    cfg.wall_r1 = (int*)malloc(max_walls * sizeof(int));
    cfg.wall_r2 = (int*)malloc(max_walls * sizeof(int));

    while (fgets(line, MAX_LINE, file)){
        // On ignore les commentaires et les lignes vides
        if (line[0] == '#' || line[0] == '\n')
            continue;

        // Suppression des \n
        line[strcspn(line, "\n")] = '\0';

        // Détections des sections
        if (strstr(line, "sources:") != NULL){
            section = 1;
            continue;
        }

        if (strstr(line, "walls:") != NULL){
            section = 2;
            continue;
        }

        // Lecture des paramètres de la grille
        if (sscanf(line, "n = %d", &cfg.n) == 1)
            continue;
        if (sscanf(line, "m = %d", &cfg.m) == 1)
            continue;
        if (sscanf(line, "h = %lf", &cfg.h) == 1)
            continue;

        // Lecture des sources
        if (section == 1){
            int i, j;
            if (sscanf(line, "%d %d",&i, &j) == 2){
                cfg.src_i[cfg.nsources] = i;
                cfg.src_j[cfg.nsources] = j;
                cfg.nsources++;
            }
        }

        // Lecture des murs
        if (section == 2){
            int c1, c2, r1, r2;
            if (sscanf(line, "%d %d %d %d", &c1, &c2, &r1, &r2)){
                cfg.wall_c1[cfg.nwalls] = c1;
                cfg.wall_c2[cfg.nwalls] = c2;
                cfg.wall_r1[cfg.nwalls] = r1;
                cfg.wall_r2[cfg.nwalls] = r2;
                cfg.nwalls++;
            }
        }
    }

    fclose(file);
    return cfg;
}


// Nettoyage de la configuration
void free_config(Config* cfg){
    free(cfg->src_i);
    free(cfg->src_j);
    free(cfg->wall_c1);
    free(cfg->wall_c2);
    free(cfg->wall_r1);
    free(cfg->wall_r2);
}


// Ajout des murs à la grille
void add_walls(EikonalGrid* g, Config cfg){
    for (int w=0; w < cfg.nwalls; w++){
        int c1 = cfg.wall_c1[w];
        int c2 = cfg.wall_c2[w];
        int r1 = cfg.wall_r1[w];
        int r2 = cfg.wall_r2[w];

        // Vérification des limites
        if (c1 < 0)
            c1 = 0;
        if (c2 >= cfg.m)
            c2 = cfg.m - 1;
        if (r1 < 0)
            r1 = 0;
        if (r2 >= cfg.n)
            r2 = cfg.n - 1;

        for (int i=r1; i <= r2; i++){
            for (int j=c1; j <= c2; j++){
                eikonal_grid_set_obstacle(g, i, j);
            }
        }
    }
}

// Suppression des sources qui sont dans un mur
void removing_sources_in_walls(EikonalGrid* g, Config* cfg){
    int valid = 0;
    for (int s=0; s < cfg->nsources; s++){
        int i = cfg->src_i[s];
        int j = cfg->src_j[s];

        if (eikonal_grid_is_obstacle(g,i,j)){
            printf("Attention: La source %d située aux coordonnées (%d,%d) est dans un mur, elle va être ignorée.\n", s, i, j);
            continue;
        }

        // On garde seulement les sources valides
        cfg->src_i[valid] = i;
        cfg->src_j[valid] = j;
        valid++;
    }
    cfg->nsources = valid;
}

// Enregistrement des coordonnées des sources
void save_sources(Config cfg){
    FILE* file = fopen("coords_source.txt", "w");
    if (!file){
        printf("Erreur: Impossible de créer coords_source.txt\n");
        return;
    }

    for (int s=0; s < cfg.nsources; s++){
        double x = cfg.src_j[s] * cfg.h;
        double y = cfg.src_i[s] * cfg.h;
        fprintf(file, "%.6f %.6f\n", x, y);   
    }
    fclose(file);
}


void set_variable_speed(EikonalGrid *g, int n, int m, double h){
    for (int i = 0; i<n; i++){
        for (int j = 0; j<m; j++){
            double x = j*h;
            double y = i*h;

            double F = 1 - 0.8 * exp(- (pow(x - 0.65, 2) + pow(y - 0.65, 2)) / 0.02);

            eikonal_grid_set_speed(g, i, j, F);
        }
    }
}


int main(int argc, char** argv){
    // Nom du fichier de configuration
    const char* config_file = "config.txt";
    if (argc  > 1)
        config_file = argv[1];

    // Lecture de la configuration
    Config cfg = read_config(config_file);

    if (cfg.n == 0 || cfg.m == 0){
        printf("Erreur: configuration invalide (n ou m pas défini)\n");
        free_config(&cfg);
        return 1;
    }

    if (cfg.h <= 0){
        printf("Erreur: configuration invalide (h pas défini)\n");
        free_config(&cfg);
        return 1;
    }

    if (cfg.nsources == 0){
        printf("Erreur: aucune source spécifiée\n");
        free_config(&cfg);
        return 1;
    }

    // Création de la grille
    EikonalGrid* g = eikonal_grid_create(cfg.n, cfg.m, cfg.h);
    if (!g){
        printf("Erreur: impossible de créer la grille\n");
        free_config(&cfg);
        return 1;
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
        return 1;
    }

    // Stockage des informations sur les sources
    save_sources(cfg);

    save_speed(g, "speed.txt");

    // Exécutuin de la FIM
    fim_solve(g, cfg.src_i, cfg.src_j, cfg.nsources, EPSILON, -1.0);

    eikonal_save_matrix(g, "matrix_fim.txt");

    // Nettoyage
    eikonal_grid_free(g);
    free_config(&cfg);
    
    return 0;
}