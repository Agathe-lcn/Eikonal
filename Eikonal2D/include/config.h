#ifndef CONFIG_H
#define CONFIG_H

#include "FIM2D.h"

#include <stdbool.h>

#define EPSILON 1e-12
#define MAX_LINE 1024

// Structure de configuration
typedef struct{
    int n;
    int m;
    double h;
    int max_depth;

    int nsources;
    int* src_i;
    int* src_j;

    int nwalls;
    int* wall_c1;   // Colonne de départ
    int* wall_c2;   // Colonne de fin
    int* wall_r1;   // Ligne de départ
    int* wall_r2;   // Ligne de fin

    bool valid;
} Config;

// Lecture du fichier de configuration
Config read_config(const char* filename);

// Nettoyage de la configuration
void free_config(Config* cfg);

// Ajout des murs à la grille
void add_walls(EikonalGrid* g, Config cfg);

// Suppression des sources qui sont dans un mur
void removing_sources_in_walls(EikonalGrid* g, Config* cfg);

// Enregistrement des coordonnées des sources
void save_sources(Config cfg);

// Définit la vitesse
void set_variable_speed(EikonalGrid* g, int n, int m, double h);

#endif