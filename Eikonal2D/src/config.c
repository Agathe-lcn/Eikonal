#include "../include/config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <stdbool.h>
#include <ctype.h>

// Vérifie si une chaîne est un entier
bool is_integer(const char* str){
    if (str == NULL || *str == '\0')
        return false;

    int i = 0;

    if (str[0] == '-')
        i = 1;

    // On vérifie que tous les caractères sont des chiffres
    while(str[i] != '\0'){
        if (!isdigit(str[i]))
            return false;
        i++;
    }

    // On vérifie qu'il y a au moins 1 chiffre
    if (i == 0 || (i == 1 && str[0] == '-'))
        return false;

    return true;
}

// Vérifie si une chaîne est un double
bool is_double(const char* str){
    if (str == NULL || *str == '\0')
        return false;

    int i = 0;
    int nb_point = 0;

    if (str[0] == '-')
        i = 1;

    while(str[i] != '\0'){
        if (str[i] == '.'){
            nb_point++;

            if (nb_point > 1)
                return false;
        }
        else if( !isdigit(str[i]))
            return false;

        i++;
    }

    // On vérifie qu'il y a au moins 1 chiffre
    if (i == 0 || (i == 1 && str[0] == '-'))
        return false;

    return true;
}

// Validation du fichier de configuration
bool validate_config(Config* cfg, const char* n_str, const char* m_str, const char* h_str, const char* max_depth_str){
    // n doit être un entier
    if (!is_integer(n_str)){
        printf("Erreur: n doit être un entier.\n");
        return false;
    }

    // m doit être un entier
    if (!is_integer(m_str)){
        printf("Erreur: m doit être un entier.\n");
        return false;
    }

    // h doit être un double
    if (!is_double(h_str)){
        printf("Erreur: h doit être un double.\n");
        return false;
    }

    // max_depth doit être un entier (s'il a été précisé dans le fichier).
    // S'il est absent (chaîne vide), on garde la valeur par défaut (-1) déjà fixée dans read_config.
    if (max_depth_str[0] != '\0' && !is_integer(max_depth_str)){
        printf("Erreur: max_depth doit être un entier.\n");
        return false;
    }


    // n doit être positif
    if (cfg->n <= 0){
        printf("Erreur: n doit être positif.\n");
        return false;
    }

    // m doit être positif
    if (cfg->m <= 0){
        printf("Erreur: m doit être positif.\n");
        return false;
    }

    // h doit être positif
    if (cfg->h <= 0.0){
        printf("Erreur: h doit être positif.\n");
        return false;
    }

    // max_depth doit être >= 0 quand il est activé, sinon -1 pour désactiver
    if (cfg->max_depth < -1 || (cfg->max_depth > -1 && cfg->max_depth < 0)){
        printf("Erreur: max_depth doit être un entier >= 0, ou -1 pour désactiver la limitation.\n");
        return false;
    }




    // Validation des sources
    for (int s=0; s < cfg->nsources; s++){
        if (cfg->src_i[s] < 0 || cfg->src_i[s] >= cfg->n){
            printf("Erreur: Coordonnée i de la source %d hors limites.\n",s);
            return false;
        }

        if (cfg->src_j[s] < 0 || cfg->src_j[s] >= cfg->m){
            printf("Erreur: Coordonnée j de la source %d hors limites.\n",s);
            return false;
        }
    }



    // Validation des murs
    for (int w=0; w < cfg->nwalls; w++){
        if(cfg->wall_c1[w] < 0 || cfg->wall_c1[w] > cfg->m){
            printf("Erreur: Coordonnée c1 du mur %d hors limites.\n",w);
            return false;
        }

        if (cfg->wall_c2[w] < 0 || cfg->wall_c2[w] > cfg->m){
            printf("Erreur: Coordonnée c2 du mur %d hors limites.\n",w);
            return false;
        }

        if (cfg->wall_r1[w] < 0 || cfg->wall_r1[w] > cfg->n){
            printf("Erreur: Coordonnée r1 du mur %d hors limites.\n",w);
            return false;
        }

        if (cfg->wall_r2[w] < 0 || cfg->wall_r2[w] > cfg->n){
            printf("Erreur: Coordonnée r2 du mur %d hors limites.\n",w);
            return false;
        }

        if (cfg->wall_c1[w] > cfg->wall_c2[w]){
            printf("Erreur: c1 > c2 pour le mur %d.\n",w);
            return false;
        }

        if (cfg->wall_r1[w] > cfg->wall_r2[w]){
            printf("Erreur: r1 > r2 pour le mur %d.\n", w);
            return false;
        }
    }

    return true;
}



// Lecture du fichier de configuration
Config read_config(const char* filename){
    Config cfg = {};
    cfg.max_depth = -1;
    FILE* file = fopen(filename, "r");
    if (!file){
        printf("Erreur: Impossible d'ouvrir %s\n", filename);
        return cfg;
    }

    char line[MAX_LINE];
    int section = 0;    // 0: aucune, 1: sources, 2: murs

    int max_sources = 10000;
    int max_walls = 1000;

    // Allocations initiales
    cfg.src_i = (int*)malloc(max_sources * sizeof(int));
    cfg.src_j = (int*)malloc(max_sources * sizeof(int));
    cfg.wall_c1 = (int*)malloc(max_walls * sizeof(int));
    cfg.wall_c2 = (int*)malloc(max_walls * sizeof(int));
    cfg.wall_r1 = (int*)malloc(max_walls * sizeof(int));
    cfg.wall_r2 = (int*)malloc(max_walls * sizeof(int));

    // Vérification des allocations
    if (!cfg.src_i || !cfg.src_j || !cfg.wall_c1 || !cfg.wall_c2 || !cfg.wall_r1 || !cfg.wall_r2){
        printf("Erreur: Echec d'allocation mémoire dans read_config.\n");
        free_config(&cfg);
        fclose(file);
        return cfg;
    }

    // Variables pour stocker les chaînes originales
    char n_str[MAX_LINE] = "";
    char m_str[MAX_LINE] = "";
    char h_str[MAX_LINE] = "";
    char max_depth_str[MAX_LINE] = "";

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

        // Lecture des paramètres n,m et h avec stockage des chaînes
        char tempo[MAX_LINE];

        // Lecture de n
        if (sscanf(line, "n = %s", tempo) == 1){
            strcpy(n_str,tempo);
            cfg.n = atoi(tempo);
            continue;
        }

        // Lecture de m
        if (sscanf(line, "m = %s", tempo) == 1){
            strcpy(m_str, tempo);
            cfg.m = atoi(tempo);
            continue;
        }

        // Lecture de h
        if (sscanf(line, "h = %s", tempo) == 1){
            strcpy(h_str, tempo);
            cfg.h = atof(tempo);
            continue;
        }

        // Lecture de max_depth
        if (sscanf(line, "max_depth = %s", tempo) == 1){
            strcpy(max_depth_str, tempo);
            cfg.max_depth = atof(tempo);
            continue;
        }

        // Lecture des sources
        if (section == 1){
            char i_str[MAX_LINE];
            char j_str[MAX_LINE];
            char extra[MAX_LINE];
            
            // Vérifier qu'il y a exactement 2 entiers
            if (sscanf(line, "%s %s %s", i_str, j_str, extra) == 2){
                // Vérifier que les coordonnées sont des entiers
                if (!is_integer(i_str) || !is_integer(j_str)) {
                    printf("Erreur: Les coordonnées des sources doivent être des entiers: %s.\n", line);
                    free_config(&cfg);
                    fclose(file);
                    cfg.valid = false;
                    return cfg;
                }

                // Réallocation dynamique si le maximum est atteint
                if (cfg.nsources >= max_sources){
                    max_sources *= 2;
                    int* new_src_i = (int*)realloc(cfg.src_i, max_sources * sizeof(int));
                    int* new_src_j = (int*)realloc(cfg.src_j, max_sources * sizeof(int));
                    if (!new_src_i || !new_src_j){
                        printf("Erreur: Echec de realloc pour les sources.\n");
                        free(new_src_i ? NULL : cfg.src_i);
                        free(new_src_j ? NULL : cfg.src_j);
                        cfg.src_i = new_src_i ? new_src_i : cfg.src_i;
                        cfg.src_j = new_src_j ? new_src_j : cfg.src_j;
                        free_config(&cfg);
                        fclose(file);
                        cfg.valid = false;
                        return cfg;
                    }
                    cfg.src_i = new_src_i;
                    cfg.src_j = new_src_j;
                }

                cfg.src_i[cfg.nsources] = atoi(i_str);
                cfg.src_j[cfg.nsources] = atoi(j_str);
                cfg.nsources++;
            } else {
                printf("Erreur: Les lignes contenant les coordonnées d'une source doivent contenir 2 entiers: %s\n", line);
                free_config(&cfg);
                fclose(file);
                cfg.valid = false;
                return cfg;
            }
        }

        // Lecture des murs
        if (section == 2){
            char c1_str[MAX_LINE];
            char c2_str[MAX_LINE];
            char r1_str[MAX_LINE];
            char r2_str[MAX_LINE];
            char extra[MAX_LINE];
            
            // Vérifier qu'il y a exactement 4 valeurs
            if (sscanf(line, "%s %s %s %s %s", c1_str, c2_str, r1_str, r2_str, extra) == 4){
                // Vérifier que les coordonnées sont des entiers
                if (!is_integer(c1_str) || !is_integer(c2_str) || !is_integer(r1_str) || !is_integer(r2_str)) {
                    printf("Erreur: Les coordonnées des murs doivent être des entiers: %s.\n", line);
                    free_config(&cfg);
                    fclose(file);
                    cfg.valid = false;
                    return cfg;
                }

                // Réallocation dynamique si le maximum est atteint
                if (cfg.nwalls >= max_walls){
                    max_walls *= 2;
                    int* new_wall_c1 = (int*)realloc(cfg.wall_c1, max_walls * sizeof(int));
                    int* new_wall_c2 = (int*)realloc(cfg.wall_c2, max_walls * sizeof(int));
                    int* new_wall_r1 = (int*)realloc(cfg.wall_r1, max_walls * sizeof(int));
                    int* new_wall_r2 = (int*)realloc(cfg.wall_r2, max_walls * sizeof(int));
                    if (!new_wall_c1 || !new_wall_c2 || !new_wall_r1 || !new_wall_r2){
                        printf("Erreur: Echec de realloc pour les murs.\n");
                        cfg.wall_c1 = new_wall_c1 ? new_wall_c1 : cfg.wall_c1;
                        cfg.wall_c2 = new_wall_c2 ? new_wall_c2 : cfg.wall_c2;
                        cfg.wall_r1 = new_wall_r1 ? new_wall_r1 : cfg.wall_r1;
                        cfg.wall_r2 = new_wall_r2 ? new_wall_r2 : cfg.wall_r2;
                        free_config(&cfg);
                        fclose(file);
                        cfg.valid = false;
                        return cfg;
                    }
                    cfg.wall_c1 = new_wall_c1;
                    cfg.wall_c2 = new_wall_c2;
                    cfg.wall_r1 = new_wall_r1;
                    cfg.wall_r2 = new_wall_r2;
                }

                cfg.wall_c1[cfg.nwalls] = atoi(c1_str);
                cfg.wall_c2[cfg.nwalls] = atoi(c2_str);
                cfg.wall_r1[cfg.nwalls] = atoi(r1_str);
                cfg.wall_r2[cfg.nwalls] = atoi(r2_str);
                cfg.nwalls++;
            } else {
                printf("Erreur: Les lignes contenant les coordonnées d'un mur doivent contenir 4 entiers: %s\n", line);
                free_config(&cfg);
                fclose(file);
                cfg.valid = false;
                return cfg;
            }
        }
    }

    fclose(file);

    if (!validate_config(&cfg, n_str, m_str, h_str, max_depth_str)) {
        cfg.valid = false;
    } else {
        cfg.valid = true;
    }

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

        for (int i=c1; i <= c2; i++){
            for (int j=r1; j <= r2; j++){
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