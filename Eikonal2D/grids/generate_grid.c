#include "../include/config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char** argv){
    // Nom du fichier de configuration
    const char* config_file = "config.txt";
    if (argc  > 1)
        config_file = argv[1];

    // Lecture de la configuration
    Config cfg = read_config(config_file);

    if (!cfg.valid){
        printf("Erreur: Configuration invalide.\n");
        free_config(&cfg);
        return EXIT_FAILURE;
    }

    if (cfg.nsources == 0){
        printf("Erreur: aucune source spécifiée.\n");
        free_config(&cfg);
        return EXIT_FAILURE;
    }

    // Création de la grille
    EikonalGrid* g = eikonal_grid_create(cfg.n, cfg.m, cfg.h);
    if (!g){
        printf("Erreur: impossible de créer la grille.\n");
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

    // Exécution de la FIM
    fim_solve(g, cfg.src_i, cfg.src_j, cfg.nsources, EPSILON, -1.0, cfg.max_depth);

    eikonal_save_matrix(g, "matrix_fim.txt");

    // Nettoyage
    eikonal_grid_free(g);
    free_config(&cfg);
    
    return 0;
}