#ifndef SOLVEEIKONAL2D_H
#define SOLVEEIKONAL2D_H

#include <float.h>

#define EIKONAL_INF DBL_MAX // Valeur infinie pour les cellules qui n'ont pas encore été atteintes
#define EIKONAL_EPS 1e-12


typedef struct{
    int n,m;    // Dimensions de la grille
    double h;   // Pas de la grille (égal pour x et y)
    double* F;  // Vitesse de propagation
    double* T;  // Temps d'arrivée
}EikonalGrid;


// Alloue et initialise une grille
EikonalGrid* eikonal_grid_create(int n, int m, double h);

// Nettoie la mémoire de la grille
void eikonal_grid_free(EikonalGrid* g);

// Définit la vitesse F (constante) pour l'ensemble de la grille
void eikonal_grid_set_speed_constant(EikonalGrid* g, double F);

// Définit la valeur F cellule par cellule à l'aide d'un tableau
void eikonal_grid_set_speed(EikonalGrid* g, int i, int j, double F);

// Obtenir la valeur de la vitesse cellule par cellule
double eikonal_grid_get_speed(const EikonalGrid* g, int i, int j);

// Marque une cellule comme un obstacle
void eikonal_grid_set_obstacle(EikonalGrid* g, int i, int j);

// Vérifie si une cellule est dans un mur
int eikonal_grid_is_obstacle(const EikonalGrid* g, int i, int j);

// Solution locale de l'équation eikonale en 2D pour la cellule (i,j)
// Renvoie la nouvelle valeur estimée de T 
double eikonal_solve_local(const EikonalGrid* g, int i, int j);

// Enregistre la matrice T dans un fichier texte (compatible avec Numpy pour la visualisation)
int eikonal_save_matrix(const EikonalGrid* g, const char* filename);

// Enregistre les valeurs de F dans un fichier texte
int save_speed(const EikonalGrid* g, const char* filename);

// Enregistre les tags associées aux sources les plus proches
int eikonal_save_tags(const EikonalGrid* g, const int* source_tag, const char* filename);



// Accès inline T(i,j)
static inline double eikonal_T(const EikonalGrid* g, int i, int j){
    if (i < 0 || j < 0 || i >= g->n || j >= g->m)
        return EIKONAL_INF;

    return g->T[i * g->m + j];
}


#endif /* SOLVEEIKONAL2D_H */