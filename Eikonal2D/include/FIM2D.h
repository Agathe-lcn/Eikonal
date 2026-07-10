#ifndef FIM_H
#define FIM_H

#include "SolveEikonal2D.h"

// Résout l'équation eikonale sur la grille g à l'aide de la méthode FIM
// src_i est le tableau des indices i des sources
// src_j est le tableau des indices j des sources
// ns est le nombre de sources
// Si max_radius est inférieur ou égal à 0, on parcourt l'intégralité de la grille; sinon, on applique le seuil

void fim_solve(EikonalGrid* g, const int* src_i, const int* src_j, int ns, double epsilon, double max_radius);


#endif /* FIM_H */