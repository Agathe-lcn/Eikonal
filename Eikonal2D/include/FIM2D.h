#ifndef FIM_H
#define FIM_H

#include "SolveEikonal2D.h"

// Solves the eikonal equation on the grid g using the FIM
// src_i is the array of sources indices i
// src_j is the array of sources indices j
// ns is the number of sources
// If max_radius is less than 0, then we traverse the entire grid; otherwise, we apply the threshold

void fim_solve(EikonalGrid* g, const int* src_i, const int* src_j, int ns, double epsilon, double max_radius);


#endif /* FIM_H */