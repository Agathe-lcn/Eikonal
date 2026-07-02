#ifndef SOLVEEIKONAL2D_H
#define SOLVEEIKONAL2D_H

#include <float.h>

#define EIKONAL_INF DBL_MAX // Infinite value for cells that have not been reached
#define EIKONAL_EPS 1e-12


typedef struct{
    int n,m;    // Grid dimensions
    double h;   // Grid step (for both x and y)
    double* F;  // Propagation speed
    double* T;  // Time-of-arrival
}EikonalGrid;


// Allocates and initializes a grid
EikonalGrid* eikonal_grid_create(int n, int m, double h);

// Clears the grid's memory
void eikonal_grid_free(EikonalGrid* g);

// Sets the speed F for the entire grid (constant)
void eikonal_grid_set_speed_constant(EikonalGrid* g, double F);

// Sets the F-value cell by cell using an array
void eikonal_grid_set_speed(EikonalGrid* g, int i, int j, double F);

// Gets the F-value cell by cell
double eikonal_grid_get_speed(const EikonalGrid* g, int i, int j);

// Mark a cell as an obstacle
void eikonal_grid_set_obstacle(EikonalGrid* g, int i, int j);

// Check if a source is inside a wall
int eikonal_grid_is_obstacle(const EikonalGrid* g, int i, int j);

// Local solution of the 2D eikonal equation for cell (i,j)
// Returns the new estimated T value (first-order upwind)
double eikonal_solve_local(const EikonalGrid* g, int i, int j);

// Saves the T matrix to a text file (Numpy-compatible for visualization)
int eikonal_save_matrix(const EikonalGrid* g, const char* filename);

// Saves the F values to a text file
int save_speed(const EikonalGrid* g, const char* filename);

// Saves the tags with the closest sources
int eikonal_save_tags(const EikonalGrid* g, const int* source_tag, const char* filename);



// Utilities

// Inline access T(i,j)
static inline double eikonal_T(const EikonalGrid* g, int i, int j){
    if (i < 0 || j < 0 || i >= g->n || j >= g->m)
        return EIKONAL_INF;

    return g->T[i * g->m + j];
}


#endif /* SOLVEEIKONAL2D_H */