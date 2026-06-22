#include "../include/SolveEikonal2D.h"
#include <stdio.h>
#include <math.h>

int main(){
    printf("\nTest SolveEikonal2D.c \n");

    // Create a 10x10 grid
    EikonalGrid* g = eikonal_grid_create(10, 10, 0.1);
    if (!g){
        printf("eikonal_grid_create failed\n");
        return 1;
    }
    printf("eikonal_grid_create OK\n");

    // Test set_speed_constant
    eikonal_grid_set_speed_constant(g, 1.0);
    if (g->F[0] == 1.0)
        printf("eikonal_grid_set_speed_constant OK\n");
    else
        printf("eikonal_grid_set_speed_constant failed\n");

    // Test solve_local on a cell
    g->T[0] = 0.0;  //source
    double T_new = eikonal_solve_local(g, 0, 1);
    double T_expected = 0.1;
    if (fabs(T_new - T_expected) < 1e-6)
        printf("eikonal_solve_local OK (T_new = %f, T_expected = %f)\n", T_new, T_expected);
    else 
        printf("eikonal_solve_local failed (T_new = %f, T_expected = %f)\n", T_new, T_expected);

    eikonal_grid_free(g);
    printf("eikonal_grid_free OK \n\n");
    return 0;
}