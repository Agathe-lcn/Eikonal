#include "../include/SolveEikonal2D.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

EikonalGrid* eikonal_grid_create(int n, int m, double h){
    EikonalGrid* g = (EikonalGrid*)malloc(sizeof(EikonalGrid));
    if(!g)
        return NULL;

    g->n = n;
    g->m = m;
    g->h = h;

    g->F = (double *)malloc(n * m * sizeof(double));
    g->T = (double *)malloc(n * m * sizeof(double));

    if (!g->F || !g->T){
        free(g->F);
        free(g->T);
        free(g);
        return NULL;
    }

    // Default initialisation: speed = 1 and infinite time
    for (int k=0; k < n*m; k++){
        g->F[k] = 0;
        g->T[k] = EIKONAL_INF;
    }

    return g;
}


void eikonal_grid_free(EikonalGrid* g){
    if (!g)
        return;

    free(g->F);
    free(g->T);
    free(g);
}


void eikonal_grid_set_speed_constant(EikonalGrid* g, double speed){
    for (int k=0; k < g->n * g->m; k++)
        g->F[k] = speed;
}


void eikonal_grid_set_speed(EikonalGrid* g, const double* F){
    memcpy(g->F, F, g->n * g->m * sizeof(double));
}


void eikonal_grid_set_obstacle(EikonalGrid* g, int i, int j){
    g->F[i * g->m + j] = 0.0;
    g->T[i * g->m + j] = EIKONAL_INF;
}


double eikonal_solve_local(const EikonalGrid* g, int i, int j){
    double F_ij = g->F[i * g->m + j];
    
    // If obstacle or zero speed -> T is infinite
    if (F_ij <= EIKONAL_EPS)
        return EIKONAL_INF;

    double h = g->h;

    double Tx = DBL_MAX;
    double Ty = DBL_MAX;

    double T;

    // Direction x
    if (i>0){
        T = g->T[(i-1) * g->m + j];
        if (T < Tx)
            Tx = T;
    }
    if (i < g->n - 1){
        T = g->T[(i+1) * g->m + j];
        if (T < Tx)
            Tx = T;
    }

    // Direction y
    if (j>0){
        T = g->T[i * g->m + j-1];
        if (T<Ty)
            Ty = T;
    }
    if (j < g->m - 1){
        T = g->T[i * g->m + j+1];
        if (T<Ty)
            Ty = T;
    }

    // Case where no neighbor is valid
    if (Tx == DBL_MAX && Ty == DBL_MAX)
        return EIKONAL_INF;

    
    double T_new;
    double hF = h/F_ij;

    if (Tx == DBL_MAX)
        // Only one neighbor available in y
        T_new = Ty + hF;
    else if (Ty == DBL_MAX)
        // Only one neighbor available in x
        T_new = Tx + hF;
    else{
        double hF2 = hF * hF;
        double diff = Tx - Ty;
        double disc = 2.0 * hF2 - diff * diff;

        double T_tempo = 0.5 * (Tx + Ty + sqrt(disc));
        if (disc >= 0 && T_tempo >= fmax(Tx, Ty))
            // Verified causal relationship
            //T_new = 0.5 * (Tx + Ty + sqrt(disc));
            T_new = T_tempo;
        else{
            if (Tx < Ty)
                T_new = Tx + hF;
            else
                T_new = Ty + hF;
        }
    }

    return T_new;
}


int eikonal_save_matrix(const EikonalGrid *g, const char *filename){
    FILE *f = fopen(filename, "w");
    if (!f)
        return -1;

    for (int i=0; i < g->n; i++){
        for (int j=0; j < g->m; j++){
            double v = g->T[i * g->m + j];
            if (v >= EIKONAL_INF)
                fprintf(f, "inf ");
            else
                fprintf(f, "%.10g",v);
            if (j < g->m - 1)
                fprintf(f, " ");
        }
        fprintf(f, "\n");
    }
    fclose(f);
    return 0;
}
