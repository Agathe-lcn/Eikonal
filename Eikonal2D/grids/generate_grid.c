#include "../include/FIM2D.h"

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdio.h>

// Wall structure
typedef struct{
    int col_start;  // Start column
    int col_end;    // End column
    int row_start;  // Start row
    int row_end;    // End row
} Wall;

// Parameter structure
typedef struct{
    int n;
    int m;
    double length; 
    int nsources;
    Wall* walls;    // Wall tables
    int nwalls;     // Number of walls
    int source_type;    // 0: random, 1: central, 2: sinusoidal, 3: circle, 4: diagonal, 5: horizontal, 6: vertical, 7: square
} Config;

// Sources structure
typedef struct{
    int* i; // Row table
    int* j; // Column table
    int count;  // Number of sources
} Sources;


/*
// Add the walls to the grid
void add_walls_to_grid(EikonalGrid* g, Config cfg){
    if (cfg.nwalls == 0)
        return;

    for (int w=0; w < cfg.nwalls; w++){
        int c1 = cfg.walls[w].col_start;
        int c2 = cfg.walls[w].col_end;
        int r1 = cfg.walls[w].row_start;
        int r2 = cfg.walls[w].row_end;

        // Check the limits
        if (c1 < 0)
            c1 = 0;
        if (c2 >= cfg.m)
            c2 = cfg.m - 1;
        if (r1 < 0)
            r1 = 0;
        if (r2 >= cfg.n)
            r2 = cfg.n - 1;

        for (int i = r1; i <= r2; i++){
            for (int j = c1; j <= c2; j++)
                eikonal_grid_set_obstacle(g, i, j);
        }
    }
}*/

// Function to check whether a point is inside a wall
int is_in_wall(int i, int j, Config cfg){
    for (int w=0; w < cfg.nwalls; w++){
        if (i >= cfg.walls[w].row_start && i <= cfg.walls[w].row_end && j >= cfg.walls[w].col_start && j >= cfg.walls[w].col_end)
            return 1;
    }
    return 0;
}
/*
// Parse the -wall arguments
int parse_walls(int argc, char** argv, Config* cfg){
    cfg->walls = NULL;
    cfg->nwalls = 0;

    for (int i = 1; i<argc; i++){
        if (strcmp(argv[i], "-wall") == 0){
            if (i + 4 < argc){
                int c1 = atoi(argv[i+1]);
                int c2 = atoi(argv[i+2]);
                int r1 = atoi(argv[i+3]);
                int r2 = atoi(argv[i+4]);

                // Ensure that c1 <= c2 and r1 <= r2
                if (c1 > c2){
                    int c_temp = c2;
                    c2 = c1;
                    c1 = c_temp;
                }
                if (r1 > r2){
                    int r_temp = r2;
                    r2 = r1;
                    r1 = r_temp;
                }

                cfg->nwalls++;
                cfg->walls = (Wall*)realloc(cfg->walls, cfg->nwalls * sizeof(Wall));
                cfg->walls[cfg->nwalls-1].col_start = c1;
                cfg->walls[cfg->nwalls-1].col_end = c2;
                cfg->walls[cfg->nwalls-1].row_start = r1;
                cfg->walls[cfg->nwalls-1].row_end = r2;

                i += 4;
            }
            else
                printf("Error: -wall requires 4 arguments");
        }
    }
    return 0;
}*/

/*
// Find a free position (not in a wall)
int find_free_position(Config* cfg, int* i, int* j){
    for (int i = 0; i < cfg.n; i++){
        for (int j = 0; j < cfg.m; j++){
            if (!is_in_wall)
        }
    }
}*/


// Generating sources by type
void generate_sources(Config cfg, int** src_i, int** src_j, int* ns){
    int n = cfg.n;
    int m = cfg.m;
    int s = cfg.nsources;

    switch(cfg.source_type){
        // O: random
        case 0:{
            *ns = s;
            *src_i = (int*)malloc((*ns) * sizeof(int));
            *src_j = (int*)malloc((*ns) * sizeof(int));
            srand(time(NULL));

            for (int k=0; k < *ns; k++){
                int valid = 0;
                int attempts = 0;
                while (!valid && attempts < 1000){
                    int ti = rand() % n;
                    int tj = rand() % m;
                    if (!is_in_wall(ti, tj, cfg)){
                        (*src_i)[k] = ti;
                        (*src_j)[k] = tj;
                        valid = 1;
                    }
                    attempts++;
                }
            }
            break;
        }

        // 1: central
        case 1:{
            *ns = 1;
            *src_i = (int*)malloc(sizeof(int));
            *src_j = (int*)malloc(sizeof(int));
            (*src_i)[0] = n/2;
            (*src_j)[0] = m/2;

            // If the source is in a wall, we move it
            if (is_in_wall((*src_i)[0], (*src_j)[0], cfg)){
                for (int j = (*src_j)[0] + 1; j < m; j++){
                    if (!is_in_wall((*src_i)[0], j, cfg)){
                        (*src_j)[0] = j;
                        break;
                    }
                }

                if (is_in_wall((*src_i)[0], (*src_j)[0], cfg)){
                    for (int j = (*src_j)[0] - 1; j>=0; j--){
                        if (!is_in_wall((*src_i)[0], j, cfg)){
                            (*src_j)[0] = j;
                            break;
                        }
                    }
                }
            }
            break;
        }

        // 2: sinusoidal
        case 2:{
            *ns = s;
            *src_i = (int*)malloc((*ns) * sizeof(int));
            *src_j = (int*)malloc((*ns) * sizeof(int));

            double amplitude = (n-1) / 4.0;
            double periode = (double)m / s;

            for (int k = 0; k < *ns; k++){
                int j = (int)(k * periode + periode/2);
                if (j >= m)
                    j = m-1;

                int i = n/2 + (int)(amplitude * sin(2 * M_PI * k/ (double)*ns));
                if (i<0)
                    i = 0;
                if (i >= n)
                    i = n-1;

                if (is_in_wall(i, j, cfg)){
                    for (int di=-1; di <= 1; di++){
                        for (int dj = -1; dj <= 1; dj++){
                            int i_temp = i + di;
                            int j_temp = j + dj;
                            if (0 <= i_temp < n && 0 <= j_temp < m && !is_in_wall(i_temp, j_temp, cfg)){
                                i = i_temp;
                                j = j_temp;
                                break;
                            }
                        }
                    }
                }
                (*src_i)[k] = i;
                (*src_j)[k] = j;
            }
            break;
        }

        // 3: circle
        case 3:{
            *ns = s;
            *src_i = (int*)malloc((*ns) * sizeof(int));
            *src_j = (int*)malloc((*ns) * sizeof(int));

            int center_i = n / 2;
            int center_j = m / 2;
            int radius;
            if (n < m)
                radius = n / 4;
            else
                radius = m / 4;

            for (int k=0; k < *ns; k++){
                double theta = 2 * M_PI * k / (double)*ns;
                int i = center_i + (int)(radius * sin(theta));
                int j = center_j + (int)(radius * cos(theta));

                if (i<0)
                    i = 0;
                if (i >= n)
                    i = n-1;
                if (j<0)
                    j = 0;
                if (j >= m)
                    j = m-1;

                if (is_in_wall(i, j, cfg)){
                    for (int di = -1; di <= 1; di++){
                        for (int dj = -1; dj <= 1; dj++){
                            int i_temp = i + di;
                            int j_temp = j + dj;

                            if (0 <= i_temp < n && 0 <= j_temp < m && !is_in_wall(i_temp, j_temp, cfg)){
                                i = i_temp;
                                j = j_temp;
                                break;
                            }
                        }
                    }
                }
                (*src_i)[k] = i;
                (*src_j)[k] = j;
            }
            break;
        }

        // 4: diagonal
        case 4:{
            *ns = s;
            *src_i = (int*)malloc((*ns) * sizeof(int));
            *src_j = (int*)malloc((*ns) * sizeof(int));

            int step_i = n / (*ns - 1);
            int step_j = m / (*ns - 1);

            for (int k = 0; k < *ns; k++){
                int i = (k + 1) * step_i;
                int j = (k + 1) * step_j;
            }
        }
    }
}