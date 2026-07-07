#include "../include/FIM2D.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EPSILON 1e-12
#define MAX_LINE 1024

// Configuration structure
typedef struct{
    int n;
    int m;
    double h;

    int nsources;
    int* src_i;
    int* src_j;

    int nwalls;
    int* wall_c1;   // Start column
    int* wall_c2;   // End column
    int* wall_r1;   // Start row
    int* wall_r2;   // End row
} Config;


// Reading the configuration file
Config read_config(const char* filename){
    Config cfg = {};
    FILE* file = fopen(filename, "r");
    if (!file){
        printf("Error: Unable to open %s\n", filename);
        return cfg;
    }

    char line[MAX_LINE];
    int section = 0;    // 0: none, 1: sources, 2: walls
    int max_sources = 200;
    int max_walls = 200;

    // Allocations
    cfg.src_i = (int*)malloc(max_sources * sizeof(int));
    cfg.src_j = (int*)malloc(max_sources * sizeof(int));
    cfg.wall_c1 = (int*)malloc(max_walls * sizeof(int));
    cfg.wall_c2 = (int*)malloc(max_walls * sizeof(int));
    cfg.wall_r1 = (int*)malloc(max_walls * sizeof(int));
    cfg.wall_r2 = (int*)malloc(max_walls * sizeof(int));

    while (fgets(line, MAX_LINE, file)){
        // Ignore comments and empty lines
        if (line[0] == '#' || line[0] == '\n')
            continue;

        // Delete the \n
        line[strcspn(line, "\n")] = '\0';

        // Detect sections
        if (strstr(line, "sources:") != NULL){
            section = 1;
            continue;
        }

        if (strstr(line, "walls:") != NULL){
            section = 2;
            continue;
        }

        // Read the grid settings
        if (sscanf(line, "n = %d", &cfg.n) == 1)
            continue;
        if (sscanf(line, "m = %d", &cfg.m) == 1)
            continue;
        if (sscanf(line, "h = %lf", &cfg.h) == 1)
            continue;

        // Read the sources
        if (section == 1){
            int i, j;
            if (sscanf(line, "%d %d",&i, &j) == 2){
                cfg.src_i[cfg.nsources] = i;
                cfg.src_j[cfg.nsources] = j;
                cfg.nsources++;
            }
        }

        // Read the walls
        if (section == 2){
            int c1, c2, r1, r2;
            if (sscanf(line, "%d %d %d %d", &c1, &c2, &r1, &r2)){
                cfg.wall_c1[cfg.nwalls] = c1;
                cfg.wall_c2[cfg.nwalls] = c2;
                cfg.wall_r1[cfg.nwalls] = r1;
                cfg.wall_r2[cfg.nwalls] = r2;
                cfg.nwalls++;
            }
        }
    }

    fclose(file);
    return cfg;
}


// Cleaning the configuration
void free_config(Config* cfg){
    free(cfg->src_i);
    free(cfg->src_j);
    free(cfg->wall_c1);
    free(cfg->wall_c2);
    free(cfg->wall_r1);
    free(cfg->wall_r2);
}


// Adding walls to the grid
void add_walls(EikonalGrid* g, Config cfg){
    for (int w=0; w < cfg.nwalls; w++){
        int c1 = cfg.wall_c1[w];
        int c2 = cfg.wall_c2[w];
        int r1 = cfg.wall_r1[w];
        int r2 = cfg.wall_r2[w];

        // Check the limits
        if (c1 < 0)
            c1 = 0;
        if (c2 >= cfg.m)
            c2 = cfg.m - 1;
        if (r1 < 0)
            r1 = 0;
        if (r2 >= cfg.n)
            r2 = cfg.n - 1;

        for (int i=r1; i <= r2; i++){
            for (int j=c1; j <= c2; j++){
                eikonal_grid_set_obstacle(g, i, j);
            }
        }
    }
}

// Removing the sources that are inside a wall
void removing_sources_in_walls(EikonalGrid* g, Config* cfg){
    int valid = 0;
    for (int s=0; s < cfg->nsources; s++){
        int i = cfg->src_i[s];
        int j = cfg->src_j[s];

        if (eikonal_grid_is_obstacle(g,i,j)){
            printf("Warning: source %d at (%d,%d) is inside a wall, it will be ignored.\n", s, i, j);
            continue;
        }

        // We keep only the valid sources
        cfg->src_i[valid] = i;
        cfg->src_j[valid] = j;
        valid++;
    }
    cfg->nsources = valid;
}

// Saving source coordinates
void save_sources(Config cfg){
    FILE* file = fopen("coords_source.txt", "w");
    if (!file){
        printf("Error: Unable to create coords_source.txt\n");
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


int main(int argc, char** argv){
    // Configuration file name
    const char* config_file = "config.txt";
    if (argc  > 1)
        config_file = argv[1];

    // Reading the configuration
    Config cfg = read_config(config_file);

    if (cfg.n == 0 || cfg.m == 0){
        printf("Error: Invalid configuration (n or m not defined)\n");
        free_config(&cfg);
        return 1;
    }

    if (cfg.h <= 0){
        printf("Error: Invalid configuration (h not defined)\n");
        free_config(&cfg);
        return 1;
    }

    if (cfg.nsources == 0){
        printf("Error: No source specified\n");
        free_config(&cfg);
        return 1;
    }

    // Creating the grid
    EikonalGrid* g = eikonal_grid_create(cfg.n, cfg.m, cfg.h);
    if (!g){
        printf("Error: Unable to create grid\n");
        free_config(&cfg);
        return 1;
    }

    // Constant speed of 1
    eikonal_grid_set_speed_constant(g, 1.0);

    // Set variable speed instead of constant speed
    //set_variable_speed(g, cfg.n, cfg.m, cfg.h);

    // Adding the walls
    add_walls(g, cfg);

    // Removing sources inside walls
    removing_sources_in_walls(g, &cfg);

    if (cfg.nsources == 0){
        printf("Error: All sources are inside walls, nothing to propagate.\n");
        eikonal_grid_free(g);
        free_config(&cfg);
        return 1;
    }

    // Storing source information
    save_sources(cfg);

    save_speed(g, "speed.txt");

    // Execution of FIM
    fim_solve(g, cfg.src_i, cfg.src_j, cfg.nsources, EPSILON, 30.0);

    eikonal_save_matrix(g, "matrix_fim.txt");

    // Cleaning
    eikonal_grid_free(g);
    free_config(&cfg);
    
    return 0;
}