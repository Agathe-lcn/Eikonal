#include "../include/FIM2D_mpi.h"
#include "../include/FIM2D_io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <limits.h>

#define EPSILON 1e-12
#define MAX_LINE 1024
#define TEST_MODE 0     // 1 pour désactiver les communications et 0 pour les activer

static int parse_named_int_exact(const char* line, const char* key, int* out_value){
    size_t key_len = strlen(key);
    const char* cursor = line;
    char* end_ptr = NULL;
    long value;

    while (*cursor == ' ' || *cursor == '\t')
        cursor++;

    if (strncmp(cursor, key, key_len) != 0)
        return 0;
    cursor += key_len;

    while (*cursor == ' ' || *cursor == '\t')
        cursor++;
    if (*cursor != '=')
        return -1;
    cursor++;

    while (*cursor == ' ' || *cursor == '\t')
        cursor++;
    if (*cursor == '\0')
        return -1;

    value = strtol(cursor, &end_ptr, 10);
    if (end_ptr == cursor)
        return -1;
    while (*end_ptr == ' ' || *end_ptr == '\t')
        end_ptr++;
    if (*end_ptr != '\0')
        return -1;
    if (value < INT_MIN || value > INT_MAX)
        return -1;

    *out_value = (int)value;
    return 1;
}

// Lecture du fichier de configuration
Config2 read_config_mpi(const char* filename){
    Config2 cfg = {};
    FILE* file = fopen(filename, "r");
    if (!file){
        printf("Erreur: Impossible d'ouvrir %s\n", filename);
        return cfg;
    }

    char line[MAX_LINE];
    int section = 0;    // 0: aucune, 1: sources, 2: murs
    int max_sources = 2000;
    int max_walls = 2000;

    // Allocations
    cfg.src_i = (int*)malloc(max_sources * sizeof(int));
    cfg.src_j = (int*)malloc(max_sources * sizeof(int));
    cfg.src_i_overlap = (int*)malloc(max_sources * sizeof(int));
    cfg.src_j_overlap = (int*)malloc(max_sources * sizeof(int));
    cfg.wall_c1 = (int*)malloc(max_walls * sizeof(int));
    cfg.wall_c2 = (int*)malloc(max_walls * sizeof(int));
    cfg.wall_r1 = (int*)malloc(max_walls * sizeof(int));
    cfg.wall_r2 = (int*)malloc(max_walls * sizeof(int));

    cfg.nsources_overlap = 0;
    cfg.valid = true;

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

        // Lecture des paramètres de la grille
        int parsed = parse_named_int_exact(line, "overlap", &cfg.overlap);
        if (parsed == 1)
            continue;
        if (parsed < 0){
            printf("Erreur: overlap doit etre un entier strictement formate dans %s.\n", filename);
            cfg.valid = false;
            break;
        }
        if (sscanf(line, "n = %d", &cfg.n) == 1)
            continue;
        if (sscanf(line, "m = %d", &cfg.m) == 1)
            continue;
        if (sscanf(line, "h = %lf", &cfg.h) == 1)
            continue;

        // Lecture des sources
        if (section == 1){
            int i, j;
            if (sscanf(line, "%d %d",&i, &j) == 2){
                if (cfg.nsources < max_sources){
                    cfg.src_i[cfg.nsources] = i;
                    cfg.src_j[cfg.nsources] = j;
                    cfg.nsources++;
                }
            }
        }

        // Lecture des murs
        if (section == 2){
            int c1, c2, r1, r2;
            if (sscanf(line, "%d %d %d %d", &c1, &c2, &r1, &r2) == 4){
                if (cfg.nwalls < max_walls){
                    cfg.wall_c1[cfg.nwalls] = c1;
                    cfg.wall_c2[cfg.nwalls] = c2;
                    cfg.wall_r1[cfg.nwalls] = r1;
                    cfg.wall_r2[cfg.nwalls] = r2;
                    cfg.nwalls++;
                }
            }
        }
    }

    fclose(file);
    return cfg;
}


// Nettoyage de la configuration
void free_config_mpi(Config2* cfg){
    free(cfg->src_i);
    free(cfg->src_j);
    free(cfg->src_i_overlap);
    free(cfg->src_j_overlap);
    free(cfg->wall_c1);
    free(cfg->wall_c2);
    free(cfg->wall_r1);
    free(cfg->wall_r2);
}


// Ajout des murs présents dans le sous-domaine du processus 
void add_walls_local(EikonalGrid* g, Config2 cfg, MPIDomain* domain){
    for (int w=0; w < cfg.nwalls; w++){
        int c1 = cfg.wall_c1[w];
        int c2 = cfg.wall_c2[w];
        int r1 = cfg.wall_r1[w];
        int r2 = cfg.wall_r2[w];

        // Vérification des limites globales
        if (c1 < 0) 
            c1 = 0;
        if (c2 >= cfg.m) 
            c2 = cfg.m - 1;
        if (r1 < 0) 
            r1 = 0;
        if (r2 >= cfg.n) 
            r2 = cfg.n - 1;

        // Vérifier si le mur est en dehors du sous-domaine
        if (r2 < domain->i_start_overlap || r1 > domain->i_end_overlap)
            continue;

        // Calculer les indices locaux
        int local_r1 = r1 - domain->i_start_overlap;
        int local_r2 = r2 - domain->i_start_overlap;

        // Ajuster les indices locaux pour rester dans la zone
        if (local_r1 < 0) 
            local_r1 = 0;
        if (local_r2 >= domain->n_overlap) 
            local_r2 = domain->n_overlap - 1;

        // Ajouter le mur dans la zone locale
        for (int i=local_r1; i <= local_r2; i++){
            for (int j=c1; j <= c2; j++){
                eikonal_grid_set_obstacle(g, i, j);
            }
        }
    }
}

// Suppression des sources qui sont dans un mur présent dans le sous-domaine du processus
int removing_sources_in_walls_local(EikonalGrid* g_processus, Config2* cfg, MPIDomain* domain){
    int valid = 0;
    int local_count = 0;

    // TEST
    //printf("========== DEBUT removing_sources_in_walls_local ==========\n");
    //printf("[Processus %d] cfg->nsources = %d\n", domain->rank, cfg->nsources);
    //printf("[Processus %d] domain->i_start_overlap=%d, domain->i_end_overlap=%d\n", domain->rank, domain->i_start_overlap, domain->i_end_overlap);


    for (int s=0; s < cfg->nsources; s++){
        int i_global = cfg->src_i[s];
        int j_global = cfg->src_j[s];

        // TEST
        //printf("[Processus %d] Source %d globale: (%d, %d)\n", domain->rank, s, i_global, j_global);

        // 1er cas: la source est en dehors du sous-domaine du processus
        if (i_global < domain->i_start_overlap || i_global > domain->i_end_overlap){

            // TEST
            //printf("[Processus %d]   -> HORS ZONE, on la garde\n", domain->rank);

            // On la conserve telle qu'elle est car elle sera traitée par un autre processus
            cfg->src_i[valid] = i_global;
            cfg->src_j[valid] = j_global;
            valid++;
            continue;
        }


        // 2e cas: la source est dans le sous-domaine du processus
        int i_local = i_global - domain->i_start_overlap;
        int j_local = j_global;

        // TEST
        //printf("[Processus %d]   -> DANS ZONE -> locale (%d, %d)\n", domain->rank, i_local, j_local);

        // On vérifie si la source est dans un mur
        if (eikonal_grid_is_obstacle(g_processus, i_local, j_local)){
            printf("Attention: La source %d située aux coordonnées (%d,%d) est dans un mur, elle être ignorée.\n", s, i_global, j_global);
            continue;
        }

        // Si la source est valide, on la conserve puis on l'ajoute aux sources locales

        // TEST
        //printf("[Processus %d]   -> VALIDE\n", domain->rank);

        cfg->src_i[valid] = i_global;
        cfg->src_j[valid] = j_global;
        valid++;
        cfg->src_i_overlap[local_count] = i_local;
        cfg->src_j_overlap[local_count] = j_local;
        local_count++;
    }
    cfg->nsources = valid;
    cfg->nsources_overlap = local_count;

    // TEST
    //printf("[Processus %d] RESULTAT: %d sources valides, %d sources locales\n", domain->rank, valid, local_count);
    //printf("========== FIN removing_sources_in_walls_local ==========\n");

    return local_count;
}


// Enregistrement des coordonnées des sources
void save_sources_mpi(Config2 cfg, int rank){
    if (rank == 0){
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
}


void save_local_result_mpi(const MPIDomain* domain, const EikonalGrid* g_processus){
    char matrix_filename[256];
    char meta_filename[256];
    snprintf(matrix_filename, sizeof(matrix_filename), "result_rank%d.txt", domain->rank);
    snprintf(meta_filename, sizeof(meta_filename), "result_rank%d_meta.txt", domain->rank);

    FILE* matrix_file = fopen(matrix_filename, "w");
    if (!matrix_file){
        printf("[Processus %d] Erreur: impossible de créer %s\n", domain->rank, matrix_filename);
        return;
    }

    for (int i = 0; i < g_processus->n; i++){
        for (int j = 0; j < g_processus->m; j++){
            fprintf(matrix_file, "%.17g", g_processus->T[i * g_processus->m + j]);
            if (j + 1 < g_processus->m)
                fputc(' ', matrix_file);
        }
        fputc('\n', matrix_file);
    }
    fclose(matrix_file);

    FILE* meta_file = fopen(meta_filename, "w");
    if (!meta_file){
        printf("[Processus %d] Erreur: impossible de créer %s\n", domain->rank, meta_filename);
        return;
    }

    fprintf(meta_file, "rank=%d\n", domain->rank);
    fprintf(meta_file, "n_global=%d\n", domain->n);
    fprintf(meta_file, "m_global=%d\n", domain->m);
    fprintf(meta_file, "h=%0.17g\n", domain->h);
    fprintf(meta_file, "overlap=%d\n", domain->overlap);
    fprintf(meta_file, "n_local=%d\n", domain->n_overlap);
    fprintf(meta_file, "top_ghost=%d\n", domain->top_ghost);
    fprintf(meta_file, "bottom_ghost=%d\n", domain->bottom_ghost);
    fprintf(meta_file, "i_start_overlap=%d\n", domain->i_start_overlap);
    fprintf(meta_file, "i_end_overlap=%d\n", domain->i_end_overlap);
    fprintf(meta_file, "i_owned_start=%d\n", domain->i_owned_start);
    fprintf(meta_file, "i_owned_end=%d\n", domain->i_owned_end);
    fprintf(meta_file, "n_owned=%d\n", domain->n_owned);
    fclose(meta_file);
}


void initialize_grid_with_sources(EikonalGrid* g_processus, Config2* cfg_processus, MPIDomain* domain, int* start){
    int n = g_processus->n;
    int m = g_processus->m;
    int ncell = n * m;

    // TEST
    //printf("========== DEBUT initialize_grid_with_sources ==========\n");
    //printf("[Processus %d] n=%d, m=%d, ncell=%d\n", domain->rank, n, m, ncell);
    
    // Initialisation: définit toutes les mailles à +inf
    for (int k = 0; k < ncell; k++)
        g_processus->T[k] = EIKONAL_INF;

    // Initialisation des sources locales
    int ns_local = cfg_processus->nsources_overlap;
    int* src_i_overlap = cfg_processus->src_i_overlap;
    int* src_j_overlap = cfg_processus->src_j_overlap;

    // TEST
    //printf("[Processus %d] ns_local = %d\n", domain->rank, ns_local);
    
    for (int s = 0; s < ns_local; s++){
        int i = src_i_overlap[s];
        int j = src_j_overlap[s];

        // TEST
        //printf("[Processus %d] Source %d: (%d, %d)\n", domain->rank, s, i, j);
        
        // Vérification des limites
        if (i < 0 || i >= n || j < 0 || j >= m){

            // TEST
            //printf("[Processus %d] ERREUR: Source %d hors limites! (%d, %d) n=%d m=%d\n", domain->rank, s, i, j, n, m);

            continue;
        }
        
        int index = i * m + j;
        g_processus->T[index] = 0.0;
        start[index] = 1;   // On marque la source comme point de départ

        // TEST
        //printf("[Processus %d] Source %d -> index=%d, T=0, start=1\n", domain->rank, s, index);
    }

    // TEST
    /*int nb_start = 0;
    int nb_zeros = 0;
    for (int k = 0; k < ncell; k++) {
        if (start[k]) nb_start++;
        if (g_processus->T[k] == 0.0) nb_zeros++;
    }
    printf("[Processus %d] FIN initialize: start=%d, T zeros=%d\n", domain->rank, nb_start, nb_zeros);
    printf("========== FIN initialize_grid_with_sources ==========\n");*/
}

// 2 paramètres : le premier correspond à l'overlap et le second au fichier de configuration
int main(int argc, char** argv){
    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    FIMIO_Init(MPI_COMM_WORLD);

    Config2 cfg_processus;
    int overlap = 0;

    if (rank == 0){

        // Nom du fichier de configuration
        const char* config_file = "config.txt";
        if (argc == 2){
            config_file = argv[1];
        }
        else if (argc > 2){
            printf("Erreur: arguments invalides. Usage: generate_grid_mpi [config.txt]\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Lecture de la configuration
        cfg_processus = read_config_mpi(config_file);

        if (!cfg_processus.valid){
            free_config_mpi(&cfg_processus);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (cfg_processus.n == 0 || cfg_processus.m == 0){
            printf("Erreur: configuration invalide (n ou m pas défini)\n");
            free_config_mpi(&cfg_processus);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (cfg_processus.h <= 0){
            printf("Erreur: configuration invalide (h pas défini)\n");
            free_config_mpi(&cfg_processus);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (cfg_processus.nsources == 0){
            printf("Erreur: aucune source spécifiée\n");
            free_config_mpi(&cfg_processus);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (cfg_processus.overlap < 1){
            printf("Erreur: configuration invalide (overlap doit etre un entier >= 1)\n");
            free_config_mpi(&cfg_processus);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        overlap = cfg_processus.overlap;
    }

    // Distribution de overlap et de certains champs de cfg_processus à tous les processus
    MPI_Bcast(&overlap, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cfg_processus.n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cfg_processus.m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cfg_processus.h, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cfg_processus.nsources, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cfg_processus.nwalls, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Allocation et diffusion des tableaux de sources
    if (rank != 0){
        cfg_processus.src_i = (int*)malloc(cfg_processus.nsources * sizeof(int));
        cfg_processus.src_j = (int*)malloc(cfg_processus.nsources * sizeof(int));
        cfg_processus.wall_c1 = (int*)malloc(cfg_processus.nwalls * sizeof(int));
        cfg_processus.wall_c2 = (int*)malloc(cfg_processus.nwalls * sizeof(int));
        cfg_processus.wall_r1 = (int*)malloc(cfg_processus.nwalls * sizeof(int));
        cfg_processus.wall_r2 = (int*)malloc(cfg_processus.nwalls * sizeof(int));
        // Allocation des tableaux locaux
        cfg_processus.src_i_overlap = (int*)malloc(cfg_processus.nsources * sizeof(int));
        cfg_processus.src_j_overlap = (int*)malloc(cfg_processus.nsources * sizeof(int));
    }

    // Distribution des autres champs de cfg_processus (sauf les champs locaux)
    MPI_Bcast(cfg_processus.src_i, cfg_processus.nsources, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cfg_processus.src_j, cfg_processus.nsources, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cfg_processus.wall_c1, cfg_processus.nwalls, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cfg_processus.wall_c2, cfg_processus.nwalls, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cfg_processus.wall_r1, cfg_processus.nwalls, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cfg_processus.wall_r2, cfg_processus.nwalls, MPI_INT, 0, MPI_COMM_WORLD);

    // Définition de la topologie
    MPIDomain* domain = topology_create(cfg_processus.n, cfg_processus.m, cfg_processus.h, overlap);
    if (!domain){
        free_config_mpi(&cfg_processus);
        FIMIO_Finalize();
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Création de la grille
    EikonalGrid* g_processus = eikonal_grid_create(domain->n_overlap, domain->m, domain->h);
    if (!g_processus){
        printf("Erreur [rang %d]: impossible de créer la grille\n", rank);
        free_config_mpi(&cfg_processus);
        topology_free(domain);
        MPI_Finalize();
        return 1;
    }

    // Vitesse constante égale à 1
    eikonal_grid_set_speed_constant(g_processus, 1.0);

    // On ajoute les murs
    add_walls_local(g_processus, cfg_processus, domain);

    // On supprime les sources dans le mur
    int ns_local = removing_sources_in_walls_local(g_processus, &cfg_processus, domain);
    cfg_processus.nsources_overlap = ns_local;

    if (cfg_processus.nsources == 0){
        printf("Erreur: toutes les sources sont dans un mur, il n'y a rien à propager.\n");
        eikonal_grid_free(g_processus);
        free_config_mpi(&cfg_processus);
        topology_free(domain);
        MPI_Finalize();
        return 1;
    }

    // Initialisation des sources
    int ncell = domain-> n_overlap * domain->m;
    int* start = (int*)calloc(ncell, sizeof(int));
    initialize_grid_with_sources(g_processus, &cfg_processus, domain, start);

    // Stockage des informations sur les sources
    save_sources_mpi(cfg_processus, rank);

    // Exécution de la FIM
    fim_solve_mpi(domain, g_processus, &cfg_processus, start, EPSILON, -1);

    // Sauvegarde du résultat local de chaque rang pour la visualisation par sous-domaine.
    save_local_result_mpi(domain, g_processus);

    // Dans le cas où on veut tester dans les communications, chaque processus sauvegarde ses résultats dans un fichier différent
    #if TEST_MODE
        char filename[256];
        snprintf(filename, sizeof(filename), "result_test_rank%d.txt", rank);
        
        FILE* file = fopen(filename, "w");
        if (file){
            int n = g_processus->n;
            int m = g_processus->m;
            for (int i = 0; i < n; i++){
                for (int j = 0; j < m; j++){
                    fprintf(file, "%.6f ", g_processus->T[i * m + j]);
                }
                fprintf(file, "\n");
            }
            fclose(file);
        }
        else {
            printf("[Processus %d] Erreur: impossible de sauvegarder %s\n", rank, filename);
        }
    #endif

    // Sauvegarde des résultats
    MPI_Info info = MPI_INFO_NULL;
    int err = FIMIO_Checkpoint("matrix_fim.txt", domain, g_processus, -1, info);
    if (err != MPI_SUCCESS && rank == 0)
        printf("Erreur lors de l'enregistrement des résultats de la FIM");

    // Nettoyage
    free(start);
    eikonal_grid_free(g_processus);
    free_config_mpi(&cfg_processus);

    FIMIO_Finalize();

    MPI_Finalize();
    
    return 0;
}