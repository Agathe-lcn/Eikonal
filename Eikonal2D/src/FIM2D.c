#include "../include/FIM2D.h"

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

// Liste doublement chaînée

typedef struct Node{
    int index;  // Indice i * m + j
    struct Node* prev;
    struct Node* next;
}Node;

typedef struct{
    Node head;
    Node** index_to_node;
    int ncell;  // ncell = n*m
    int size;   // Nombre d'éléments dans la liste
}NodeList;

static NodeList* list_create(int ncell){
    NodeList* l = (NodeList *)malloc(sizeof(NodeList));
    l->index_to_node = (Node **)calloc(ncell, sizeof(Node *));
    l->ncell = ncell;
    l->size = 0;
    l->head.prev = &l->head;
    l->head.next = &l->head;
    return l;
}

static void list_free(NodeList *l){
    for (int k=0; k < l->ncell; k++){
        if (l->index_to_node[k])
            free(l->index_to_node[k]);
    }
    free(l->index_to_node);
    free(l);
}

// Insert index into the list before the 'node' node
static void list_insert_before(NodeList* l, int index, Node* node){
    // ALready present
    if (l->index_to_node[index])
        return;

    Node* new = (Node *)malloc(sizeof(Node));
    new->index = index;
    new->next = node;
    new->prev = node->prev;
    node->prev->next = new;
    node->prev = new;
    l->index_to_node[index] = new;
    l->size++;
}

// Insére un indice au début de la liste
static void list_push_front(NodeList* l, int index){
    list_insert_before(l, index, l->head.next);
}

// Insére un indice à la fin de la liste
static void list_push_back(NodeList* l, int index){
    list_insert_before(l, index, &l->head);
}

// Supprime le noeud qui contient l'indice 'index' de la liste
static void list_remove(NodeList* l, int index){
    Node* node = l->index_to_node[index];
    if (!node)
        return;

    node->prev->next = node->next;
    node->next->prev = node->prev;
    free(node);
    l->index_to_node[index] = NULL;
    l->size--;
}

// Récupère et renvoie le premier élément de la liste
static int list_pop_front(NodeList* l){
    if (l->size == 0)
        return -1;
    
    Node* first = l->head.next;
    int index = first->index;
    list_remove(l, index);
    return index;
}

static inline int list_contains(const NodeList* l, int index){
    return l->index_to_node[index] != NULL;
}

static inline int list_is_empty(const NodeList* l){
    return l->size == 0;
}

// Vérifie si un voisin existe et peut être amélioré
static inline int is_neighbor_updatable(const EikonalGrid* g, int index, double T_new, double epsilon){
    if (index < 0)
        return 0;

    double T_old = g->T[index];

    if (T_old - T_new > epsilon)
        return 1;

    return 0;
}
/*
// Ajoute un voisin à la liste s'il n'y est pas déjà
static void add_neighbor_if_needed(NodeList* list, const EikonalGrid* g, int index, double T_new, Node* current_node, double epsilon){
    if (index < 0)
        return;

    if (!list_contains(list, index)){
        double T_old = g->T[index];
        if (T_old - T_new > epsilon)
            list_insert_before(list, index, current_node);
    }
}*/


// Vérifie si une cellule est dans le rayon
static int is_in_radius(const int* src_i, const int* src_j, int ns, int i, int j, int source_tag, double max_radius){
    if (max_radius <= 0.0)
        return 1;

    double di = (double)(i - src_i[source_tag]);
    double dj = (double)(j - src_j[source_tag]);
    double dist = sqrt(di * di + dj * dj);
    return (dist <= max_radius);
}


// Compte le nombre de cellules dans le rayon
static int count_cells_in_radius(const EikonalGrid* g, const int* src_i, const int* src_j, int ns, double max_radius) {
    if (max_radius <= 0.0)
        return g->n * g->m;
    
    int n = g->n;
    int m = g->m;
    int count = 0;
    int radius_cells = (int)ceil(max_radius);
    
    int* visited = (int*)calloc(n * m, sizeof(int));
    if (!visited)
        return g->n * g->m;
    
    for (int s = 0; s < ns; s++) {
        int i0 = src_i[s];
        int j0 = src_j[s];
        
        for (int di = -radius_cells; di <= radius_cells; di++) {
            for (int dj = -radius_cells; dj <= radius_cells; dj++) {
                int i = i0 + di;
                int j = j0 + dj;
                
                if (i >= 0 && i < n && j >= 0 && j < m) {
                    double dist = sqrt((double)(di * di + dj * dj));
                    if (dist <= max_radius) {
                        int index = i * m + j;
                        if (!visited[index]) {
                            visited[index] = 1;
                            count++;
                        }
                    }
                }
            }
        }
    }
    
    free(visited);
    return count;
}


// Détermine le tag de la source qui a mis à jour la valeur T en (i,j)
static int find_tag(const EikonalGrid* g, const int* source_tag, int i, int j, int n, int m){
    int neighbors[4][2] = {{i-1, j}, {i+1, j}, {i, j-1}, {i, j+1}};
    int tag = -1;
    double T_min = EIKONAL_INF;

    for (int k = 0; k < 4; k++){
        int ni = neighbors[k][0];
        int nj = neighbors[k][1];

        if (ni < 0 || ni >= n || nj < 0 || nj >= m)
            continue;

        int index_neighbor = ni * m + nj;
        if (source_tag[index_neighbor] < 0)
            continue;

        if (g->T[index_neighbor] < T_min){
            T_min = g->T[index_neighbor];
            tag = source_tag[index_neighbor];
        }
    }

    return tag;
}


void fim_solve(EikonalGrid* g, const int* src_i, const int* src_j, int ns, double epsilon, double max_radius) {
    if (!g || !src_i || !src_j || ns <= 0)
        return;

    clock_t start_time = clock();

    int n = g->n;
    int m = g->m;
    int ncell = n * m;

    // Alloue et initialise le tag
    int* source_tag = (int*)malloc(ncell * sizeof(int));
    if (!source_tag){
        printf("Erreur: Impossible d'allouer de la mémoire pour les tags\n");
        return;
    }
    for (int k = 0; k < ncell; k++)
        source_tag[k] = -1;

    // Initialisation: définit toutes les celulles à +inf
    for (int k=0; k < ncell; k++)
        g->T[k] = EIKONAL_INF;

    // Initialisation des sources
    for (int s=0; s < ns; s++){
        int i = src_i[s];
        int j = src_j[s];
        if (i >= 0 && i < n && j >= 0 && j < m){
            int index = i * m + j;
            g->T[index] = 0.0;
            source_tag[index] = s;
        }
    }

    // Création de la Narrowband
    NodeList* narrow = list_create(ncell);
    if (!narrow) {
        free(source_tag);
        return;
    }

    // Ajoute les voisins des sources à la Narrowband
    for (int s=0; s < ns; s++){
        int i = src_i[s];
        int j = src_j[s];
        if (i < 0 || i >= n || j < 0 || j >= m)
            continue;

        // Visite des 4 voisins
        int neighbors[4][2] = {{i-1, j}, {i+1, j}, {i, j-1}, {i, j+1}};
        for (int k=0; k<4; k++){
            int ni = neighbors[k][0];
            int nj = neighbors[k][1];

            if (ni >= 0 && ni < n && nj >= 0 && nj < m){
                int index = ni * m + nj;

                // Vérifie si le voisin se trouve à l'intérieur du cercle de rayon max_seuil et de centre s
                if (!is_in_radius(src_i, src_j, ns, ni, nj, s, max_radius))
                    continue;

                if (!list_contains(narrow, index)){
                    double T_new = eikonal_solve_local(g, ni, nj);
                    if (T_new < g->T[index] - 1e-12){
                        g->T[index] = T_new;
                        source_tag[index] = s;
                        list_push_back(narrow, index);
                    }
                }
            }
        }
    }


    int compute_all = (max_radius <= 0.0);
    int total_cells = 0;
    int cells_processed = 0;

    // Compte le nombre de cellules dans le rayon
    if (!compute_all) 
        total_cells = count_cells_in_radius(g, src_i, src_j, ns, max_radius);

    // Boucle principale
    while(!list_is_empty(narrow)) {
        // Enlève le premier élément de la liste
        int index = list_pop_front(narrow);
        if (index < 0)
            continue;

        int i = index / m;
        int j = index % m;

        // Vérifie le rayon
        if (!compute_all && !is_in_radius(src_i, src_j, ns, i, j, source_tag[index], max_radius))
            continue; 

        // Incrémente le compteur
        if (!compute_all && (g->T[index] > 0 || source_tag[index] < 0)) {
            cells_processed++;
        }

        double T_old = g->T[index];
        double T_new = eikonal_solve_local(g, i, j);
        double diff = fabs(T_new - T_old);

        if (diff <= epsilon) {
            // La cellule a convergé: on la fige et ses voisins susceptibles d'être améliorés sont ajoutés à la liste
            int neighbors[4][2] = {{i-1, j}, {i+1, j}, {i, j-1}, {i, j+1}};
            for (int k=0; k < 4; k++){
                int ni = neighbors[k][0];
                int nj = neighbors[k][1];
                if (ni >= 0 && ni < n && nj >= 0 && nj < m){
                    int index_neighbor = ni * m + nj;

                    // Vérifie le rayon
                    if (!compute_all && !is_in_radius(src_i, src_j, ns, ni, nj, source_tag[index], max_radius))
                        continue;

                    // Vérifie si l evoisin peut être amélioré
                    double T_neighbor_new = eikonal_solve_local(g, ni, nj);
                    if (T_neighbor_new < g->T[index_neighbor] - 1e-12){
                        g->T[index_neighbor] = T_neighbor_new;
                        source_tag[index_neighbor] = find_tag(g, source_tag, ni, nj, n, m);
                        if (!list_contains(narrow, index_neighbor))
                            list_push_back(narrow, index_neighbor);
                    }
                }
            }
        } else {
            // La cellule n'a pas convergé: mise à jour de sa valeur
            g->T[index] = T_new;
            source_tag[index] = find_tag(g, source_tag, i, j, n, m);

            // La cellule est réinsérée dans la liste (elle sera recalculée)
            list_push_front(narrow, index);

            // Visite des 4 voisins pour les ajouter s'ils peuvent être améliorés
            int neighbors[4][2] = {{i-1, j}, {i+1, j}, {i, j-1}, {i, j+1}};
            for (int k=0; k < 4; k++){
                int ni = neighbors[k][0];
                int nj = neighbors[k][1];
                if (ni >= 0 && ni < n && nj >= 0 && nj < m){
                    int index_neighbor = ni * m + nj;

                    // Vérifie le rayon
                    if (!compute_all && !is_in_radius(src_i, src_j, ns, ni, nj, source_tag[index], max_radius))
                        continue;

                    // Vérifie si le voisin peut être amélioré
                    double T_neighbor_new = eikonal_solve_local(g, ni, nj);
                    if (T_neighbor_new < g->T[index_neighbor] - 1e-12){
                        g->T[index_neighbor] = T_neighbor_new;
                        source_tag[index_neighbor] = find_tag(g, source_tag, ni, nj, n, m);
                        if (!list_contains(narrow, index_neighbor))
                            list_push_back(narrow, index_neighbor);
                    }
                }
            }
        }

        // Vérifie les conditions d'arrêt
        if (!compute_all && cells_processed >= total_cells)
            break;
    }

    // Enregistrement des tags
    eikonal_save_tags(g, source_tag, "source_tags.txt");

    clock_t end_time = clock();
    double time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

    if (compute_all)
        printf("FIM terminée (sans seuil) en %.6f secondes\n", time);
    else
        printf("FIM terminée (avec seuil): %d cellules traitées en %.6f secondes\n", cells_processed, time);

    // Nettoyage
    list_free(narrow);
    free(source_tag);
}