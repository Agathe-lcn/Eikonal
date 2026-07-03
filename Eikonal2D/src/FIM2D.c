#include "../include/FIM2D.h"

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

// Double-linked list for narrowband

typedef struct Node{
    int index;  // Index i * m + j
    struct Node* prev;
    struct Node* next;
}Node;

typedef struct{
    Node head;
    Node** index_to_node;
    int ncell;  // ncell = n*m
    int size;   // Number of items currently in the list
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

// Insert index at the beginning of the list
static void list_push_front(NodeList* l, int index){
    list_insert_before(l, index, l->head.next);
}

// Insert index at the end of the list
static void list_push_back(NodeList* l, int index){
    list_insert_before(l, index, &l->head);
}

// Removes the node with index from the list
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

// Retrieves and returns the first item in the list
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

// Check whether a neighbor exists and can be improved
static inline int is_neighbor_updatable(const EikonalGrid* g, int index, double T_new, double epsilon){
    if (index < 0)
        return 0;

    double T_old = g->T[index];

    if (T_old - T_new > epsilon)
        return 1;

    return 0;
}

// Add a neighbor to the list if it isn't already on it
static void add_neighbor_if_needed(NodeList* list, const EikonalGrid* g, int index, double T_new, Node* current_node, double epsilon){
    if (index < 0)
        return;

    if (!list_contains(list, index)){
        double T_old = g->T[index];
        if (T_old - T_new > epsilon)
            list_insert_before(list, index, current_node);
    }
}


void fim_solve(EikonalGrid* g, const int* src_i, const int* src_j, int ns, double epsilon){
    if (!g || !src_i || !src_j || ns <= 0)
        return;

    int n = g->n;
    int m = g->m;
    int ncell = n * m;

    // Tag allocation and initialization
    int* source_tag = (int*)malloc(ncell * sizeof(int));
    if (!source_tag){
        printf("Error: Unable to allocate memory for source tags\n");
        return;
    }
    for (int k = 0; k < ncell; k++)
        source_tag[k] = -1;

    // Initialization: set all cells to +inf
    for (int k=0; k < ncell; k++)
        g->T[k] = EIKONAL_INF;

    // Initializing sources
    for (int s=0; s < ns; s++){
        int i = src_i[s];
        int j = src_j[s];
        if (i >= 0 && i < n && j >= 0 && j < m){
            int index = i * m + j;
            g->T[index] = 0.0;
            source_tag[index] = s;
        }
    }

    // Creation of Narrow Band
    NodeList* narrow = list_create(ncell);
    if (!narrow)
        return;

    // Add the neighbors of the sources to the active list
    for (int s=0; s < ns; s++){
        int i = src_i[s];
        int j = src_j[s];
        if (i < 0 || i >= n || j < 0 || j >= m)
            continue;

        // The 4 neighbors tour
        int neighbors[4][2] = {{i-1, j}, {i+1, j}, {i, j-1}, {i, j+1}};
        for (int k=0; k<4; k++){
            int ni = neighbors[k][0];
            int nj = neighbors[k][1];

            if (ni >= 0 && ni < n && nj >= 0 && nj < m){
                int index = ni * m + nj;
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

    // Main loop
    int iterations = 0;
    while(!list_is_empty(narrow) && iterations < ncell * 10){
        iterations ++;

        // Remove the first item of the list
        int index = list_pop_front(narrow);
        if (index < 0)
            continue;

        int i = index/m;
        int j = index%m;
        double T_old = g->T[index];
        double T_new = eikonal_solve_local(g, i, j);
        double diff = fabs(T_new - T_old);

        if (diff <= epsilon){
            // The cell has converged: we freeze it and its neighbors that can be improved are added to the list
            int neighbors[4][2] = {{i-1, j}, {i+1, j}, {i, j-1}, {i, j+1}};
            for (int k=0; k < 4; k++){
                int ni = neighbors[k][0];
                int nj = neighbors[k][1];
                if (ni >= 0 && ni < n && nj >= 0 && nj < m){
                    int index_neighbor = ni * m + nj;

                    // Check if the neighbor can be improved
                    double T_neighbor_new = eikonal_solve_local(g, ni, nj);
                    if (T_neighbor_new < g->T[index_neighbor] - 1e-12){
                        g->T[index_neighbor] = T_neighbor_new;
                        source_tag[index_neighbor] = source_tag[index];
                        if (!list_contains(narrow, index_neighbor))
                            list_push_back(narrow, index_neighbor);
                    }
                }
            }
        }

        else{
            // The cell has not converged: update its value
            g->T[index] = T_new;

            // The cell is reinserted into the list (it will be recalculated)
            list_push_front(narrow, index);

            // We go through the four neighbors to add them if they can be upgraded
            int neighbors[4][2] = {{i-1, j}, {i+1, j}, {i, j-1}, {i, j+1}};
            for (int k=0; k < 4; k++){
                int ni = neighbors[k][0];
                int nj = neighbors[k][1];
                if (ni >= 0 && ni < n && nj >= 0 && nj < m){
                    int index_neighbor = ni * m + nj;

                    // Check if the neighbor can be improved
                    double T_neighbor_new = eikonal_solve_local(g, ni, nj);
                    if (T_neighbor_new < g->T[index_neighbor] - 1e-12){
                        g->T[index_neighbor] = T_neighbor_new;
                        source_tag[index_neighbor] = source_tag[index];
                        if (!list_contains(narrow, index_neighbor))
                            list_push_back(narrow, index_neighbor);
                    }
                }
            }
        }
    }

    // Saving tags
    eikonal_save_tags(g, source_tag, "source_tags.txt");

    // Cleaning
    list_free(narrow);
}