#include "../include/FIM2D.h"

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

// Fonctions pour la liste doublement chaînée

NodeList* list_create(int ncell){
    NodeList* l = (NodeList *)malloc(sizeof(NodeList));
    l->index_to_node = (Node **)calloc(ncell, sizeof(Node *));
    l->ncell = ncell;
    l->size = 0;
    l->head.prev = &l->head;
    l->head.next = &l->head;
    return l;
}

void list_free(NodeList *l){
    for (int k=0; k < l->ncell; k++){
        if (l->index_to_node[k])
            free(l->index_to_node[k]);
    }
    free(l->index_to_node);
    free(l);
}

void list_insert_before(NodeList* l, int index, Node* node){
    // Déjà présent
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

void list_push_front(NodeList* l, int index){
    list_insert_before(l, index, l->head.next);
}

void list_push_back(NodeList* l, int index){
    list_insert_before(l, index, &l->head);
}

void list_remove(NodeList* l, int index){
    Node* node = l->index_to_node[index];
    if (!node)
        return;

    node->prev->next = node->next;
    node->next->prev = node->prev;
    free(node);
    l->index_to_node[index] = NULL;
    l->size--;
}

int list_pop_front(NodeList* l){
    if (l->size == 0)
        return -1;
    
    Node* first = l->head.next;
    int index = first->index;
    list_remove(l, index);
    return index;
}

int list_contains(const NodeList* l, int index){
    return l->index_to_node[index] != NULL;
}

int list_is_empty(const NodeList* l){
    return l->size == 0;
}

int is_neighbor_updatable(const EikonalGrid* g, int index, double T_new, double epsilon){
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

int find_tag(const EikonalGrid* g, const int* source_tag, int i, int j, int n, int m){
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


void fim_solve(EikonalGrid* g, const int* src_i, const int* src_j, int ns, double epsilon, int max_depth) {
    if (!g || !src_i || !src_j || ns <= 0)
        return;

    int n = g->n;
    int m = g->m;
    int ncell = n * m;
    int use_depth_limit = (max_depth >= 0);

    // Alloue et initialise le tag
    int* source_tag = (int*)malloc(ncell * sizeof(int));
    int* cell_depth = (int*)malloc(ncell * sizeof(int));
    if (!source_tag || !cell_depth){
        printf("Erreur: Impossible d'allouer de la mémoire pour les tags\n");
        free(source_tag);
        free(cell_depth);
        return;
    }
    for (int k = 0; k < ncell; k++){
        source_tag[k] = -1;
        cell_depth[k] = -1;
    }

    // Initialisation: définit toutes les mailles à +inf
    for (int k=0; k < ncell; k++)
        g->T[k] = EIKONAL_INF;

    // Création de la Narrowband
    NodeList* narrow = list_create(ncell);
    if (!narrow) {
        free(source_tag);
        free(cell_depth);
        return;
    }

    // Initialisation des sources et ajout de leurs voisins à la narrowband
    for (int s=0; s < ns; s++){
        int i = src_i[s];
        int j = src_j[s];
        if (i < 0 || i >= n || j < 0 || j >= m)
            continue;
        int index = i * m + j;
        g->T[index] = 0.0;
        source_tag[index] = s;
        cell_depth[index] = 0;

        int neighbors[4][2] = {{i-1, j}, {i+1, j}, {i, j-1}, {i, j+1}};
        for (int k=0; k<4; k++){
            int ni = neighbors[k][0];
            int nj = neighbors[k][1];

            if (ni >= 0 && ni < n && nj >= 0 && nj < m){
                int index = ni * m + nj;
                int neighbor_depth = 1;

                // Vérifie si le voisin se trouve à l'intérieur du cercle de rayon max_seuil et de centre s
                if (use_depth_limit && neighbor_depth > max_depth)
                    continue;

                if (cell_depth[index] < 0 || cell_depth[index] > neighbor_depth)
                    cell_depth[index] = neighbor_depth;
                if (source_tag[index] < 0)
                    source_tag[index] = s;

                if (!list_contains(narrow, index)){
                    list_push_back(narrow, index);
                }
            }
        }
    }


    // Boucle principale
    while(!list_is_empty(narrow)) {
        // Enlève le premier élément de la liste
        int index = list_pop_front(narrow);
        if (index < 0)
            continue;

        int i = index / m;
        int j = index % m;

        if (use_depth_limit && (cell_depth[index] < 0 || cell_depth[index] > max_depth))
            continue;

        double T_old = g->T[index];
        g->T[index] = eikonal_solve_local(g, i, j);
        double diff = fabs(g->T[index] - T_old);
        source_tag[index] = find_tag(g, source_tag, i, j, n, m);


        if (diff <= epsilon) {
            int current_cell_depth = cell_depth[index];
            if (use_depth_limit && current_cell_depth >= max_depth)
                continue;

            // La maille a convergé: on la fige et ses voisins susceptibles d'être améliorés sont ajoutés à la liste
            int neighbors[4][2] = {{i-1, j}, {i+1, j}, {i, j-1}, {i, j+1}};
            for (int k=0; k < 4; k++){
                int ni = neighbors[k][0];
                int nj = neighbors[k][1];
                if (ni >= 0 && ni < n && nj >= 0 && nj < m){
                    int index_neighbor = ni * m + nj;
                    int neighbor_depth = current_cell_depth + 1;

                    // Vérifie le rayon
                    if (use_depth_limit && neighbor_depth > max_depth)
                        continue;

                    // Vérifie si le voisin peut être amélioré
                    double T_neighbor_new = eikonal_solve_local(g, ni, nj);
                    if (T_neighbor_new < g->T[index_neighbor]){
                        g->T[index_neighbor] = T_neighbor_new;
                        source_tag[index_neighbor] = find_tag(g, source_tag, ni, nj, n, m);
                        if (source_tag[index_neighbor] < 0)
                            source_tag[index_neighbor] = source_tag[index];
                        if (cell_depth[index_neighbor] < 0 || cell_depth[index_neighbor] > neighbor_depth)
                            cell_depth[index_neighbor] = neighbor_depth;
                        if (!list_contains(narrow, index_neighbor))
                            list_push_back(narrow, index_neighbor);
                    }
                }
            }
        }
        else
            list_push_front(narrow,index);
    }

    // Enregistrement des tags
    eikonal_save_tags(g, source_tag, "source_tags.txt");

    // Nettoyage
    list_free(narrow);
    free(source_tag);
    free(cell_depth);
}