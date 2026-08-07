#ifndef FIM2D_H
#define FIM2D_H

#include "SolveEikonal2D.h"

typedef struct Node{
    int index;  // Indice i * m + j
    struct Node* prev;
    struct Node* next;
}Node;

typedef struct NodeList{
    Node head;
    Node** index_to_node;
    int ncell;  // ncell = n*m
    int size;   // Nombre d'éléments dans la liste
}NodeList;

// Création de la liste
NodeList* list_create(int ncell);

// Libération de la mémoire de la liste
void list_free(NodeList* l);

// Ajout d'un élément avant le noeud 'node'
void list_insert_before(NodeList* l, int index, Node* node);

// Ajout d'un élément au début de la liste
void list_push_front(NodeList* l, int index);

// Ajout d'un élément à la fin de la liste
void list_push_back(NodeList* l, int index);

// Supprime le noeud qui contient l'indice 'index' de la liste
void list_remove(NodeList* l, int index);

// Récupère et renvoie le premier élément de la liste
int list_pop_front(NodeList* l);

// Teste si la liste possède un noeud qui contient l'indice 'index'
int list_contains(const NodeList* l, int index);

// Test si la liste est vide
int list_is_empty(const NodeList* l);

// Vérifie si un voisin existe et peut être amélioré
int is_neighbor_updatable(const EikonalGrid* g, int index, double T_new, double epsilon);

// Détermine le tag de la source qui a mis à jour la valeur T en (i,j)
int find_tag(const EikonalGrid* g, const int* source_tag, int i, int j, int n, int m);




// Résout l'équation eikonale sur la grille g à l'aide de la méthode FIM
// src_i est le tableau des indices i des sources
// src_j est le tableau des indices j des sources
// ns est le nombre de sources
// Si max_depth est >= 0, on limite la propagation à max_depth voisinages 4-connexes depuis les sources.
void fim_solve(EikonalGrid* g, const int* src_i, const int* src_j, int ns, double epsilon, int max_depth);


#endif /* FIM2D_H */