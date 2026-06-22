#include <stdio.h>
#include <stdlib.h>

#include "../include/SolveEikonal2D.h"

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




int main(){
    printf("\nTest double-linked list\n");

    NodeList* l = list_create(100);
    if (!l){
        printf("list_create failed\n");
        return 1;
    }
    printf("list_create OK\n");

    list_push_back(l, 42);
    if (l->size == 1 && list_contains(l, 42))
        printf("list_push_back OK\n");
    else
        ("list_push_back failed\n");

    int index = list_pop_front(l);
    if (index == 42 && l->size == 0)
        printf("list_pop_front OK\n");
    else
        printf("lists_pop_front failed\n");

    list_free(l);
    printf("list_free OK\n\n");
    return 0;
}
