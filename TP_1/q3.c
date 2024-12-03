#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*
    STRING PROCESSING ALGS
*/

// Function to compute the maximum overlap using a naive method
int compute_max_overlap_naive(const char *u, const char *v) {
    int max_overlap = 0;
    int u_len = strlen(u);
    int v_len = strlen(v);

    for (int i = 1; i <= u_len && i <= v_len; ++i) {
        if (strncmp(u + u_len - i, v, i) == 0) {
            max_overlap = i;
        }
    }
    return max_overlap;
}

/**
 * @brief Computes the maximum overlap between two strings.
 *
 * This function calculates the maximum overlap between the suffix of the first string (u)
 * and the prefix of the second string (v). It uses the Knuth-Morris-Pratt (KMP) algorithm
 * to find the longest prefix which is also a suffix in the concatenated string "v#u".
 *
 * @param u The first input string.
 * @param v The second input string.
 * @return The length of the maximum overlap between the suffix of u and the prefix of v.
 */
int compute_max_overlap(const char *u, const char *v) {
    int u_len = strlen(u);
    int v_len = strlen(v);
    int string_len = v_len + u_len + 1;

    // Concatenate v, separator (#), and u
    char *string = (char *)malloc((string_len + 1) * sizeof(char));
    sprintf(string, "%s#%s", v, u);

    int *lps = (int *)calloc(string_len, sizeof(int));
    int length = 0;

    for (int i = 1; i < string_len; ++i) {
        while (length > 0 && string[i] != string[length]) {
            length = lps[length - 1];
        }
        if (string[i] == string[length]) {
            ++length;
        }
        lps[i] = length;
    }

    int max_overlap = lps[string_len - 1];

    free(lps);
    free(string);

    return max_overlap;
}

/*
    GRAPH RELATED WHAT-NOTS
*/

typedef struct AdjListNode {
    int dest;
    int overlap;
    struct AdjListNode* next;
} AdjListNode;

typedef struct AdjList {
    AdjListNode* head;
} AdjList;

typedef struct Graph {
    int num_fragments;
    AdjList* array;
} Graph;

AdjListNode* newAdjListNode(int dest, int overlap) {
    AdjListNode* newNode = (AdjListNode*)malloc(sizeof(AdjListNode));
    newNode->dest = dest;
    newNode->overlap = overlap;
    newNode->next = NULL;
    return newNode;
}

// Function to create a graph of given number of fragments
Graph* createGraph(int num_fragments) {
    Graph* graph = (Graph*)malloc(sizeof(Graph));
    graph->num_fragments = num_fragments;
    graph->array = (AdjList*)malloc(num_fragments * sizeof(AdjList));
    for (int i = 0; i < num_fragments; ++i) {
        graph->array[i].head = NULL;
    }
    return graph;
}

// Function to add an edge to the graph
void addEdge(Graph* graph, int src, int dest, int overlap) {
    AdjListNode* newNode = newAdjListNode(dest, overlap);
    newNode->next = graph->array[src].head;
    graph->array[src].head = newNode;
}

/**
 * @brief Prints the graph structure along with its associated fragments.
 *
 * This function takes a graph and an array of string fragments, and prints
 * the graph structure along with the corresponding fragments. The number of
 * fragments is also provided as an input.
 *
 * @param graph Pointer to the Graph structure to be printed.
 * @param fragments Array of string fragments associated with the graph.
 * @param num_fragments Number of string fragments in the fragments array.
 */
void print_graph(Graph* graph, char **fragments, int num_fragments) {
    for (int i = 0; i < num_fragments; ++i) {
        AdjListNode* pCrawl = graph->array[i].head;
        printf("Fragment %d (%s) overlaps:\n", i, fragments[i]);
        while (pCrawl) {
            printf(" -> Fragment %d (%s) with overlap %d\n", pCrawl->dest, fragments[pCrawl->dest], pCrawl->overlap);
            pCrawl = pCrawl->next;
        }
        printf("\n");
    }
}

/**
 * @brief Builds an overlap graph from a set of string fragments.
 *
 * This function constructs an overlap graph where each node represents a string fragment
 * and edges represent the overlap between fragments as determined by the provided overlap function.
 *
 * @param fragments An array of string fragments to be used for building the overlap graph.
 * @param num_fragments The number of string fragments in the array.
 * @param overlap_function A pointer to a function that takes two string fragments as input
 *                         and returns an integer representing the overlap between them.
 * @return A pointer to the constructed Graph.
 */
Graph* build_overlap_graph(char **fragments, int num_fragments, int (*overlap_function)(const char *, const char *)) {
    Graph* graph = createGraph(num_fragments);
    for (int i = 0; i < num_fragments; ++i) {
        for (int j = 0; j < num_fragments; ++j) {
            if (i != j) {
                int overlap = overlap_function(fragments[i], fragments[j]);
                addEdge(graph, i, j, overlap);
            }
        }
    }
    return graph;
}

int main() {
    // Example usage
    char *fragments[] = {"TAG", "CTA", "ACT"};
    int num_fragments = sizeof(fragments) / sizeof(fragments[0]);
    Graph* g = build_overlap_graph(fragments, num_fragments, compute_max_overlap);
    print_graph(g, fragments, num_fragments);
    return 0;
}