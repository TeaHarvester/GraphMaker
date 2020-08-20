#ifndef GRAPH
#define GRAPH
#include<vector>
#include "sparsematrix.h"

class Graph
{
    public: 
    unsigned int dimension;
    unsigned int n_edges;
    unsigned int n_communities;
    unsigned int max_degree;
    SparseMatrix *adjacency_matrix;
    std::vector<unsigned int>* degree;
    std::vector<unsigned int> *true_communities;

    void GenerateLFRGraph(unsigned int dim,
                          unsigned int k_min,
                          unsigned int k_max,
                          float power1,
                          unsigned int n_comm,
                          float power2,
                          float mixing_parameter);

    Graph();
    ~Graph();
};

#endif