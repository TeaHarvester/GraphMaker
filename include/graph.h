#ifndef GRAPH
#define GRAPH
#include<vector>
#include "sparsematrix.h"

class Graph
{
    public: unsigned int dimension;
    public: SparseMatrix *adjacency_matrix;
    std::vector<unsigned int> *true_communities;

    public: void GenerateLFRGraph(int dim,
                                 int k_min,
                                 int k_max,
                                 float power1,
                                 int n_comm,
                                 float power2,
                                 float mixing_parameter
                                );

    Graph();
    ~Graph();
};

#endif