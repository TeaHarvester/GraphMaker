#ifndef GRAPH
#define GRAPH
#include "sparsematrix.h"

class Graph
{
    public: unsigned int dimension;
    public: SparseMatrix *adjacency_matrix;

    public: void GenerateLFRGraph(int dim,
                                 int k_min,
                                 int k_max,
                                 float power1
                                //  int n_comm,
                                //  int s_min,
                                //  float power2,
                                //  int seed,
                                //  float noise
                                );

    Graph();
    ~Graph();
};

#endif