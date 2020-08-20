#ifndef GRAPH
#define GRAPH
#include<vector>
#include "sparsematrix.h"

class Graph
{
    public: 
    SparseMatrix *adjacency_matrix;
    unsigned int dimension;
    unsigned int n_edges;
    unsigned int n_communities;
    unsigned int n_comm_detected;
    unsigned int max_degree;
    std::vector<unsigned int>* degree;
    std::vector<unsigned int> *true_communities;
    std::vector<unsigned int> *detected_communities;

    void GenerateLFRGraph(unsigned int dim,
                          unsigned int k_min,
                          unsigned int k_max,
                          float power1,
                          unsigned int n_comm,
                          float power2,
                          float mixing_parameter);

    // detect communities using Louvain algorithm
    void Louvain(Graph& G, unsigned int recur_counter);
    private:
    bool LouvainOptimise(Graph& G);
    float LouvainGetModularity(Graph &G, unsigned int vertex_1, unsigned int community);
    Graph LouvainAggregate(Graph& G);
    void GetDegree();

    public:
    Graph();
    Graph(SparseMatrix*& M);
    ~Graph();
};

#endif