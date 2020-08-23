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
    float GetMixingParameter(bool ground_truth);
    float GetMixingParameter(bool ground_truth, unsigned int community);
    float GetMixingParameter(bool ground_truth, unsigned int community_1, unsigned int community_2);
    std::vector<unsigned int> GetCommunityMembers(bool ground_truth, unsigned int community);

    private:
    void OptimiseLFRGraph(float target_mixing_parameter);
    bool LouvainOptimise(Graph& tmp_graph);
    float LouvainGetModularity(Graph&tmp_graph, unsigned int vertex_1, unsigned int community);
    Graph* LouvainAggregate(Graph& tmp_graph);
    void GetDegree();

    public:
    Graph();
    Graph(const SparseMatrix& S);
    Graph(const Graph& G);
    ~Graph();
};

#endif