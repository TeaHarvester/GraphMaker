#include<iostream>
#include<cmath>
#include<random>
#include "graph.h"

void Graph::GenerateLFRGraph(unsigned int dim,
                             unsigned int k_min,
                             unsigned int k_max,
                             float a,
                             unsigned int n_comm,
                             float b,
                             float mixing_parameter)
{
    // generate graph with random degree for each vertex based on a power law
    // 
    // dim: number of vertices
    // k_min: minimum vertex degree
    // k_max: maximum vertex degree
    // a: power that alters the distribution of node degrees
    // n_comm: number of true communities
    // b: power that alters the distribution of community sizes
    // mixing_parameter: controls the modularity/amount of noise in the graph
    //
    // k ~ Distr(k_min, k_max, a)
    // ideally:
    // pmf = x^(-a) / SUM_X[x^(-a)] 
    // support: x c [k_min, k_max]

    // work around the default uniform PMF
    // substitute new support with the approximation:
    // x_min = k_min^(-1/a)
    // x_max = k_max^(-1/a)

    // generate random number between x_min and x_max and raise by -a
    // tis a silly approximation

    dimension = dim;
    n_communities = n_comm;

    if (adjacency_matrix != NULL)
    {
        delete adjacency_matrix;
    }

    adjacency_matrix = new SparseMatrix(dimension);

    if (true_communities != NULL)
    {
        delete true_communities;
    }

    true_communities = new std::vector<unsigned int>(dimension);

    float x_min = pow(float(k_min),-1/a);
    float x_max = pow(float(k_max),-1/a);
    float c_min = 0;
    float c_max = pow(float(n_comm),1/b);
    std::random_device rd;
    std::uniform_real_distribution<> Ddist(x_min, x_max);
    std::uniform_real_distribution<> Cdist(c_min, c_max);
    std::vector<std::vector<unsigned int>> available_vertices(dimension, std::vector<unsigned int>(3));

    for (unsigned int i = 0; i < dimension; ++i)
    {
        int k = (int)pow(Ddist(rd), -a);
        int c = (int)pow(Cdist(rd), b);
        available_vertices[i][0] = i;
        available_vertices[i][1] = k * mixing_parameter;
        available_vertices[i][2] = k * (1 - mixing_parameter);
        (*true_communities)[i] = c;
    }

    srand(2948372);
    unsigned int n_iters = 0;

    while (available_vertices.size() > 1)
    {
        unsigned int index_1 = 0;
        unsigned int index_2 = 0;

        while (index_1 == index_2)
        {
            index_1 = rand() % available_vertices.size();
            index_2 = rand() % available_vertices.size();
        } 

        unsigned int vertex_1 = available_vertices[index_1][0];
        unsigned int vertex_2 = available_vertices[index_2][0];
        unsigned int valent_index = 2;

        if (n_iters > dimension)
        {
            for (unsigned int i = 0; i < available_vertices[index_1][1] + available_vertices[index_1][2]; ++i)
            {   
                start:
                unsigned int n_iters_2 = 0;

                while (vertex_1 == vertex_2 || adjacency_matrix->IsAdjacent(vertex_1, vertex_2))
                {
                    if (n_iters_2 > dimension)
                    {
                        std::cout << "Vertex " << vertex_1 << " is saturated; skipping" << std::endl;
                        ++i;
                        goto start; 
                    }

                    vertex_2 = rand() % dimension;
                    ++n_iters_2;
                } 

                std::cout << "no suitable matches; connecting vertices " << vertex_1 << " and " << vertex_2 << std::endl;

                adjacency_matrix->AddConnection(vertex_1, vertex_2, 1);
                adjacency_matrix->AddConnection(vertex_2, vertex_1, 1);
                available_vertices.erase(available_vertices.cbegin() + index_1);
                n_iters_2 = 0;
                continue;
            }
        }

        if ((*true_communities)[vertex_1] == (*true_communities)[vertex_2])
        {
            valent_index = 1;
        }

        if (available_vertices[index_1][valent_index] > 0 &&
            available_vertices[index_1][valent_index] > 0 &&
            !adjacency_matrix->IsAdjacent(vertex_1, vertex_2))
        {
            adjacency_matrix->AddConnection(vertex_1, vertex_2, 1);
            adjacency_matrix->AddConnection(vertex_2, vertex_1, 1);

            --available_vertices[index_1][valent_index];
            --available_vertices[index_2][valent_index];

            if (available_vertices[index_1][1] + available_vertices[index_1][2] == 0)
            {
                available_vertices.erase(available_vertices.cbegin() + index_1);

                if (index_2 > index_1)
                {
                    --index_2;
                } 
            }

            if (available_vertices[index_2][1] + available_vertices[index_2][2] == 0)
            {
                available_vertices.erase(available_vertices.cbegin() + index_2);
            }

            n_iters = 0;
            continue;
        }

        ++n_iters;
    }

    n_edges = adjacency_matrix->GetEdges().size();
}

Graph::Graph()
:
dimension(0),
n_edges(0),
n_communities(0),
adjacency_matrix(NULL),
true_communities(NULL)
{}

Graph::~Graph()
{
    if (adjacency_matrix != NULL)
    {
        delete adjacency_matrix;
    }

    if (true_communities != NULL)
    {
        delete true_communities;
    }
}