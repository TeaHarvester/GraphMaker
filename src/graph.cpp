#include<iostream>
#include<cmath>
#include<random>
#include "graph.h"

void Graph::GenerateLFRGraph(int dim,
                             int k_min,
                             int k_max,
                             float a
                            //  int n_comm,
                            //  int s_min,
                            //  float b,
                            //  int seed,
                            //  float noise
                            )
{
    // generate graph with random degree for each vertex based on a power law
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

    if (adjacency_matrix != NULL)
    {
        delete adjacency_matrix;
    }

    //SparseMatrix new_adjacency_matrix(dimension);
    adjacency_matrix = new SparseMatrix(dimension);

    float x_min = pow(float(k_min),-1/a);
    float x_max = pow(float(k_max),-1/a);
    std::random_device rd;
    std::uniform_real_distribution<> dist(x_min, x_max);

    std::vector<unsigned int> valency(dim);

    for (int i = 0; i < dimension; ++i)
    {
        int k = (int)pow(dist(rd), -a);
        valency[i] = k;
    }

    srand(2948372);

    for (int i = 0; i < dimension; ++i)
    {
        int iters = 0;

        while (valency[i] > 0)
        {
            int paired_vertex = rand() % dimension;

            if (paired_vertex != i)
            {
                if (valency[paired_vertex] > 0 || iters > dimension)
                {
                    adjacency_matrix->AddConnection(i, paired_vertex, 1);
                    adjacency_matrix->AddConnection(paired_vertex, i, 1);
                    --valency[i];
                    iters = 0;

                    if (valency[paired_vertex] == 0)
                    {
                        std::cout << "Couldn't find a partner for vertex " << i << ", pairing regardless" << std::endl; 
                        continue;  
                    }

                    std::cout << "pairing " << i << "with " << paired_vertex << std::endl;

                    --valency[paired_vertex];
                }
            }     

            ++iters;
        }
    }

    return;
}

Graph::Graph()
:
dimension(0),
adjacency_matrix(NULL)
{}

Graph::~Graph()
{
    if (adjacency_matrix != NULL)
    {
        delete adjacency_matrix;
    }
}