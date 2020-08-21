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
        available_vertices[i][1] = (unsigned int)((float)k * mixing_parameter);
        available_vertices[i][2] = k - available_vertices[i][1];
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
        unsigned int& av_1 = available_vertices[index_1][0];
        unsigned int& av_2 = available_vertices[index_1][1];
        unsigned int& av_3 = available_vertices[index_1][2];

        unsigned int valent_index = 2;

        if (n_iters > dimension)
        {
            unsigned int n_iters_2 = 0;

            for (unsigned int i = 0; i < available_vertices[index_1][1] + available_vertices[index_1][2]; ++i)
            {   
                std::cout << "i: " << i << std::endl;
                while (vertex_1 == vertex_2 || adjacency_matrix->IsAdjacent(vertex_1, vertex_2))
                {
                    if (n_iters_2 > dimension)
                    {
                        std::cout << "Vertex " << vertex_1 << " is saturated; skipping" << std::endl;
                        goto end;
                    }

                    vertex_2 = rand() % dimension;
                    ++n_iters_2;
                } 

                std::cout << "no suitable matches; connecting vertices " << vertex_1 << " and " << vertex_2 << std::endl;

                adjacency_matrix->AddConnection(vertex_1, vertex_2, 1);
                adjacency_matrix->AddConnection(vertex_2, vertex_1, 1);
            }
                end:
                available_vertices.erase(available_vertices.cbegin() + index_1);
                n_iters_2 = 0;
                continue;
        }

        if ((*true_communities)[vertex_1] == (*true_communities)[vertex_2])
        {
            valent_index = 1;
        }

        if (available_vertices[index_1][valent_index] > 0 &&
            available_vertices[index_2][valent_index] > 0 &&
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
    GetDegree();
}

void Graph::Louvain(Graph& G, unsigned int recur_counter)
{
	if (recur_counter == 0)
	{
		if (detected_communities != NULL)
		{
			delete detected_communities;
		}

		detected_communities = new std::vector<unsigned int>(dimension);
		n_comm_detected = dimension;

		for (unsigned int i = 0; i < dimension; ++i)
		{
			(*detected_communities)[i] = i;
		}
	}

    bool optimal_partition = false;
	unsigned int optim_iter = 0;
	
	while (optimal_partition == false)
	{
		optimal_partition = LouvainOptimise(G);
		++optim_iter;
		std::cout << std::endl << "partitioning... ";
	}

	for (unsigned int i = 0; i < this->dimension; ++i)
	{
		for (unsigned int j = 0; j < G.dimension; ++j)
		{
			if ((*this->detected_communities)[i] == j)
			{
				this->detected_communities[i] = G.detected_communities[j];
				break;
			}
		}
	}

    this->n_comm_detected = G.n_comm_detected;

    Graph aggregate_graph = LouvainAggregate(G);
    
	std::cout << std::endl << std::endl << "aggregating... ";

	if (recur_counter > 0)
	{
        G.~Graph();
	}

	if (optim_iter > 1)
	{
		Louvain(aggregate_graph, recur_counter + 1);
	}

	else
	{
		std::cout << std::endl << "number of communities detected: " << n_comm_detected;
	}
}

bool Graph::LouvainOptimise(Graph& G)
{
	std::vector<unsigned int> community_sizes(G.n_comm_detected);
	bool maximum_modularity = true;
	// Q: is modularity;
    // delta_Q: difference in Q after assigning a node to a community
	float delta_Q;
	float delta_Q_max;
	unsigned int optimal_community;

	// Records the size of each community
	for (unsigned int i = 0; i < G.n_comm_detected; ++i)
	{
		unsigned int iter = 0;

		for (unsigned int j = 0; j < G.dimension; ++j)
		{
			if ((*G.detected_communities)[j] == i)
			{
				++iter;
			}
			
			community_sizes[i] = iter;
		}
	}

	for (unsigned int i = 0; i < G.dimension; ++i)
	{
		delta_Q_max = 0;

		for (unsigned int j = 0; j < G.n_comm_detected; ++j)
		{
			delta_Q = LouvainGetModularity(G, i, j);

			if (delta_Q > delta_Q_max)
			{
				delta_Q_max = delta_Q;
				optimal_community = j;
			}
		}

		if (delta_Q_max > 0)
		{
			unsigned int drained_community = (*G.detected_communities)[i];

			(*G.detected_communities)[i] = optimal_community;
			++community_sizes[optimal_community];
			--community_sizes[drained_community];

			if (community_sizes[drained_community] == 0)
			{
				community_sizes.erase(community_sizes.cbegin() + drained_community);
				--G.n_comm_detected;

				for (unsigned int j = 0; j < G.dimension; ++j)
				{
					if ((*G.detected_communities)[j] > drained_community)
					{
						--(*G.detected_communities)[j];
					}
				}
				maximum_modularity = false;
			}
		}
	}

	return maximum_modularity;
}

float Graph::LouvainGetModularity(Graph &G, unsigned int vertex_1, unsigned int community)
{
    // Change in modularity moving vertex_i into community c
	// delta_Q = (k_i_c / two_m) - (sigma_total * k_i / 2 * m ^ 2)

	float m = 0;
	float k_i_c = 0;
	float k_i = (float)(*G.degree)[vertex_1];
	float sigma_in = 0;
	float sigma_total = 0;
	float sigma_total_in = 0;
	float sigma_total_out = 0;
	int c_size = 0;

	for (unsigned int i = 0; i < G.dimension; ++i)
	{
		if ((*G.detected_communities)[i] == community)
		{
			++c_size;
			k_i_c += G.adjacency_matrix->GetEdgeWeight(vertex_1, i);
		}

		for (unsigned int j = 0; j < G.dimension; ++j)
		{
            float edge_weight = G.adjacency_matrix->GetEdgeWeight(i, j);

			if (j < i)
			{
				continue;
			}

			if (j == i)
			{
				m += edge_weight / 2;
			}

			else
			{
				m += edge_weight;
			}

			if ((*G.detected_communities)[i] == community)
			{
				sigma_total += edge_weight;
				sigma_total_out += edge_weight;
				sigma_total_in += edge_weight;

				if ((*G.detected_communities)[j] == community)
				{
					sigma_in += edge_weight;
				}
			}
		}
	}

	return (k_i_c / (2.0f * m)) - ((sigma_total * k_i) / (2.0f * pow(m, 2.0f)));
}

Graph Graph::LouvainAggregate(Graph& G)
{
    unsigned int new_dimension = G.n_comm_detected;
	SparseMatrix aggregate(new_dimension);

	for (unsigned int i = 0; i < G.dimension; ++i)
	{
		for (unsigned int j = 0; j < G.dimension; ++j)
		{
			//w_new[communities[i]][communities[j]] += w[i][j];
            unsigned int vertex_1 = (*G.detected_communities)[i];
            unsigned int vertex_2 = (*G.detected_communities)[j];
            aggregate.AddConnection(vertex_1, vertex_2, G.adjacency_matrix->GetEdgeWeight(i, j));
		}
	}

	for (unsigned int i = 0; i < new_dimension; ++i)
	{
		(*G.detected_communities)[i] = i;
        G.detected_communities->pop_back();
	}

    SparseMatrix* S;
    S = &aggregate;

	return Graph(S);
}

void Graph::GetDegree()
{
    if (degree != NULL)
    {
        delete degree;
    }

    degree = new std::vector<unsigned int>(dimension, 0.0f);

    for (unsigned int i = 0; i < dimension; ++i)
    {
        (*degree)[i] = adjacency_matrix->GetDegree(i);
        if ((*degree)[i] > max_degree)
        {
            max_degree = (*degree)[i];
        }
    }
}

Graph::Graph()
:
dimension(0),
n_edges(0),
n_communities(0),
n_comm_detected(0),
max_degree(0),
degree(NULL),
true_communities(NULL),
detected_communities(NULL)
{
    adjacency_matrix = NULL;
}

Graph::Graph(SparseMatrix*& M)
:
adjacency_matrix(M),
dimension(M->dimension),
n_edges(M->col.size()),
n_communities(0),
n_comm_detected(0),
max_degree(0),
degree(NULL),
true_communities(NULL),
detected_communities(NULL)
{
    GetDegree();
}

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

    if (degree != NULL)
    {
        delete true_communities;
    }
}