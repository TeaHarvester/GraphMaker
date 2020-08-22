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

    if (detected_communities != NULL)
    {
        delete detected_communities;
    }

    detected_communities = new std::vector<unsigned int>(dimension); 

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

        unsigned int valent_index = 2;

        if (n_iters > dimension)
        {
            unsigned int n_iters_2 = 0;

            for (unsigned int i = 0; i < available_vertices[index_1][1] + available_vertices[index_1][2]; ++i)
            {   
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

    detected_communities = new std::vector<unsigned int>(dimension);
	n_comm_detected = dimension;

    for (unsigned int i = 0; i < dimension; ++i)
    {
        (*detected_communities)[i] = i;
    }

    GetDegree();
}

void Graph::Louvain(Graph& G, unsigned int recur_counter)
{   
    Graph* tmp_graph = new Graph(G);

    if (recur_counter > 0)
    {
        delete &G;
    }

    else
    {
        std::cout << std::endl << "running Louvain community detection: " << std::endl;
    }
    

    bool optimal_partition = false;
	unsigned int optim_iter = 0;

	while (optimal_partition == false)
	{
        std::cout << "partitioning... " << std::endl;
		optimal_partition = LouvainOptimise(*tmp_graph);
		++optim_iter;
	}
    
    for (unsigned int i = 0; i < this->dimension; ++i)
    {
        for (unsigned int j = 0; j < tmp_graph->dimension; ++j)
        {
            if ((*(this->detected_communities))[i] == j)
            {
                (*(this->detected_communities))[i] = (*(tmp_graph->detected_communities))[j];

                break;
            }
        }
    }

    n_comm_detected = tmp_graph->n_comm_detected;

    std::cout << std::endl << "aggregating... " << std::endl;
    Graph* aggregate_graph = LouvainAggregate(*tmp_graph);

    delete &tmp_graph;

	if (optim_iter > 1)
	{
		Louvain(*aggregate_graph, recur_counter + 1);
	}

	else
	{
		std::cout << std::endl << "number of communities detected: " << n_comm_detected << std::endl;
	}
}

bool Graph::LouvainOptimise(Graph& tmp_graph)
{
	std::vector<unsigned int> community_sizes(tmp_graph.n_comm_detected);
	bool maximum_modularity = true;
	// Q: is modularity;
    // delta_Q: difference in Q after assigning a node to a community
	float delta_Q;
	float delta_Q_max;
	unsigned int optimal_community;

	// Records the size of each community
	for (unsigned int i = 0; i < tmp_graph.n_comm_detected; ++i)
	{
		unsigned int iter = 0;

		for (unsigned int j = 0; j < tmp_graph.dimension; ++j)
		{
			if ((*tmp_graph.detected_communities)[j] == i)
			{
				++iter;
        	}
			
			community_sizes[i] = iter;
		}
	}

	for (unsigned int i = 0; i < tmp_graph.dimension; ++i)
	{
		delta_Q_max = 0;

		for (unsigned int j = 0; j < tmp_graph.n_comm_detected; ++j)
		{
			delta_Q = LouvainGetModularity(tmp_graph, i, j);

			if (delta_Q > delta_Q_max)
			{
				delta_Q_max = delta_Q;
				optimal_community = j;
			}
		}

		if (delta_Q_max > 0)
		{
			unsigned int drained_community = (*tmp_graph.detected_communities)[i];

			(*tmp_graph.detected_communities)[i] = optimal_community;
			++community_sizes[optimal_community];
			--community_sizes[drained_community];

			if (community_sizes[drained_community] == 0)
			{
				community_sizes.erase(community_sizes.cbegin() + drained_community);
				--tmp_graph.n_comm_detected;

				for (unsigned int j = 0; j < tmp_graph.dimension; ++j)
				{
					if ((*tmp_graph.detected_communities)[j] > drained_community)
					{
						--(*tmp_graph.detected_communities)[j];
					}
				}

				maximum_modularity = false;
            }
		}
	}

	return maximum_modularity;
}

float Graph::LouvainGetModularity(Graph& tmp_graph, unsigned int vertex_1, unsigned int community)
{
    // Change in modularity moving vertex_i into community c
	// delta_Q = (k_i_c / two_m) - (sigma_total * k_i / 2 * m ^ 2)

	float m = 0;
	float k_i_c = 0;
	float k_i = (float)(*tmp_graph.degree)[vertex_1];
	float sigma_in = 0;
	float sigma_total = 0;
	float sigma_total_in = 0;
	float sigma_total_out = 0;
	int c_size = 0;

	for (unsigned int i = 0; i < tmp_graph.dimension; ++i)
	{
		if ((*tmp_graph.detected_communities)[i] == community)
		{
			++c_size;
			k_i_c += tmp_graph.adjacency_matrix->GetEdgeWeight(vertex_1, i);
		}

		for (unsigned int j = 0; j < tmp_graph.dimension; ++j)
		{
            float edge_weight = tmp_graph.adjacency_matrix->GetEdgeWeight(i, j);

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

			if ((*tmp_graph.detected_communities)[i] == community)
			{
				sigma_total += edge_weight;
				sigma_total_out += edge_weight;
				sigma_total_in += edge_weight;

				if ((*tmp_graph.detected_communities)[j] == community)
				{
					sigma_in += edge_weight;
				}
			}
		}
	}

	return (k_i_c / (2.0f * m)) - ((sigma_total * k_i) / (2.0f * pow(m, 2.0f)));
}

Graph* Graph::LouvainAggregate(Graph& tmp_graph)
{
    unsigned int new_dimension = tmp_graph.n_comm_detected;
	SparseMatrix aggregate(new_dimension);

	for (unsigned int i = 0; i < tmp_graph.dimension; ++i)
	{
		for (unsigned int j = 0; j < tmp_graph.dimension; ++j)
		{
            unsigned int vertex_1 = (*tmp_graph.detected_communities)[i];
            unsigned int vertex_2 = (*tmp_graph.detected_communities)[j];
            aggregate.AddConnection(vertex_1, vertex_2, tmp_graph.adjacency_matrix->GetEdgeWeight(i, j));
		}
	}

	return new Graph(aggregate);
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

Graph::Graph(const SparseMatrix& S)
:
dimension(S.dimension),
n_edges(S.col.size()),
n_communities(0),
n_comm_detected(S.dimension),
max_degree(0),
degree(NULL),
true_communities(NULL),
detected_communities(NULL)
{
    adjacency_matrix = new SparseMatrix(S);
	n_comm_detected = dimension;
    detected_communities = new std::vector<unsigned int>(dimension);

    for (unsigned int i = 0; i < dimension; ++i)
    {
        (*detected_communities)[i] = i;
    }

    GetDegree();
}

Graph::Graph(const Graph& G)
:
dimension(G.dimension),
n_edges(G.adjacency_matrix->col.size()),
n_communities(G.n_communities),
n_comm_detected(G.n_comm_detected),
max_degree(0),
degree(NULL),
true_communities(NULL),
detected_communities(NULL)
{
    adjacency_matrix = new SparseMatrix(*G.adjacency_matrix);

    if (G.true_communities != NULL)
    {
        true_communities = new std::vector<unsigned int>(*G.true_communities);
    }

    if (G.detected_communities != NULL)
    {
        detected_communities = new std::vector<unsigned int>(*G.detected_communities);
    }

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