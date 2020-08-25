#include<iostream>
#include<cmath>
#include<random>
#include <fstream>
#include <sstream>
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
        available_vertices[i][1] = (unsigned int)((float)k * (1 - mixing_parameter));
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

    while (!OptimiseLFRGraph(mixing_parameter))
    {
        std::cout << "optimising LFR graph..." << std::endl;
    }
}

bool Graph::OptimiseLFRGraph(float target_mixing_parameter)
{
    std::vector<float> community_mix_params(n_communities);
    float initial_mixing_parameter =  GetMixingParameter(true);
    float current_mixing_parameter = initial_mixing_parameter;
    unsigned int community_iterator = 0;

    while(community_iterator < n_communities)
    {
        start:
        float deviation = target_mixing_parameter - current_mixing_parameter;
        float sign_coeff = floor(deviation) == -1.0f ? -1.0f : 1.0f;

        // select vertices to disconnect/reconnect

        std::vector<unsigned int> stray_vertices = GetCommunityMembers(true, community_iterator);
        unsigned int s_vertex_1;
        unsigned int s_vertex_2;
        unsigned int p_vertex;

        std::random_device rd;
        std::uniform_int_distribution<> Sdist(0, stray_vertices.size() - 1);
        bool s_chosen = false;
        bool p_chosen = false;
        unsigned int s_iters = 0;
        unsigned int p_iters = 0;

        while(!s_chosen)
        {
            s_chosen = false;
            p_chosen = false;

            if (s_iters > std::pow((float)dimension, 2.0f))
            {
                std::cout << "community " << community_iterator << " is saturated" << std::endl << std::endl;
                ++community_iterator;

                if (community_iterator < n_communities)
                {   
                    goto start;
                }

                else 
                {
                    goto end;
                }
            }

            s_vertex_1 = stray_vertices[Sdist(rd)];
            s_vertex_2 = stray_vertices[Sdist(rd)];

            if (s_vertex_1 == s_vertex_2)
            {
                continue;
            }

            if (sign_coeff == -1 && !adjacency_matrix->IsAdjacent(s_vertex_1, s_vertex_2))      
            {

                s_chosen = true;
                std::vector<unsigned int> adjacent_1 = adjacency_matrix->GetAdjacentVertices(s_vertex_1);

                for (unsigned int i = 0; i < adjacent_1.size(); ++i)
                {
                    p_vertex = adjacent_1[i];

                    if ((*true_communities)[p_vertex] != community_iterator &&
                        (*degree)[p_vertex] > 1)
                    {
                        p_chosen = true;
                        s_iters = 0;
                    }
                }
            }

            else if (sign_coeff == 1 && 
                     adjacency_matrix->IsAdjacent(s_vertex_1, s_vertex_2) &&
                     (*degree)[s_vertex_1] > 1 &&
                     (*degree)[s_vertex_2] > 1)
            {
                s_chosen = true;

                std::uniform_int_distribution<> Pdist(0, dimension - 1);

                while (p_chosen == false)
                {
                    p_chosen = false;

                    if (p_iters > dimension)
                    {
                        s_chosen = false;
                        p_iters = 0;
                        break;
                    }

                    p_vertex = Pdist(rd);

                    if (!adjacency_matrix->IsAdjacent(s_vertex_1, p_vertex) &&
                        (*true_communities)[p_vertex] != community_iterator)
                    {
                        p_chosen = true;
                        s_iters = 0;
                        p_iters = 0;
                    }
                }
            }

            ++s_iters;
        }

        if (sign_coeff == -1)
        {
            adjacency_matrix->EraseConnection(s_vertex_1, p_vertex);
            adjacency_matrix->EraseConnection(p_vertex, s_vertex_1);
            adjacency_matrix->AddConnection(s_vertex_1, s_vertex_2, 1.0f);
            adjacency_matrix->AddConnection(s_vertex_2, s_vertex_1, 1.0f);
        }

        else
        {
            adjacency_matrix->EraseConnection(s_vertex_1, s_vertex_2);
            adjacency_matrix->EraseConnection(s_vertex_2, s_vertex_1);
            adjacency_matrix->AddConnection(s_vertex_1, p_vertex, 1.0f);
            adjacency_matrix->AddConnection(p_vertex, s_vertex_1, 1.0f);
        }

        std::cout << "rewiring..." << std::endl << std::endl;

        float new_mixing_parameter = GetMixingParameter(true);
        
        if (std::pow(target_mixing_parameter - new_mixing_parameter, 2.0f) >= std::pow(deviation, 2.0f))
        {
            if (sign_coeff == -1)
            {
                adjacency_matrix->EraseConnection(s_vertex_1, s_vertex_2);
                adjacency_matrix->EraseConnection(s_vertex_2, s_vertex_1);
                adjacency_matrix->AddConnection(s_vertex_1, p_vertex, 1.0f);
                adjacency_matrix->AddConnection(p_vertex, s_vertex_1, 1.0f);
            }

            else
            {
                adjacency_matrix->EraseConnection(s_vertex_1, p_vertex);
                adjacency_matrix->EraseConnection(p_vertex, s_vertex_1);
                adjacency_matrix->AddConnection(s_vertex_1, s_vertex_2, 1.0f);
                adjacency_matrix->AddConnection(s_vertex_2, s_vertex_1, 1.0f);
            }

            std::cout << "mixing parameter converged at " << new_mixing_parameter <<
             " for community " << community_iterator << std::endl << std::endl;

             ++community_iterator;
        }

        current_mixing_parameter = new_mixing_parameter;
        GetDegree();
    }

    end:

    if (std::pow(target_mixing_parameter - current_mixing_parameter, 2.0f) < 
        std::pow(target_mixing_parameter - initial_mixing_parameter, 2.0f))
        {
            return false;
        }

    return true;
}

float Graph::GetMixingParameter(bool ground_truth)
{
    const std::vector<unsigned int>& communities = ground_truth ? *true_communities : *detected_communities;

    float external_degree = 0.0f;
    float total_degree = 0.0f;

    for (unsigned int i = 0; i < dimension; ++i)
    {
        std::vector<unsigned int> adjacent_vertices = adjacency_matrix->GetAdjacentVertices(i);

        for (unsigned int j = 0; j < adjacent_vertices.size(); ++j)
        {
            float edge_weight = adjacency_matrix->GetEdgeWeight(i, adjacent_vertices[j]);

            if (communities[i] != communities[adjacent_vertices[j]])
            {
                external_degree += edge_weight;
            }

            total_degree += edge_weight;
        }
    }

    std::cout << "mixing parameter (" << (ground_truth ? "ground truth)" : "detected)") << " communities: "
    << external_degree / total_degree << std::endl;

    return external_degree / total_degree;
}

std::vector<unsigned int> Graph::GetCommunityMembers(bool ground_truth, unsigned int community)
{
    const std::vector<unsigned int>& communities = ground_truth ? *true_communities : *detected_communities;
    std::vector<unsigned int> members;

    for (unsigned int i = 0; i < dimension; ++i)
    {
        if (communities[i] == community)
        {
            members.push_back(i);
        }
    }

    return members;
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
        GetMixingParameter(false);
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

Graph::Graph(std::string file_path)
:
dimension(0),
n_edges(0),
n_communities(0),
max_degree(0),
degree(NULL),
true_communities(NULL),
detected_communities(NULL)
{
    std::ifstream i_stream;

	i_stream.open(file_path);

	if (!i_stream) 
    {
		std::cerr << "Unable to open file datafile";
		std::exit(1);   // call system to stop
	}

    std::string text_line;
	int line_counter = -1;

	while (std::getline(i_stream, text_line))
	{
		unsigned int char_counter = 0;
		unsigned int datum;
		std::string datum_str = "\0";

		// Extract words from lines
		while (char_counter < text_line.length())
		{
			// if character is not (space) or NULL, add character to 'word'
			if (text_line[char_counter] != 32)
			{
				datum_str.push_back(text_line[char_counter]);
				char_counter++;
			}
			
			else if (text_line[char_counter] == 32)
			{
                std::stringstream(datum_str) >> datum;

                if (line_counter == -1)
                {
                    dimension = datum;
                    adjacency_matrix = new SparseMatrix(dimension);
                    n_communities = dimension;
                    n_comm_detected = dimension;
                    break;
                }

				adjacency_matrix->AddConnection(line_counter, datum - 1, 1.0f);
				datum_str = "\0";
				++char_counter;
			}
		}

		++line_counter;
	}

    n_edges = adjacency_matrix->GetEdges().size();

    true_communities = new std::vector<unsigned int>(dimension);
    detected_communities = new std::vector<unsigned int>(dimension);

    for (unsigned int i = 0; i < dimension; ++i)
    {
        (*true_communities)[i] = i;
        (*detected_communities)[i] = i;
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