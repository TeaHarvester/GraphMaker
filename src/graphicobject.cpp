#include<iostream>
#include<cmath>
#include<vector>
#include"graphicobject.h"

void GraphicObject::Init()
{
    n_vertices = source_graph->dimension;
    n_indices = source_graph->n_edges * 2;

    vertex_array = new float[n_vertices * 7];
    index_array = new unsigned int[n_indices * 2];

    const std::vector<unsigned int>& communities = *(source_graph->true_communities);

    WriteVertexClusters(communities);
    WriteColours(communities);
    WriteIndexArray();
}

void GraphicObject::WriteVertexClusters(const std::vector<unsigned int>& communities)
{
    // arrange communities in an outer circle around the origin
    // arrange vertices in a circle within their respective communities
    // assign coordinates for vertices
    // param1: scales community radii according to the community size
    // param2: lower_bound of community radius

    float pi = 3.141592f;
    const unsigned int n_communities = source_graph->n_communities;
    const float global_radius = 0.55f;
    const float global_angular_increment = 2.0f * pi / (float)n_communities;
    const float param1 = 1.0f;
    const float param2 = 0.3f;

    // membership: N_communities x N_members

    std::vector<std::vector<unsigned int>> membership(n_communities, std::vector<unsigned int>());

    for (unsigned int i = 0; i < n_vertices; ++i)
    {
        membership[communities[i]].push_back(i);
    }

    for (unsigned int i = 0; i < n_communities; ++i)
    {

        const float centre_x = global_radius * std::sin(float(i) * global_angular_increment);
		const float centre_y = global_radius * std::cos(float(i) * global_angular_increment);
        const float radius = ((param1 * membership[i].size()) / (n_vertices * (float)n_communities)) + param2;
        const float angular_increment = 2.0f * pi / membership[i].size();

        for (unsigned int j = 0; j < membership[i].size(); ++j)
        {
            const unsigned int vertex_ptr = membership[i][j];
            
            vertex_array[vertex_ptr * 7] = (radius * std::sin(j * angular_increment)) + centre_x;
            vertex_array[(vertex_ptr * 7) + 1] = (radius * std::cos(j * angular_increment)) + centre_y;
            vertex_array[(vertex_ptr * 7) + 2] = 0.0f;
        }
    }
}

void GraphicObject::WriteColours(const std::vector<unsigned int>& communities)
{
    const float& dimension = (float)source_graph->dimension;
    const unsigned int& n_communities = source_graph->n_communities;
    const float colour_increment = 3.0f / dimension;
    float colour_mixer[3] = {1.0f, 1.0f, 0.0f};
    unsigned int n_iters = 0;

    for (unsigned int i = 0; i < n_communities; ++i)
    {
        // Vary colours in YMCK space
        if (n_iters < dimension / 3)
        {
           colour_mixer[0] = 1.0f - n_iters * colour_increment;
           colour_mixer[1] = 1.0f;
           colour_mixer[2] = n_iters * colour_increment;
        }

        else if (n_iters < 2 * dimension / 3)
        {   
           colour_mixer[0] = n_iters * colour_increment - 1.0f;
           colour_mixer[1] = 2.0f - n_iters * colour_increment;
           colour_mixer[2] = 1.0f;
        }

        else
        {   
           colour_mixer[0] = 1.0f;
           colour_mixer[1] = n_iters * colour_increment - 2.0f;
           colour_mixer[2] = 3.0f - n_iters * colour_increment;
        }

        for (unsigned int j = 0; j < dimension; ++j)
        {
            if (i == communities[j])
            {
                vertex_array[7 * j + 3] = colour_mixer[0];
                vertex_array[7 * j + 4] = colour_mixer[1];
                vertex_array[7 * j + 5] = colour_mixer[2];
                vertex_array[7 * j + 6] = 1.0f;
                ++n_iters;
            }
        }
    }
}

void GraphicObject::WriteIndexArray()
{
    const std::vector<std::vector<unsigned int>> edges = source_graph->adjacency_matrix->GetEdges();

    for (unsigned int i = 0; i < edges.size(); ++i)
    {
        const unsigned int index = i * 2;
        index_array[index] = edges[i][0];
        index_array[index + 1] = edges[i][1];
    }
}

GraphicObject::GraphicObject(const Graph& graph)
:
source_graph(&graph)
{
    Init();
}

GraphicObject::~GraphicObject()
{
    delete vertex_array;
    delete index_array;
}