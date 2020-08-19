#include<iostream>
#include<cmath>
#include<vector>
#include"graphicobject.h"

void GraphicObject::WriteVertexClusters()
{
    // arrange communities in an outer circle around the origin
    // arrange vertices in a circle within their respective communities
    // assign coordinates for vertices
    // param1: scales community radii according to the community size
    // param2: lower_bound of community radius

    float pi = 3.141592f;
    const unsigned int n_vertices = source_graph->dimension;
    const unsigned int n_communities = source_graph->n_communities;
    std::vector<unsigned int>& communities = *(source_graph->true_communities);

    const float global_radius = 0.8f;
    const float global_angular_increment = 2.0f * pi / (float)n_communities;
    const float param1 = 2.5f;
    const float param2 = 0.2f;

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
            const unsigned int v = membership[i][j];
            
            vertices[v * 7] = (radius * std::sin(j * angular_increment)) + centre_x;
            vertices[(v * 7) + 1] = (radius * std::cos(j * angular_increment)) + centre_y;
            vertices[(v * 7) + 2] = 0.0f;
        }
    }
}

void GraphicObject::WriteIndices()
{
    const std::vector<std::vector<unsigned int>> edges = source_graph->adjacency_matrix->GetEdges();

    for (unsigned int i = 0; i < edges.size(); ++i)
    {
        const unsigned int index = i * 2;
        indices[index] = edges[i][0];
        indices[index + 1] = edges[i][1];

        std::cout << indices[index] << ", " << indices[index + 1] << std::endl;
    }
}

GraphicObject::GraphicObject(Graph &graph)
:
source_graph(&graph)
{
    vertices = new float[source_graph->dimension * 7];
    indices = new unsigned int[source_graph->n_edges * 2];
}

GraphicObject::~GraphicObject()
{
    delete vertices;
    delete indices;
}