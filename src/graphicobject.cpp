#include<iostream>
#include<cmath>
#include<vector>
#include"graphicobject.h"

void GraphicObject::Init()
{
    n_vertices = source_graph->dimension;
    n_indices = source_graph->n_edges * 2;

    std::vector<unsigned int>& true_communities = *source_graph->true_communities;
    std::vector<unsigned int>& detected_communities = *source_graph->detected_communities;
    const unsigned int& n_communities = source_graph->n_communities;
    const unsigned int& n_comm_detected = source_graph->n_comm_detected;

    vertex_array = new float[n_vertices * 7];
    index_array = new unsigned int[n_indices * 2];

    //WriteVertexRings(detected_communities, n_comm_detected);
    WriteVertexRings(detected_communities, n_comm_detected);
    WriteColours(detected_communities, n_comm_detected);
    WriteIndexArray();
}

void GraphicObject::WriteVertexRings(const std::vector<unsigned int>& communities, unsigned int n_communities)
{
    // arrange communities in an outer circle around the origin
    // arrange vertices in a circle within their respective communities
    // assign coordinates for vertices
    // param1: scales community radii according to the community size
    // param2: lower_bound of community radius

    const float pi = 3.141592f;
    const float param1 = 0.8f;
    const float param2 = 0.95f;

    // membership: N_communities x N_members

    std::vector<std::vector<unsigned int>> membership(n_communities, std::vector<unsigned int>());

    for (unsigned int i = 0; i < n_vertices; ++i)
    {
        membership[communities[i]].push_back(i);
    }

    for (unsigned int i = 0; i < n_communities; ++i)
    {
        const float frac_coeff = (float)membership[i].size() / (float)source_graph->dimension;
        const float global_radius = std::pow((1.0f - frac_coeff), 2.0f) * param1;
        const float global_angular_increment = 2.0f * pi / (float)n_communities;
        const float centre_x = global_radius * std::sin(float(i) * global_angular_increment);
		const float centre_y = global_radius * std::cos(float(i) * global_angular_increment);
        const float radius = std::pow(frac_coeff * param2, 0.5f) * param2;
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

void GraphicObject::WriteVertexDynamics()
{
    // arrange vertices based on Hooke's law
    // Fij = -k (M(ri - rj) - r0) * unit(ri - rj)
    // vij = - k * (M(ri - rj) - r0))^2 * unit(ri - rj) / 2
    // sij = - k * (M(ri - rj) - r0))^2 * unit(ri - rj) * dt / 2

    // initialise lattice structure

    const float separation2D = 2.0f / (std::sqrt((float)n_vertices) + 1.0f); 
    const unsigned int n_vertices_line = (int)(std::sqrt((float)n_vertices) + 1.0f);

    for (unsigned int i = 0; i < n_vertices_line; ++i)
    {
        float x_coord = -1.0f + (((float)i + 1.0f) * separation2D);
        float z_coord = 0.0f;

        for (unsigned int j = 0; j < n_vertices_line; ++j)
        {
            unsigned int vertex = (i * n_vertices_line) + j;

            if (vertex >= n_vertices)
            {
                break;
            }

            unsigned int vertex_ptr = vertex * 7;
            float y_coord = -1.0f + (((float)j + 0.5f) * separation2D) ;

            vertex_array[vertex_ptr] = x_coord;
            vertex_array[vertex_ptr + 1] = y_coord;
            vertex_array[vertex_ptr + 2] = z_coord;
        }
    }

    std::vector<std::vector<unsigned int>> edges = source_graph->adjacency_matrix->GetEdges();
    
    // simulate dynamics and reset velocity to 0 every timestep
    // add forcefield to the screen edges
    const unsigned int n_timesteps = 10000;
    const float dt = 0.01;
    const float r0 = 0.2;
    const float k = 0.0001f;
    std::vector<float> forces_x(n_vertices);
    std::vector<float> forces_y(n_vertices);
    unsigned int t = 0;

    while (t < n_timesteps)
    {
        // get forces felt by every vertex 
        for (unsigned int i = 0; i < n_vertices; ++i)
        {
            unsigned int vertex_i_ptr = i * 7;
            float& r_i_x = vertex_array[vertex_i_ptr];
            float& r_i_y = vertex_array[vertex_i_ptr + 1];

            forces_x[i] = 0.0f;
            forces_y[i] = 0.0f;

            std::vector<unsigned int> adjacent = source_graph->adjacency_matrix->GetAdjacentVertices(i);

            for (unsigned int j = 0; j < adjacent.size(); ++j)
            {
                unsigned int vertex_j_ptr = j * 7;
                float& r_j_x = vertex_array[vertex_j_ptr];
                float& r_j_y = vertex_array[vertex_j_ptr + 1];

                forces_x[i] -= ((r_i_x - r_j_x) - r0) * k;
                forces_y[i] -= ((r_i_y - r_j_y) - r0) * k;
            }
        } 

        // affect forces on every vertex
        for (unsigned int i = 0; i < n_vertices; ++i)
        {
            unsigned int vertex_ptr = i * 7;
            float& r_x = vertex_array[vertex_ptr];
            float& r_y = vertex_array[vertex_ptr + 1];

            float s_x  = - pow(forces_x[i], 2.0f) * dt / (2.0f * k);
            float s_y  = - pow(forces_y[i], 2.0f) * dt / (2.0f * k);

            if (forces_x[i] < 0)
            {
                s_x *= -1;
            }

            if (forces_y[i] < 0)
            {
                s_y *= -1;
            }

            r_x += s_x ;
            r_y += s_y;

            float orbit = VectorMagnitude(r_x, r_y);

            if (orbit > 0.9f)
            {
                r_x /= (orbit + 0.1f);
                r_y /= (orbit + 0.1f);
            }
        }

        ++t;
    }
}

void GraphicObject::WriteColours(const std::vector<unsigned int>& communities, unsigned int n_communities)
{
    const float& dimension = (float)source_graph->dimension;
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

float GraphicObject::VectorMagnitude(float i, float j)
{
    return std::sqrt(std::pow(i, 2.0f) + std::pow(j, 2.0f));
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