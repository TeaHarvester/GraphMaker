#ifndef GRAPHIC_OBJECT
#define GRAPHIC_OBJECT

#include<iostream>
#include"graph.h"

// object to be rendered by the render function
// vertices: N x 7 array of position (xyz) and colours (RGBa)
// indices: V x 2 array for drawing lines between vertices

class GraphicObject
{
    public:
    unsigned int n_vertices;
    unsigned int n_indices; 
    float* vertex_array;
    unsigned int* index_array;
    const Graph* source_graph;

    void Init();

    private:
    void WriteVertexRings(const std::vector<unsigned int>& communities, unsigned int n_com);
    void WriteVertexDynamics();
    void WriteColours(const std::vector<unsigned int>& communities, unsigned int n_com);
    void WriteIndexArray();
    float VectorMagnitude(float i, float j);

    public:
    GraphicObject(const Graph& graph);
    ~GraphicObject();
};

#endif