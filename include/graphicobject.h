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
    void WriteVertexClusters(const std::vector<unsigned int>& communities);
    void WriteColours(const std::vector<unsigned int>& communities);
    void WriteIndexArray();

    public:
    GraphicObject(const Graph& graph);
    ~GraphicObject();
};

#endif