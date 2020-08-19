#ifndef GRAPHIC_OBJECT
#define GRAPHIC_OBJECT

#include<iostream>
#include"graph.h"

// object to be rendered by the render function
// vertices: N x 7 array of position (xyz) and colours (RGBa)
// indices: V x 2 array for drawing lines between vertices

struct GraphicObject
{
    float* vertices;
    unsigned int* indices;
    Graph* source_graph;

    void WriteVertexClusters();
    void WriteIndices();
    GraphicObject(Graph &graph);
    ~GraphicObject();
};

#endif