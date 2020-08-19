#include<iostream>
#include<vector>
#include"graphicobject.h"

GraphicObject::GraphicObject(Graph &graph)
{
    vertices = new float[graph.dimension * 7];
    indices = new unsigned int[graph.n_edges];
}

GraphicObject::~GraphicObject()
{
    delete vertices;
    delete indices;
}