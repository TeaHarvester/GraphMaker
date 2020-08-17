#include<iostream>
#include"graph.h"

int main()
{
    Graph testgraph; 
    testgraph.GenerateLFRGraph(20, 2, 10, 2, 3, 2, 0.8);
    testgraph.adjacency_matrix->Print();

    return 0;
}