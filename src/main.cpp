#include<iostream>
#include"graph.h"

int main()
{
    Graph morphy; 
    morphy.GenerateLFRGraph(20, 2, 10, 2);
    morphy.adjacency_matrix->Print();
    // for (int i = 0; i < morphy.dimension; ++i)
    // {
    //     std::cout << morphy.adjacency_matrix->col[i] << std::endl;
    // }

    return 0;
}