#ifndef SPARSEMATRIX
#define SPARSEMATRIX
#include<vector>

struct SparseMatrix
{
    unsigned int dimension;
    std::vector<float> val;
    std::vector<unsigned int> col;
    std::vector<unsigned int> rowptr;

    std::vector<unsigned int> GetAdjacentVertices(const unsigned int vertex);
    bool IsAdjacent(const unsigned int vertex1, const unsigned int vertex2);
    void AddConnection(unsigned int row, unsigned int column, float value);
    void Print();
    SparseMatrix(unsigned int dim);
};

#endif