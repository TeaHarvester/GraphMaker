#ifndef SPARSEMATRIX
#define SPARSEMATRIX
#include<vector>

struct SparseMatrix
{
    unsigned int dimension;
    std::vector<float> val;
    std::vector<unsigned int> col;
    std::vector<unsigned int> rowptr;

    void AddConnection(const unsigned int row, const unsigned int column, float value);
    void EraseConnection(const unsigned int row, const unsigned int column);
    const std::vector<unsigned int> GetAdjacentVertices(const unsigned int vertex);
    const std::vector<std::vector<unsigned int>> GetEdges();
    bool IsAdjacent(const unsigned int vertex1, const unsigned int vertex2);
    float GetEdgeWeight(unsigned int vertex_1, unsigned int vertex_2);
    unsigned int GetDegree(unsigned int vertex);
    void Print();
    SparseMatrix(unsigned int dim);
    SparseMatrix(const SparseMatrix& S);
};

#endif