#ifndef SPARSEMATRIX
#define SPARSEMATRIX
#include<vector>

struct SparseMatrix
{
    unsigned int dimension;
    std::vector<float> val;
    std::vector<unsigned int> col;
    std::vector<unsigned int> rowptr;

    private: 
    const std::vector<unsigned int> GetAdjacentVertices(const unsigned int vertex);

    public:
    const std::vector<std::vector<unsigned int>> GetEdges();
    bool IsAdjacent(const unsigned int vertex1, const unsigned int vertex2);
    void AddConnection(const unsigned int row, const unsigned int column, float value);
    float GetEdgeWeight(unsigned int vertex_1, unsigned int vertex_2);
    void Print();
    SparseMatrix(unsigned int dim);
};

#endif