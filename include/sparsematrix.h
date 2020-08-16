#ifndef SPARSEMATRIX
#define SPARSEMATRIX
#include<vector>

struct SparseMatrix
{
    unsigned int dimension;
    std::vector<float> val;
    std::vector<unsigned int> col;
    std::vector<unsigned int> rowptr;

    SparseMatrix(unsigned int dim);
    void AddConnection(unsigned int row, unsigned int column, float value);
    void Print();
};

#endif