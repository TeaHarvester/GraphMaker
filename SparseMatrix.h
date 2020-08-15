#IFNDEF SPARSEMATRIX
#DEfiNE SPARSEMATRIX
#include<vector>

struct SparseMatrix
{
    std::vector<float> val;
    std::vector<unsigned int> col;
    std::vector<unsigned int> rowptr;
}

#ENDIF