#include<iostream>
#include"sparsematrix.h"

void SparseMatrix::AddConnection(unsigned int row, unsigned int column, float value)
{
    unsigned int element_counter = 0;

    for (unsigned int i = 0; i < row; ++i)
    {
        element_counter += rowptr[i + 1] - rowptr[i];
    }

    std::vector<unsigned int>::iterator col_iter;
    col_iter = col.begin() + element_counter;
    col.insert(col_iter, column);

    std::vector<float>::iterator val_iter;
    val_iter = val.begin() + element_counter; 
    val.insert(val_iter, value);

    for (unsigned int i = row; i < dimension; ++i)
    {
        ++rowptr[i + 1];
    }
}

void SparseMatrix::Print()
{
    std::cout << "i j v q" << std::endl;
    unsigned int col_iter = 0;
    unsigned int row_iter = 0;

    for (unsigned int i = 0; i < dimension; ++i)
    {
        unsigned int row_counter = rowptr[i + 1] - rowptr[i];

        for (unsigned int r = 0; r < row_counter; ++r)
        {
            std::cout << i << " " << col[col_iter] << " " << val[col_iter] << " ";
            
            if (row_iter < dimension + 1)
            {
                std::cout << rowptr[row_iter] << std::endl;
            }
            
            else 
            {
                std::cout << std::endl;
            }

            ++col_iter;
            ++row_iter;
        }
    }      

    std::cout << std::endl;
}

SparseMatrix::SparseMatrix(unsigned int dim)
: 
dimension(dim),
val(),
col(),
rowptr(dim + 1, 0)
{}