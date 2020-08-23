#include<iostream>
#include"sparsematrix.h"

void SparseMatrix::AddConnection(unsigned int row, unsigned int column, float value)
{
    if (value == 0)
    {
        return;
    }

    if (IsAdjacent(row, column))
    {
        unsigned int adjacent_ptr = rowptr[row];
        unsigned int n_adjacent = rowptr[row + 1] - adjacent_ptr;

        for (unsigned int i = adjacent_ptr; i < adjacent_ptr + n_adjacent; ++i)
        {
            if (col[i] == column)
            {
                ++val[i];
                return;
            }
        }
    }

    col.insert(col.begin() + rowptr[row], column);
    val.insert(val.begin() + rowptr[row], value);

    for (unsigned int i = row; i < dimension; ++i)
    {
        ++rowptr[i + 1];
    }
}

void SparseMatrix::EraseConnection(const unsigned int row, const unsigned int column)
{
    if (!IsAdjacent(row, column))
    {
        std::cout << "no edge to erase" << std::endl;
    }

    for (unsigned int i = rowptr[row]; i < rowptr[row + 1]; ++i)
    {
        if (col[i] == column)
        {
            col.erase(col.begin() + i);
            val.erase(val.begin() + i);
        }
    }

    for (unsigned int i = row; i < dimension; ++i)
    {
        --rowptr[i + 1];
    }
}

const std::vector<unsigned int> SparseMatrix::GetAdjacentVertices(const unsigned int vertex)
{
    unsigned int adjacent_ptr = rowptr[vertex];
    unsigned int n_adjacent = rowptr[vertex + 1] - adjacent_ptr;
    std::vector<unsigned int> output;

    for (unsigned int i = adjacent_ptr; i < adjacent_ptr + n_adjacent; ++i)
    {
        output.push_back(col[i]);
    }

    return output;
}

const std::vector<std::vector<unsigned int>> SparseMatrix::GetEdges()
{
    std::vector<std::vector<unsigned int>> output(col.size(), std::vector<unsigned int>(2));

    unsigned int n_iters = 0;

    for (unsigned int i = 0; i < dimension; ++i)
    {
        std::vector<unsigned int> adjacent_vertices = GetAdjacentVertices(i);

        for (unsigned int j = 0; j < adjacent_vertices.size(); ++j)
        {
            output[n_iters][0] = i;
            output[n_iters][1] = adjacent_vertices[j];
            ++n_iters;
        }
    }

    return output;
}

bool SparseMatrix::IsAdjacent(const unsigned int vertex1, const unsigned int vertex2)
{
    std::vector<unsigned int> adjacent_vertices = GetAdjacentVertices(vertex1);

    for (unsigned int i = 0; i < adjacent_vertices.size(); ++i)
    {
        if (adjacent_vertices[i] == vertex2)
        {
            return true;
        }
    }

    return false;
}

float SparseMatrix::GetEdgeWeight(unsigned int vertex_1, unsigned int vertex_2)
{
    unsigned int adjacent_ptr = rowptr[vertex_1];
    unsigned int n_adjacent = rowptr[vertex_1 + 1] - adjacent_ptr;
    std::vector<unsigned int> output;

    for (unsigned int i = adjacent_ptr; i < adjacent_ptr + n_adjacent; ++i)
    {
        if (col[i] == vertex_2)
        {
            return val[i];
        }
    }

    return 0.0f;
} 

unsigned int SparseMatrix::GetDegree(unsigned int vertex)
{
    std::vector<unsigned int> adjacent_vertices = GetAdjacentVertices(vertex);

    return adjacent_vertices.size();
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

SparseMatrix::SparseMatrix(const SparseMatrix& S) 
:
dimension(S.dimension),
val(S.val),
col(S.col),
rowptr(S.rowptr)
{}