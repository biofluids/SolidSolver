#include <vector>
#include <set>
#include <iostream>

extern "C"
{
    void getCRSInfo1(int nsd, int nn, int nel, int nen, int* arr, int* num);
    void getCRSInfo2(int nsd, int nn, int nel, int nen, int* arr, int* rowPtrs, int* colInds);
}
void getNeighbors(int nsd, int nn, int& no_nonzeros, const std::vector<std::vector<int> >& connect, \
    std::vector<std::set<int> >& neighborNodes);
void getCRS(int nsd, int nn, std::vector<std::set<int> >& neighborNodes,\
    std::vector<int>& row_ptr, std::vector<int>& col_ind);

void getCRSInfo1(int nsd, int nn, int nel, int nen, int* arr, int* num)
{
    // step 1, query the number of nonzeros
    std::vector<std::vector<int> > connect;
    for (int i = 0; i < nel; ++i)
    {
        std::vector<int> eleConnect;
        for (int j = 0; j < nen; ++j)
            eleConnect.push_back(arr[i*nen+j]);
        connect.push_back(eleConnect);
    }
    int no_nonzeros = 0;
    std::vector<std::set<int> > neighborNodes; // Storing node numbers beginning with 1
    getNeighbors(nsd, nn, no_nonzeros, connect, neighborNodes);
    *num = no_nonzeros;
}

void getCRSInfo2(int nsd, int nn, int nel, int nen, int* arr, int* rowPtrs, int* colInds)
{
    std::vector<std::vector<int> > connect;
    for (int i = 0; i < nel; ++i)
    {
        std::vector<int> eleConnect;
        for (int j = 0; j < nen; ++j)
            eleConnect.push_back(arr[i*nen+j]);
        connect.push_back(eleConnect);
    }
    int no_nonzeros = 0;
    std::vector<std::set<int> > neighborNodes; // Storing node numbers beginning with 1
    getNeighbors(nsd, nn, no_nonzeros, connect, neighborNodes);

    std::vector<int> row_ptr, col_ind;
    getCRS(nsd, nn, neighborNodes, row_ptr, col_ind);

    for (int i = 0; i < row_ptr.size(); ++i)
        rowPtrs[i] = row_ptr[i];
    for (int i = 0; i < col_ind.size(); ++i)
        colInds[i] = col_ind[i];
}

void getNeighbors(int nsd, int nn, int& no_nonzeros, const std::vector<std::vector<int> >& connect, \
    std::vector<std::set<int> >& neighborNodes)
{
    neighborNodes.clear();
    int nel = connect.size();
    int nen = connect[0].size();
    no_nonzeros = 0;
    for (int n = 0; n < nn; ++n)
    {
        std::set<int> my_neighborNodes;
        for (int ele = 0; ele < nel; ++ele)
        {
            for (int i = 0; i < nen; ++i)
            {
                if (connect[ele][i] == n + 1)
                {
                    my_neighborNodes.insert(connect[ele].begin(), connect[ele].end());
                    break;
                }
            }
        }
        neighborNodes.push_back(my_neighborNodes);
        no_nonzeros += my_neighborNodes.size()*nsd;
    }
    no_nonzeros = no_nonzeros*nsd; // For each node, every dof is identical.
}

void getCRS(int nsd, int nn, std::vector<std::set<int> >& neighborNodes,\
    std::vector<int>& row_ptr, std::vector<int>& col_ind)
{
    int count = 0; // the number of nonzeros
    // First nn*nsd rows
    for (int n = 0; n < nn; ++n)
    {
        // Deal with nsd rows corresponding to one node
        for (int i = 0; i < nsd; ++i)
        {
            row_ptr.push_back(count);
            // Stress terms
            for (auto j : neighborNodes[n])
            {
                for (int k = 0; k < nsd; ++k)
                {
                    col_ind.push_back(nsd*(j-1) + k);
                    count++;
                }
            }
        }
    }
    row_ptr.push_back(count);
    // Increment col_ind and row_ptr for FORTRAN code to use
    for (auto& i : col_ind) i++;
    for (auto& i : row_ptr) i++;
}
