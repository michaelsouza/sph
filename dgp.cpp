#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <algorithm>
#include "dgp.h"
#include "vec.h"

class edge_t
{
    public:
    int m_i;
    int m_j;
    double m_l;
    double m_u;

    edge_t(int i, int j, double l, double u)
    {
        m_i = i;
        m_j = j;
        m_l = l;
        m_u = u;
    }
};

void dgp_t::convertCOO2CSR()
{
    sortEdges();

    /// converting to CSR
    std::vector<int> coo_i(m_vec_i);

    m_vec_i.resize(m_nnodes+1);
    m_i = &(m_vec_i[0]);

    /// reset i
    for (int k = 0; k <= m_nnodes; k++)
    {
        m_i[k] = 0;
    }

    /// count the number of edges per row
    for (int k = 0; k < m_nedges; k++)
    {
        m_i[coo_i[k] + 1]++;
    }

    /// cumulative sum
    for (int k = 2; k <= m_nnodes; k++)
    {
        m_i[k] += m_i[k - 1];
    }
}

// sort using bubble sort
void dgp_t::sortEdges()
{
    std::vector<edge_t> edges;
    for (int k = 0; k < m_nedges; ++k )
    {
        edges.push_back(edge_t(m_i[k], m_j[k], m_l[k], m_u[k]));
    }

    std::sort(edges.begin(), edges.end(),
              [](edge_t &a, edge_t &b) {
                  if (a.m_i > b.m_i)
                      return false;
                  if (a.m_i == b.m_i && a.m_j > b.m_j)
                      return false;
                  return true;
              });

    for (int k = 0; k < m_nedges; ++k )
    {
        m_i[k] = edges[k].m_i;
        m_j[k] = edges[k].m_j;
        m_l[k] = edges[k].m_l;
        m_u[k] = edges[k].m_u;
    }
}

void dgp_t::neighs(int k, int &n, int **neighs)
{
    n = m_i[k + 1] - m_i[k];
    (*neighs) = &(m_j[m_i[k]]);
}

void dgp_t::print()
{
    int i, k;
    printf("-------------------------\n");
    printf("Number of nodes: %d\n", m_nnodes);
    printf("Number of edges: %d\n", m_nedges);
    printf("-------------------------\n");
    printf("nid   i  j  dij  lij  uij\n");
    for (i = 0; i < m_nnodes; i++)
    {
        for (k = m_i[i]; k < m_i[i + 1]; k++)
        {
            printf("[%2u] %2u %2u %3.2f %3.2f\n", k, i, m_j[k], m_l[k], m_u[k]);
        }
    }

    printf("-------------------------\n");
    printf("    x       y       z\n");
    if (m_x)
    {
        for (i = 0; i < m_nnodes; i++)
        {
            printf("% .5f, % .5f, % .5f\n", m_x[3 * i], m_x[3 * i + 1],
                   m_x[3 * i + 2]);
        }
    }
    printf("-------------------------\n");
}

bool dgp_t::vals(int i, int j, double *dij, double *lij, double *uij)
{
    int k;
    /// swap
    if (i > j)
    {
        k = i;
        i = j;
        j = k;
    }
    /// find values
    for (k = m_i[i]; k < m_i[i + 1]; k++)
    {
        if (m_j[k] == j)
        {            
            if (lij) *lij = m_l[k];
            if (uij) *uij = m_u[k];
            return true;
        }
    }
    return false;
}

bool dgp_t::feasible(double *x, int i, double xtol)
{
    int k, j;
    double err_max = 0, *xi, *xj, dij, lij, uij, err_ij;

    xi = x + 3 * i;
    for (k = m_i[i]; k < m_i[i + 1]; k++)
    {
        if ((j = m_j[k]) > i)
            break;
        xj = x + 3 * j;
        dij = vec_dist(xi, xj);
        err_ij = 0.0;
        if (dij < (lij = m_l[k]))
        {
            err_ij = (lij - dij) / lij;
        }
        else if (dij > (uij = m_u[k]))
        {
            err_ij = (dij - uij) / uij;
        }

        if (err_max < err_ij)
        {
            err_max = err_ij;
        }
    }
    return err_max < xtol;
}

bool dgp_t::feasible(double *x, double xtol)
{
    for (int i = 0; i < m_nnodes; i++)
    {
        if (!feasible(x, i, xtol))
        {
            return false;
        }
    }
    return true;
}
