#pragma once

#include <string>
#include <stdexcept>
#include <vector>
#include <fstream>
#include <sstream>

class dgp_t
{
 public:
	int m_nnodes;
	int m_nedges;
	int *m_i; // edge index
	int *m_j; // edge index
	double *m_l;
	double *m_u;
	double *m_x; // solution

	dgp_t(int nnodes, int nnedges, int *i, int *j, double *l, double *u, double *x)
	{
		m_nnodes = nnodes;
		m_nedges = nnedges;
		m_i = i;
		m_j = j;
		m_l = l;
		m_u = u;

		convertCOO2CSR();
	}

	dgp_t(std::string fn_dgp)
	{		
		printf("Reading file %s\n", fn_dgp.c_str());				
		std::ifstream fid(fn_dgp);

		int i, j;
		double l, u;
		m_nnodes = 0;

		for (std::string row; getline( fid, row); )
		{
			std::stringstream ss(row);
			ss >> i >> j >> l >> u;

			// append (i, j)
			m_vec_i.push_back(i - 1); // convert from one to zero-based
			m_vec_j.push_back(j - 1); 
			m_vec_l.push_back(l);
			m_vec_u.push_back(u);

			// G is undirected: append (j, i)
			m_vec_i.push_back(j - 1);
			m_vec_j.push_back(i - 1); 
			m_vec_l.push_back(l);
			m_vec_u.push_back(u);

			if (m_nnodes < i)
				m_nnodes = i;
			if (m_nnodes < j)
				m_nnodes = j;
		}		
		m_nedges = (int) m_vec_i.size();

		m_i = &(m_vec_i[0]);
		m_j = &(m_vec_j[0]);
		m_l = &(m_vec_l[0]);
		m_u = &(m_vec_u[0]);

		convertCOO2CSR();
	}

	void print();

	bool vals(int i, int j, double *dij, double *lij, double *uij);

	void neighs(int k, int &n, int **neighs);

	void convertCOO2CSR();

	void sortEdges();

	bool feasible(double *x, int i, double xtol);

	bool feasible(double *x, double xtol);

private:
	std::vector<double> m_vec_u;
	std::vector<double> m_vec_l;
	std::vector<int> m_vec_i;
	std::vector<int> m_vec_j;
};
