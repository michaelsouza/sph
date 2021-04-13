#include <algorithm>
#include <cstdio>
#include <random>
#include <gsl/gsl_multimin.h>

#include "dgp.h"
#include "vec.h"
#include "fobj_sph.h"

typedef struct
{
   double tau;
   double lambda;
   double rho;
   dgp_t *G;
} gsl_params_t;

void gsl_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *g)
{
   gsl_params_t *p = (gsl_params_t *)params;
   *f = fobj_sph(p->lambda, p->tau, *(p->G), x->data, g->data);
}

void gsl_df(const gsl_vector *x, void *params, gsl_vector *g)
{
   gsl_params_t *p = (gsl_params_t *)params;
   fobj_sph(p->lambda, p->tau, *(p->G), x->data, g->data);
}

double gsl_f(const gsl_vector *x, void *params)
{
   gsl_params_t *p = (gsl_params_t *)params;
   return fobj_sph(p->lambda, p->tau, *(p->G), x->data, NULL);
}

double init_tau(dgp_t &G)
{
   std::vector<double> tau_all(G.m_nedges);

   for (int i = 0; i < G.m_nedges; i++)
   {
      tau_all[i] = (G.m_l[i] + G.m_u[i]) / 2.0;
   }

   std::sort(tau_all.begin(), tau_all.end());
    
   return tau_all[(int)(0.75 * G.m_nedges)];;
}

double gsl_vector_norm(gsl_vector *v)
{
   double v_nrm = 0.0;
   for (size_t i = 0; i < v->size; i++)
   {
      v_nrm += v->data[i] * v->data[i];
   }
   return sqrt(v_nrm);
}

void init_sol(dgp_t &G, double *x)
{
   double dij, lij, uij, *xi, *xj;
   double y[3], y_nrm;
   int *neighs, nneighs, j;
   std::vector<int> s(G.m_nnodes);
   std::vector<bool> b(G.m_nnodes, false);

   // s = [0, 1, 2, ..., nnodes-1]
   for (int k = 0; k < G.m_nnodes; k++)
   {
      s[k] = k;
   }   
   std::random_shuffle(s.begin(), s.end());

   // x = zeros(3 * nnodes)
   for (int k = 0; k < (3 * G.m_nnodes); ++k)
   {
      x[k] = 0.0;
   }

   std::default_random_engine rndEngine;
   std::normal_distribution<double> gaussian(0.0, 1.0);

   for (int k = 0; k < G.m_nnodes; k++)
   {
      int i = s[k];
      xi = &(x[3 * i]);
      G.neighs(i, nneighs, &neighs);
      for (int jj = 0; jj < nneighs; jj++)
      {
         j = neighs[jj];
         if (b[j])
         {
            continue;
         }
         G.vals(i, j, NULL, &lij, &uij);
         dij = (lij + uij) / 2.0;
         xj = &(x[3 * j]);
         // y is randomly distributed inside the sphere((0,0,0),1)
         vec_set(y, gaussian(rndEngine), gaussian(rndEngine), gaussian(rndEngine));
         y_nrm = vec_norm(y);
         // xj = dij * (y / y_nrm) + xi  (randomly distributed over sphere(xi, dij))
         vec_axpby(dij / y_nrm, y, 1.0, xi, xj);
      }
   }
}

void sph_print_iter(int i, gsl_params_t &p, double f, double fs, gsl_vector *g, int k)
{
   if ((i % 10) == 0)
   {
      printf(" iter |    tau   |    f     |    fs    |  nrm(g)  | #f_calls \n");
   }
   printf("%5d | %3.2e | %3.2e | %3.2e | %3.2e | %5d\n", i, p.tau, f, fs, gsl_vector_norm(g), k);
}

bool sph(dgp_t &G, double ftol, double &f, double *x, bool verbose)
{
   gsl_params_t gsl_params;
   gsl_params.lambda = 0.5;
   gsl_params.rho = 0.99;
   gsl_params.tau = init_tau(G);
   gsl_params.G = &G;

   gsl_vector *y = gsl_vector_alloc(3 * G.m_nnodes);
   gsl_vector *g = gsl_vector_alloc(3 * G.m_nnodes);

   init_sol(G, y->data);

   const gsl_multimin_fdfminimizer_type *fdfmin_method = gsl_multimin_fdfminimizer_vector_bfgs2;
   gsl_multimin_fdfminimizer *fdfmin;

   gsl_multimin_function_fdf gsl_fobj;
   gsl_fobj.n = 3 * G.m_nnodes;
   gsl_fobj.f = gsl_f;
   gsl_fobj.df = gsl_df;
   gsl_fobj.fdf = gsl_fdf;
   gsl_fobj.params = (void *)(&gsl_params);

   double fs;
   bool solved = false;
   int maxit = 1000, k = 0, status;
   for (int i = 0; i < maxit; i++)
   {
      fs = fobj_sph(gsl_params.lambda, gsl_params.tau, G, y->data, g->data);
      f = fobj_sph(gsl_params.lambda, 1E-16, G, y->data, NULL);

      if (verbose)
      {
         sph_print_iter(i, gsl_params, f, fs, g, k);
      }

      if (f < ftol)
      {
         solved = true;
         break;
      }
      
      gsl_params.tau *= gsl_params.rho;
      fdfmin = gsl_multimin_fdfminimizer_alloc(fdfmin_method, 3 * G.m_nnodes);
      gsl_multimin_fdfminimizer_set(fdfmin, &gsl_fobj, y, 0.01, 1E-3);

      // call local optimization
      for (k = 0, status = GSL_CONTINUE; k < 1000 && status == GSL_CONTINUE; k++)
      {
         status = gsl_multimin_fdfminimizer_iterate(fdfmin);
         if (status)
         {
            break;
         }
         status = gsl_multimin_test_gradient(fdfmin->gradient, 1E-8);
      }

      // copy from fdfmin to y
      for (int j = 0; j < 3 * G.m_nnodes; j++)
      {
         y->data[j] = fdfmin->x->data[j];
      }

      gsl_multimin_fdfminimizer_free(fdfmin);
   }

   // copy solution from y to x
   for (int i = 0; i < 3 * G.m_nnodes; i++)
   {
      x[i] = y->data[i];
   }

   gsl_vector_free(y);

   return solved;
}
