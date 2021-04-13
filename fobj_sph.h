#pragma once

#include "dgp.h"
#include "vec.h"

inline double theta(double tau, double x[3], double g[3])
{
   double f = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + tau * tau);
   vec_scale(1 / f, x, g);
   return f;
}

inline double phi(double lambda, double tau, double x, double &g)
{
   double z = sqrt(lambda * lambda * x * x + tau * tau);
   g = lambda + (x * lambda * lambda) / z;
   return lambda * x + z;
}

double fobj_sph(double lambda, double tau, dgp_t &G, double *x, double *g)
{
   double f = 0.0;
   if (g != NULL)
   {
      for (int k = 0; k < (3 * G.m_nnodes); k++)
      {
         g[k] = 0.0;
      }
   }
   double lk, uk, *xi, *xj, vk[3], gk[3], *gi, *gj;
   double f_tta_k, f_ph_lk, f_ph_uk, g_ph_lk, g_ph_uk;
   for (int i = 0; i < G.m_nnodes; i++)
   {
      for (int k = G.m_i[i]; k < G.m_i[i + 1]; k++)
      {
         int j = G.m_j[k];
         // G is symmetric, but tril(G) is enough.
         if (j > i)
            break;
         lk = G.m_l[k];
         uk = G.m_u[k];
         xi = &(x[3 * i]);
         xj = &(x[3 * j]);
         gi = g != NULL ? &(g[3 * i]) : NULL;
         gj = g != NULL ? &(g[3 * j]) : NULL;
         vec_sub(xi, xj, vk);

         f_tta_k = theta(tau, vk, gk);
         f_ph_lk = phi(lambda, tau, lk - f_tta_k, g_ph_lk);
         f_ph_uk = phi(lambda, tau, f_tta_k - uk, g_ph_uk);

         f += f_ph_lk + f_ph_uk;
         
         // printf("%d %d %g %g %g %g %g\n", i+1, j+1, lk, uk, f, f_ph_lk, f_ph_uk);
         
         if (g != NULL)
         {
            vec_scale(g_ph_uk - g_ph_lk, gk);
            vec_add(gk, gi);
            vec_sub(gk, gj);
         }
      }
   }
   return f;
}
