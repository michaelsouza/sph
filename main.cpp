#include "dgp.h"
#include "option_handler.h"
#include "sph.h"
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>


int main( int argc, char* argv[] ) {
   option_handler options( argc, argv );
   char fn_dgp[ 256 ];   // path of dgp file
   double ftol = 1E-2;   // feasibility tolerance
   int numTests = 10;    // number of initial points
   bool verbose = false; // verbose mode

   if ( !options.get_value( "-dgp", fn_dgp ) ) {
      throw std::runtime_error(
          "the option -dgp <dgp file path> must be provided" );
   }
   options.get_value( "-ftol", ftol );
   options.get_value( "-verbose", verbose );
   options.get_value( "-ntests", numTests );

   if ( verbose ) {
      printf( "ftol .... %g\n", ftol );
      printf( "ntests .. %d\n", numTests );
   }

   dgp_t G( fn_dgp );
   G.show();

   double* x = new double[ 3 * G.m_nnodes ];
   double f;

   int numSolvedTests = 0;

   printf( "tid\t fx\t maxRelRes\t tElapsed(secs)\n" );
   printf( "-----------------------------------------\n" );
   for ( int i = 0; i < numTests; i++ ) {
	  double tic = omp_get_wtime();
      bool solved = sph( G, ftol, f, x, verbose );
	  double toc = omp_get_wtime() - tic;
      printf( "%4d %10.3E %10.3E %12.2f\n", i + 1, f, G.maxRelRes( x ), toc );
      if ( solved ) {
         numSolvedTests++;
      }
   }
   printf( "-----------------------------------------\n" );

   printf( "nSolvedTests ... %d out of %d\n", numSolvedTests, numTests );

   delete[] x;
   return EXIT_SUCCESS;
}
