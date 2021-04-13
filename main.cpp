#include <cstdlib>
#include <cstdio>
#include <stdexcept>
#include "option_handler.h"
#include "dgp.h"
#include "sph.h"

int main(int argc, char *argv[])
{
	option_handler options(argc, argv);
	char fn_dgp[256];		// path of dgp file
	double ftol = 1E-2;		// feasibility tolerance
	int numTests = 10; // number of initial points
	bool verbose = true;	// verbose mode

	if (!options.get_value("-dgp", fn_dgp))
	{
		throw std::runtime_error("the option -dgp <dgp file path> must be provided");
	}
	options.get_value("-ftol", ftol);
	options.get_value("-verbose", verbose);
	options.get_value("-ntests", numTests);

	if (verbose)
	{
		printf("ftol .... %g\n", ftol);
		printf("ntests .. %d\n", numTests);
	}

	dgp_t G(fn_dgp);

	double *x = new double[3 * G.m_nnodes];
	double f;

	int numSolvedTests = 0;

	for (int i = 0; i < numTests; i++)
	{
		printf("test: %d out of %d\n", i, numTests);
		if (sph(G, ftol, f, x, verbose))
		{
			numSolvedTests++;
		}
	}

	printf("nSolvedTests ... %d out of %d\n", numSolvedTests, numTests);

	delete[] x;
	return EXIT_SUCCESS;
}
