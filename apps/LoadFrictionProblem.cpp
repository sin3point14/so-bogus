
/*
 * Any copyright is dedicated to the Public Domain.
 * http://creativecommons.org/publicdomain/zero/1.0/
*/

#include "Interfaces/MecheInterface.hpp"

#include <Eigen/Core>

#include <iostream>
#include <cstdlib>

void usage(const char* name) {
	std::cout << "Usage: " << name << " dataFile [options] \n "
	          << "Options: \n"
	          << " -a int \t algorithm id in [0,4] \n"
	          << " -T int \t max number of threads (or 0 for OMP_MAX_THREADS) \n"
	          << " -m int \t max number of iterations \n"
	          << " -c int \t max number of Cadoux fixed-point iterations (0 means direct solving) \n"
	          << " -g int \t Projected Gradient variant in [0,4] \n"
	          << " -t real \t tolerance \n"
	          << " -r real \t problem regularization \n"
	          << " -p real \t (ADMM only) projection step size\n"
	          << " -f real \t (ADMM only) fixed-point step size\n"
	          << " -s bool \t if true, solve SOCQP instead of DCFP (static problem) \n"
	          << " -i bool \t if true, use infinity norm instead of l2 \n"
	          << " -o bool \t if true, use the old (<1.4) file format\n"
	          << " -z bool \t if true, GS will try to start at r=zero\n"
	          << " -k int  \t GS sleeping iterations \n"
	          << std::endl ;

}

int main( int argc, const char* argv[] )
{

	bogus::MecheFrictionProblem::Options options ;
	options.tryZeroAsWell = false ;

	const char* file = BOGUS_NULL_PTR(const char) ;
	double problemRegularization = 0. ;
	bool staticPb = false ;
	bool old = false ;

	for( int i = 1 ; i < argc ; ++i )
	{
		if( argv[i][0] == '-' ){
			switch(argv[i][1]) {
			case 'h':
				usage( argv[0]) ;
				return 0;
			case 'T':
				if( ++i == argc ) break ;
				options.maxThreads = std::atoi( argv[i] ) ;
				break ;
			case 'm':
				if( ++i == argc ) break ;
				options.maxIters = std::atoi( argv[i] ) ;
				break ;
			case 'c':
				if( ++i == argc ) break ;
				options.cadouxIters = std::atoi( argv[i] ) ;
				break ;
			case 't':
				if( ++i == argc ) break ;
				options.tolerance = std::strtod( argv[i], NULL ) ;
				break ;
			case 's':
				if( ++i == argc ) break ;
				staticPb = (bool) std::atoi( argv[i] ) ;
				break ;
			case 'i':
				if( ++i == argc ) break ;
				options.useInfinityNorm = (bool) std::atoi( argv[i] ) ;
				break ;
			case 'a':
				if( ++i == argc ) break ;
				options.algorithm =   (bogus::MecheFrictionProblem::Algorithm)
				        std::atoi( argv[i] ) ;
				break ;
			case 'g':
				if( ++i == argc ) break ;
				options.pgVariant =   (bogus::projected_gradient::Variant)
				        std::atoi( argv[i] ) ;
				break ;
			case 'r':
				if( ++i == argc ) break ;
				problemRegularization = std::strtod( argv[i], NULL ) ;
				break ;
			case 'o':
				if( ++i == argc ) break ;
				old = (bool) std::atoi( argv[i] ) ;
				break ;
			case 'p':
				if( ++i == argc ) break ;
				options.admmProjStepSize = std::strtod( argv[i], NULL ) ;
				break ;
			case 'f':
				if( ++i == argc ) break ;
				options.admmFpStepSize = std::strtod( argv[i], NULL ) ;
				break ;
			case 'k':
				if( ++i == argc ) break ;
				options.gsSkipIters = std::atoi( argv[i] ) ;
				break ;
			case 'z':
				if( ++i == argc ) break ;
				options.tryZeroAsWell = (bool) std::atoi( argv[i] ) ;
				break ;
			}
		} else {
			file = argv[i] ;
		}
	}

	if( !file )
	{
		std::cerr << " Please provide a problem data file " << std::endl ;
		usage(argv[0]) ;
		return 1 ;
	}


	bogus::MecheFrictionProblem mfp ;

	double * r = BOGUS_NULL_PTR(double) ;
	if( mfp.fromFile( file, r, old ) )
	{

		if( options.maxThreads > 1 )
			options.gsColoring = true ;

		if( !staticPb ) {
			options.gsRegularization = problemRegularization ;
			problemRegularization = 0 ;
		}


		const double res = mfp.solve( r, NULL, options, staticPb, problemRegularization ) ;
		std::cout << "Residual:\t "     << res << std::endl ;
		std::cout << "Solver timer:\t " << mfp.lastSolveTime() << " seconds" << std::endl ;

		delete[] r ;

		return 0 ;
	}

	std::cerr << " Could not load " << file << std::endl ;

	return 1 ;

}

