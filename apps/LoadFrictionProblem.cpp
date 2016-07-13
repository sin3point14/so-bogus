
/*
 * Any copyright is dedicated to the Public Domain.
 * http://creativecommons.org/publicdomain/zero/1.0/
*/

#include "Interfaces/MecheInterface.hpp"

#include <Eigen/Core>

#include <iostream>
#include <cstdlib>

int main( int argc, const char* argv[] )
{

	if( argc < 2 )
	{
		std::cerr << " Please provide a problem data file " << std::endl ;
		std::cerr << " Syntax: " << argv[0] << " dataFile "
		          << " [ deterministic ] [ tol ] [ maxIters ] [ staticPb ] [ regul ] [ useInfNorm ] [algorithm] [cadouxIterations]"
		          << std::endl ;
		return 1 ;
	}

	bogus::MecheFrictionProblem mfp ;

	double * r = BOGUS_NULL_PTR(double) ;
	if( mfp.fromFile( argv[1], r ) )
	{

		bogus::MecheFrictionProblem::Options options ;

		options.maxThreads    = argc > 2 ? std::atoi( argv[2] ) : 0 ;
		if( options.maxThreads > 1 )
			options.gsColoring = true ;

		options.tolerance     = argc > 3 ? std::strtod( argv[3], NULL ) : 0 ;
		options.maxIters      = argc > 4 ? std::atoi( argv[4] ) : 0 ;

		const bool staticPb   = argc > 5 ? std::atoi( argv[5] ) : 0 ;
		double problemRegularization  = argc > 6 ? std::strtod( argv[6], NULL ) : 0 ;

		if( !staticPb ) {
			options.gsRegularization = problemRegularization ;
			problemRegularization = 0 ;
		}

		options.useInfinityNorm = argc > 7 ? std::atoi( argv[7] ) : 0 ;
		options.algorithm =   (bogus::MecheFrictionProblem::Algorithm)
		        ( argc > 8 ? std::atoi( argv[8] ) : 0 ) ;
		options.cadouxIters    = argc > 9 ? std::atoi( argv[9] ) : 0 ;

		mfp.solve( r, NULL, staticPb, problemRegularization, options ) ;
		std::cout << "Solver timer: " << mfp.lastSolveTime() << " seconds" << std::endl ;


		delete[] r ;

		return 0 ;
	}

	std::cerr << " Could not load " << argv[1] << std::endl ;

	return 1 ;

}

