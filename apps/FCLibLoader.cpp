
/*
 * Any copyright is dedicated to the Public Domain.
 * http://creativecommons.org/publicdomain/zero/1.0/
*/

#include "FCLibSolver.hpp"
#include <bogus/Extra/SOC/SOCLaw.impl.hpp>

extern "C"
{
#include <fclib.h>
}

#include <cstdlib>

#include <iostream>
#include <string>
#include <fstream>
#include <sys/stat.h>



void usage(const char* name)
{
	std::cout << "Usage: " << name << " dataFile [options] \n "
	          << "Options: \n"
	          << " -a int \t algorithm id in [0,1] \n"
	          << " -T int \t max number of threads (or 0 for OMP_MAX_THREADS) \n"
	          << " -m int \t max number of iterations \n"
	          << " -c int \t max number of Cadoux fixed-point iterations (0 means direct solving) \n"
	          << " -t real \t tolerance \n"
	          << " -i bool \t if true, use infinity norm instead of l2 \n"
	          << " -r real \t proximal regularization for GS algorithm \n"
	          << " -g int \t projected gradiant variant in [0,4] \n"
	          << " -v bool \t if true, print residual at each iteration \n"
	          << std::endl ;

}

namespace bogus{
namespace fclib {

// reimplementation of fcmerc. merti function to avoid the accuracy loss
// from the difference of ordering in the matrix-vector product W*r
template < unsigned Dim >
static double merit1( const Eigen::VectorXd&r, const Eigen::VectorXd&u,
                      const Eigen::VectorXd&b,
                      const double *mu )
{
	const int n = r.rows() / Dim ;
	typedef SOCLaw< Dim, double, false > SOC ;
	SOC socLaw( n, mu ) ;

	double err = 0 ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for reduction( +:err )
#endif
	for( int i = 0 ; i < n ; ++i ) {
		typename SOC::Traits::Vector x = r.segment<Dim>( Dim*i ) - u.segment<Dim>( Dim*i ) ;
		x[0] -= mu[i] * u.segment<Dim-1>( Dim*i + 1 ).norm() ;
		socLaw.projectOnConstraint( i, x ) ;
		const double lerr = (x - r.segment<Dim>( Dim*i )).squaredNorm() ;
		err += lerr ;
	}

	// fcmer use sqrt(b.norm())
	return std::sqrt(err) / ( 1 + b.norm() ) ;
}

template< unsigned Dimension, typename EigenDerived >
static double solve( const fclib_local* problem,
                     const Eigen::SparseMatrixBase< EigenDerived >& ei_W,
                     Options& options,
                     Eigen::VectorXd &r, Eigen::VectorXd &u,
                     Stats& stats )
{
	bogus::DualFrictionProblem< Dimension > dual ;
	bogus::convert( ei_W, dual.W, Dimension, Dimension ) ;

	dual.W.prune( options.tolerance ) ;
	dual.W.cacheTranspose();

	dual.b = Eigen::VectorXd::Map( problem->q, problem->W->n ) ;
	dual.mu = Eigen::VectorXd::Map( problem->mu, problem->W->n/Dimension ) ;

	options.tolerance *= ( 1 + dual.b.norm() )/problem->W->n ;
	solveDual( dual, options, r, u, stats ) ;
	stats.error *= problem->W->n/( 1 + dual.b.norm() ) ;

	return merit1< Dimension >( r, u, dual.b, problem->mu ) ;

}


} //fclib
} //bogus


int main( int argc, const char* argv[] )
{
#ifdef BOGUS_WITH_EIGEN_STABLE_SPARSE_API

	bogus::fclib::Options options ;
	bool verbose = false ;

	const char* file = BOGUS_NULL_PTR(const char) ;

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
			case 'i':
				if( ++i == argc ) break ;
				options.useInfinityNorm = (bool) std::atoi( argv[i] ) ;
				break ;
			case 'a':
				if( ++i == argc ) break ;
				options.algorithm =   (bogus::fclib::Options::Algorithm)
				        std::atoi( argv[i] ) ;
				break ;
			case 'g':
				if( ++i == argc ) break ;
				options.pgVariant =   (bogus::projected_gradient::Variant)
				        std::atoi( argv[i] ) ;
				break ;
			case 'r':
				if( ++i == argc ) break ;
				options.gsRegularization = std::strtod( argv[i], NULL ) ;
				break ;
			case 'v':
				if( ++i == argc ) break ;
				verbose = (bool) std::atoi( argv[i] ) ;
				break ;
			}
		} else {
			file = argv[i] ;
		}
	}

	if( options.maxThreads > 1 )
		options.gsColoring = true ;

	if( !file )
	{
		std::cerr << " Please provide a problem data file " << std::endl ;
		usage(argv[0]) ;
		return 1 ;
	}

	struct fclib_local *problem = NULL ;
	struct fclib_solution *solution = NULL ;
	struct fclib_solution *guesses = NULL ;
	int n_guesses = 0 ;

	problem = fclib_read_local (file);

	if( problem )
	{
		std::cout << "Successfully loaded problem " << file << std::endl ;
		std::cout << " Name: " << problem->info->title << std::endl ;
		std::cout << " Descr: " << problem->info->description << std::endl ;
		std::cout << " Info: " << problem->info->math_info << std::endl ;

		if( problem->V ) std::cout << " Has V " << std::endl ;
		if( problem->R ) std::cout << " Has R " << std::endl ;

		const unsigned d = problem->spacedim  ;
		const unsigned n = problem->W->m / d  ;

		if( problem->W && !problem->V && !problem->R )
		{
			std::cout << " Pure " << d << "D Coulomb friction problem with "
			          << n << " contacts " << std::endl ;

			if( problem->W->nz == -2 )
			{
				std::cout << " Compressed row storage " << problem->W->m << " / " << problem->W->n << " / " << problem->W->nzmax << std::endl ;

				Eigen::SparseMatrix< double, Eigen::RowMajor > ei_W ;
				ei_W.resize( problem->W->m, problem->W->n );
				ei_W.resizeNonZeros( problem->W->nzmax ) ;

				memcpy( ei_W.outerIndexPtr(), problem->W->p, ( problem->W->m+1 ) * sizeof( int ) ) ;
				memcpy( ei_W.innerIndexPtr(), problem->W->i, problem->W->nzmax * sizeof( int ) ) ;
				memcpy( ei_W.valuePtr(), problem->W->x, problem->W->nzmax * sizeof( double ) ) ;

				if( 0 == problem->W->p[ problem->W->m ] )
				{
					std::cout << " /!\\ Malformed spase matrix ; entering repair mode " << std::endl ;
					assert( problem->W->nzmax == problem->W->n * problem->W->m ) ;

					ei_W.outerIndexPtr()[ 0 ] = 0 ;
					for( int row = 0 ; row < problem->W->m ; ++row )
					{
						const int start = ei_W.outerIndexPtr()[ row ] ;
						ei_W.outerIndexPtr()[ row+1 ] = start + problem->W->n ;
						for( int col = 0 ; col < problem->W->n ; ++col )
						{
							ei_W.innerIndexPtr()[ start + col ] = col ;
						}
					}

				}

				double res = -1. ;
				Eigen::VectorXd r, u ;
				r.setZero( problem->W->n ) ;

				bogus::fclib::Stats stats (verbose) ;

				if( problem->spacedim == 3 )
				{
					res = bogus::fclib::solve< 3u >( problem, ei_W, options, r, u, stats ) ;
				} else {
					res = bogus::fclib::solve< 2u >( problem, ei_W, options, r, u, stats ) ;
				}

				std::cout << " => Err: \t" << stats.error << std::endl ;
				std::cout << " => Iters: \t"<< stats.nIters << std::endl ;
				std::cout << " => Time: \t"<< stats.time << std::endl ;

				fclib_solution sol ;
				sol.v = NULL ;
				sol.u = u.data();
				sol.r = r.data() ;
				sol.l = NULL ;

				std::cout << " => .. FCLib Merit1: " << fclib_merit_local( problem, MERIT_1, &sol ) << std::endl ;
				std::cout << " => ..  True Merit1: " << res << std::endl ;

				//

				std::string fname ( file ) ;
				const std::size_t dir_pos = fname.rfind( '/' ) ;
				const std::string dir_name = fname.substr( 0, dir_pos ) + "/fix" ;

				struct stat info ;
				int exists = stat( dir_name.c_str(), &info ) ;

				if( exists == 0 && 0 != ( info.st_mode & S_IFDIR ) )
				{
					std::string outfname = dir_name + "/" + fname.substr(dir_pos+1) ;

					if( ! std::ifstream( outfname.c_str() ) )
					{

						std::cout << "Re-writing problem to " << outfname << std::endl ;

						ei_W.prune(1.e-12) ;
						problem->W->nzmax = ei_W.nonZeros() ;
						memcpy( problem->W->p, ei_W.outerIndexPtr(), ( problem->W->m+1 ) * sizeof( int ) ) ;
						memcpy( problem->W->i, ei_W.innerIndexPtr(), problem->W->nzmax * sizeof( int ) ) ;
						memcpy( problem->W->x, ei_W.valuePtr(), problem->W->nzmax * sizeof( double ) ) ;

						if (fclib_write_local (problem, outfname.c_str()))
						{
							if (fclib_write_solution (&sol, outfname.c_str() ))
							{
								std::cout << "Ok." << std::endl ;
							}
						}
					}
				}

			}
		}

	}

	fclib_delete_local (problem);
	if( solution ) fclib_delete_solutions (solution, 1);
	if( guesses ) fclib_delete_solutions (guesses, n_guesses);


	return 0 ;

#else
	(void) argc, (void) argv ;

	std::cerr<< "FCLibLoader requires Eigen's stable sparse API ( >= 3.1 )" << std::endl ;
	return -1 ;
#endif


}
