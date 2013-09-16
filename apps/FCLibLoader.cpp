
/*
 * Any copyright is dedicated to the Public Domain.
 * http://creativecommons.org/publicdomain/zero/1.0/
*/

extern "C"
{
#include <fclib.h>
}

#include <Core/Block.impl.hpp>
#include <Core/Block.io.hpp>
#include <Core/BlockSolvers/ProjectedGradient.impl.hpp>
#include <Core/BlockSolvers/GaussSeidel.impl.hpp>
#include <Interfaces/FrictionProblem.hpp>

#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sys/stat.h>

//#define USE_PG_FOR_CADOUX


static const char* g_meth  = 0;

void ackCurrentResidual( unsigned GSIter, double err )
{
	std::cout << " .. " << g_meth << ": " <<  GSIter << " ==> " << err << std::endl ;
}

template< unsigned Dimension, typename EigenDerived >
static double solve( const fclib_local* problem, const Eigen::SparseMatrixBase< EigenDerived >& ei_W,
					 Eigen::VectorXd &r, Eigen::VectorXd &u, const bool useCadoux )
{
	bogus::DualFrictionProblem< Dimension > dual ;
	bogus::convert( ei_W, dual.W, Dimension, Dimension ) ;
	dual.W.prune( 1.e-12 ) ;
	dual.W.cacheTranspose();

	dual.b = Eigen::VectorXd::Map( problem->q, problem->W->n ) ;
	dual.mu = Eigen::VectorXd::Map( problem->mu, problem->W->n/Dimension ) ;

	typename bogus::DualFrictionProblem< Dimension >::GaussSeidelType gs ;
	gs.setTol( 1.e-12 ) ;
	gs.setAutoRegularization( 1.e-5 ) ;

	bogus::Signal< unsigned, double > callback ;
	callback.connect( &ackCurrentResidual );

	double res = -1 ;

	if( useCadoux )
	{
		g_meth = "Cadoux" ;
#ifdef USE_PG_FOR_CADOUX
		typename bogus::DualFrictionProblem< Dimension >::ProjectedGradientType pg ;
		pg.setTol( 1.e-12 ) ;
		pg.setMaxIters( 50 ) ;
		res = dual.solveCadoux( pg, r.data(), 500, &callback ) ;
#else
		gs.setMaxIters( 100 ) ;
		res = dual.solveCadoux( gs, r.data(), 500, &callback ) ;
#endif

	} else {
		g_meth = "GS" ;
		gs.setMaxIters( 1000 ) ;
		gs.callback().connect( callback );
		res = dual.solveWith( gs, r.data() ) ;
	}

	u = dual.W * r + dual.b ;

	return res ;
}

int main( int argc, const char* argv[] )
{
#ifdef BOGUS_WITH_EIGEN_STABLE_SPARSE_API

  if( argc < 2 )
  {
	std::cerr << " Please provide a problem data file " << std::endl ;
	return 1 ;
  }

  const bool useCadoux = argc > 2 && 0 == strncmp( argv[2], "-c", 3 ) ;

  const char* file = argv[1] ;

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

			  /*
			  Eigen::MappedSparseMatrix< double, Eigen::RowMajor > ei_W
					  ( problem->W->m, problem->W->n,
						problem->W->nzmax, problem->W->p, problem->W->i,
						problem->W->x ) ;
			  */
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

			  if( problem->spacedim == 3 )
			  {
				  res = solve< 3u >( problem, ei_W, r, u, useCadoux ) ;
			  } else {
				  res = solve< 2u >( problem, ei_W, r, u, useCadoux ) ;
			  }

			  std::cout << " => Res: " << res << std::endl ;
			  //			  std::cout << r.transpose() << std::endl ;

			  fclib_solution sol ;
			  sol.v = NULL ;
			  sol.u = u.data();
			  sol.r = r.data() ;
			  sol.l = NULL ;


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

//			  std::cout << fclib_merit_local( problem, MERIT_2, &sol ) << std::endl ;

			  // Just for fun
			  // Eigen::SparseMatrix< double > ei_W_back ;
			  // bogus::convert( dual.W, ei_W_back ) ;
			  // std::cout << ei_W_back << std::endl ;


		  }
	  }

	  //solution = fclib_read_solution (file);
	  //guesses = fclib_read_guesses (file, &n_guesses);
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
