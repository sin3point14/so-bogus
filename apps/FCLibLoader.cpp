
extern "C"
{
#include <fclib.h>
}

#include <Eigen/Sparse>

#include <Core/Block.impl.hpp>
#include <Core/Block.io.hpp>
#include <Core/BlockSolvers.impl.hpp>
#include <Interfaces/FrictionProblem.hpp>

#include <iostream>
#include <cstdlib>

void ackCurrentResidual( unsigned GSIter, double err )
{
	std::cout << " .. GS: " << GSIter << " ==> " << err << std::endl ;
}

template< unsigned Dimension, typename EigenDerived >
static double solve( const fclib_local* problem, const Eigen::SparseMatrixBase< EigenDerived >& ei_W,
                     Eigen::VectorXd &r, Eigen::VectorXd &u )
{
	bogus::DualFrictionProblem< Dimension > dual ;
	bogus::convert( ei_W, dual.W ) ;
	dual.W.cacheTranspose();

	dual.b = Eigen::VectorXd::Map( problem->q, problem->W->n ) ;
	dual.mu = problem->mu ;

	typename bogus::DualFrictionProblem< Dimension >::GaussSeidelType gs ;
	gs.setMaxIters( 1000 ) ;
	gs.setTol( 1.e-12 ) ;
	gs.callback().connect( &ackCurrentResidual );
	double res = dual.solveWith( gs, r.data() ) ;

	u = dual.W * r + dual.b ;
	
	return res ;
}

int main( int argc, const char* argv[] )
{

  if( argc < 2 )
  {
    std::cerr << " Please provide a problem data file " << std::endl ;
    return 1 ;
  }

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
//				  std::cout << ei_W << std::endl ;
			  }


			  double res = -1. ;
			  Eigen::VectorXd r, u ;
			  r.setZero( problem->W->n ) ;

			  if( problem->spacedim == 3 )
			  {
				  res = solve< 3u >( problem, ei_W, r, u ) ;
			  } else {
				  res = solve< 2u >( problem, ei_W, r, u ) ;
			  }
			  
			  std::cout << " => Res: " << res << std::endl ;
			  //			  std::cout << r.transpose() << std::endl ;

			  fclib_solution sol ;
			  sol.v = NULL ;
			  sol.u = u.data();
			  sol.r = r.data() ;
			  sol.v = NULL ;

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


  return 1 ;

}
