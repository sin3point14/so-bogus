
extern "C"
{
#include <fclib.h>
}

#include <Eigen/Sparse>

#include <iostream>
#include <cstdlib>

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

	  if( problem->W && !problem->V )
	  {
		  std::cout << " Pure " << d << "D Coulomb friction problem with "
		  			<< n << " contacts " << std::endl ;

		  if( problem->W->nz == -2 )
		  {
			  std::cout << " Compressed row storage " << problem->W->nz << std::endl ;

			  Eigen::SparseMatrix< double, Eigen::RowMajor > ei_W ;
			  ei_W.resize( problem->W->m, problem->W->n );
			  ei_W.resizeNonZeros( problem->W->nzmax ) ;

			  memcpy( ei_W.outerIndexPtr(), problem->W->p, problem->W->m+1 * sizeof( int ) ) ;
			  memcpy( ei_W.innerIndexPtr(), problem->W->i, problem->W->nzmax * sizeof( int ) ) ;
			  memcpy( ei_W.valuePtr(), problem->W->x, problem->W->nzmax * sizeof( double ) ) ;



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
