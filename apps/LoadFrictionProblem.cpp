
#include "Interfaces/MecheInterface.hpp"

#include <Eigen/Core>

#include <iostream>
#include <cstdlib>

int main( int argc, const char* argv[] )
{

  if( argc < 2 )
  {
    std::cerr << " Please provide a problem data file " << std::endl ;
    return 1 ;
  }

  bogus::MecheFrictionProblem mfp ;

  double * r = NULL ;
  if( mfp.fromFile( argv[1], r ) )
  {

    const int deterministic = argc > 2 ? std::atoi( argv[2] ) : 0 ;
    const double tol        = argc > 3 ? std::strtod( argv[3], NULL ) : 0 ;
    const int maxIters      = argc > 4 ? std::atoi( argv[4] ) : 0 ;
    const int staticPb      = argc > 5 ? std::atoi( argv[5] ) : 0 ;

    mfp.solve( r, NULL, deterministic, tol, maxIters, staticPb ) ;

    delete[] r ;
  
    return 0 ;
  } 
    
  std::cerr << " Could not load " << argv[1] << std::endl ;

  return 1 ;

}

