#ifndef BOGUS_LINEAR_SOLVER_HPP
#define BOGUS_LINEAR_SOLVER_HPP

#include "NumTraits.hpp"

#include <Eigen/LU>

namespace bogus {

template < unsigned Dimension, typename Scalar >
struct LinearSolver
{
  typedef MatrixTraits< Dimension, Scalar > Traits ;
  typedef Eigen::FullPivLU< typename Traits::Matrix > FactType ;
  typedef Eigen::internal::solve_retval< FactType, typename Traits::Vector > ReturnType ; // Why did I decide not to use c++11, again ?

  static typename Traits::Vector solve( const typename Traits::Matrix &A, const typename Traits::Vector &b )
  {
	return A.fullPivLu().solve( b ) ;
  }
} ;

}

#endif
