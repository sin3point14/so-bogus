#ifndef BOGUS_LOCAL_SOC_SOLVER_HPP
#define BOGUS_LOCAL_SOC_SOLVER_HPP

#include "../Utils/NumTraits.hpp"

namespace bogus {

template< unsigned Dimension, typename Scalar, bool DeSaxceCOV >
struct LocalSOCSolver
{

  typedef MatrixTraits< Dimension, Scalar > Traits ;
  typedef typename Traits::Vector Vector ;
  typedef typename Traits::Matrix Matrix ;

  static Scalar solve(
          const typename Traits::Matrix &A,
          const typename Traits::Vector &b,
          typename Traits::Vector &x,
          const Scalar mu, const Scalar tol
          ) ;

} ;

}

#endif
