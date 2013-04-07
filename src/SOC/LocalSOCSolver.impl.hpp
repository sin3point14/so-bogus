#ifndef BOGUS_LOCAL_SOC_SOLVER_IMPL_HPP
#define BOGUS_LOCAL_SOC_SOLVER_IMPL_HPP

#include "LocalSOCSolver.hpp"

namespace bogus {

template< unsigned Dimension, typename Scalar, bool DeSaxceCOV >
Scalar LocalSOCSolver< Dimension, Scalar, DeSaxceCOV >::solve(
          const typename Traits::Matrix &A,
          const typename Traits::Vector &b,
          typename Traits::Vector &x,
          const Scalar mu, const Scalar tol
          )
{
    (void) A ;
    (void) b ;
    (void) x ;
    (void) mu ;
    (void) tol ;
    return -1 ;
}

}

#endif
