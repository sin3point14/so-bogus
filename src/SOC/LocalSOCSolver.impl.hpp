#ifndef BOGUS_LOCAL_SOC_SOLVER_IMPL_HPP
#define BOGUS_LOCAL_SOC_SOLVER_IMPL_HPP

#include "LocalSOCSolver.hpp"
#include "FischerBurmeister.hpp"

#include "../Utils/NonSmoothNewton.hpp"
#include "../Utils/NonSmoothNewton.impl.hpp"

namespace bogus {

template< unsigned Dimension, typename Scalar, bool DeSaxceCOV >
Scalar LocalSOCSolver< Dimension, Scalar, DeSaxceCOV >::solve(
          const typename Traits::Matrix &A,
          const typename Traits::Vector &b,
          typename Traits::Vector &x,
          const Scalar mu, const Scalar tol
          )
{
    typedef FischerBurmeister< Dimension, Scalar, DeSaxceCOV > FBFunc ;
    FBFunc fb( mu, A, b ) ;
    NonSmoothNewton< FBFunc > nsNewton( fb, tol )  ;

    return nsNewton.solve( x ) ;
}

}

#endif
