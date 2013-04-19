#ifndef BOGUS_LOCAL_SOC_SOLVER_HPP
#define BOGUS_LOCAL_SOC_SOLVER_HPP

#include "../SecondOrder.fwd.hpp"
#include "../Utils/NumTraits.hpp"
#include "../Utils/EigenMatrixTraits.hpp"

namespace bogus {


template< unsigned Dimension, typename Scalar, bool DeSaxceCOV,
		   local_soc_solver::Strategy Strat = local_soc_solver::Hybrid  >
struct LocalSOCSolver
{

  typedef MatrixTraits< Dimension, Scalar > Traits ;
  typedef typename Traits::Vector Vector ;
  typedef typename Traits::Matrix Matrix ;

  static Scalar solve(
		  const typename Traits::Matrix &A,
		  const typename Traits::Vector &b,
		  typename Traits::Vector &x,
		  const Scalar mu, const Scalar tol,
		  const Scalar scaling = 1
		  ) ;

} ;

}

#endif
