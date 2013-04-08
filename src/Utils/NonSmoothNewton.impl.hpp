#ifndef BOGUS_NS_NETWON_IMPL_HPP
#define BOGUS_NS_NETWON_IMPL_HPP

#include "NumTraits.hpp"
#include "NonSmoothNewton.hpp"

#include <Eigen/LU>

#include <iostream>

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

template < typename NSFunction >
typename NonSmoothNewton< NSFunction >::Scalar NonSmoothNewton<NSFunction>::solve(
  Vector& x ) const
{
  static const Scalar sigma2 = 1.e-4 ;
  static const Scalar alpha = .5 ;

  typedef LinearSolver< Traits::dimension, Scalar > LinearSolver ;

  Vector F ;
  m_func.compute( x, F ) ;
  const Scalar Phi_init = F.squaredNorm() ;

  if( Phi_init < m_tol ) return Phi_init ;

  Scalar Phi_best ;
  Vector x_best = Vector::Zero() ;

  m_func.compute( x_best, F ) ;
  const Scalar Phi_zero = F.squaredNorm() ;

  if( Phi_zero < Phi_init ) {
	Phi_best = Phi_zero ;
	x = x_best ;
  } else {
	Phi_best = Phi_init ;
	x_best = x ;
  }

  if( Phi_zero < m_tol ) return Phi_zero ;

  Matrix dF_dx ;
  Vector dPhi_dx, dx ;

  for( unsigned iter = 0 ; iter < m_maxIters ; ++iter )
  {
	m_func.computeJacobian( x, F, dF_dx ) ;
	const Scalar Phi = F.squaredNorm() ;

	if( Phi < m_tol ) return Phi ;
	if( Phi < Phi_best ) {
	  Phi_best = Phi ;
	  x_best = x ;
	}
	dPhi_dx = dF_dx.transpose() * x ;

	dx = - LinearSolver::solve( dF_dx, F ) ;
	const Scalar proj = dx.dot( dPhi_dx ) ;

	if( proj > 0 || proj * proj < sigma2 * dx.squaredNorm() * dPhi_dx.squaredNorm() )
	{
		dx *= alpha ;
	}

	x += dx ;

  }

  x = x_best ;
  return Phi_best ;

}

}


#endif

