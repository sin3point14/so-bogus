#ifndef BOGUS_LOCAL_SOC_SOLVER_IMPL_HPP
#define BOGUS_LOCAL_SOC_SOLVER_IMPL_HPP

#include "LocalSOCSolver.hpp"
#include "FischerBurmeister.hpp"
#include "FischerBurmeister.impl.hpp"

#include "../Utils/NonSmoothNewton.hpp"
#include "../Utils/NonSmoothNewton.impl.hpp"

#include "../Utils/Polynomial.hpp"
#include "../Utils/Polynomial.impl.hpp"

#include "../Utils/EigenLinearSolvers.hpp"

#define BOGUS_PURE_ENUMERATIVE

namespace bogus {


template< unsigned Dimension, typename Scalar, bool DeSaxceCOV, local_soc_solver::Strategy Strat >
Scalar LocalSOCSolver< Dimension, Scalar, DeSaxceCOV, Strat >::solve(
		  const typename Traits::Matrix &A,
		  const typename Traits::Vector &b,
		  typename Traits::Vector &x,
		  const Scalar mu, const Scalar tol, const Scalar scaling
		  )
{
	typedef FischerBurmeister< Dimension, Scalar, DeSaxceCOV > FBFunc ;
	FBFunc fb( mu, A, b, scaling ) ;
	NonSmoothNewton< FBFunc > nsNewton( fb, tol )  ;

	return nsNewton.solve( x ) ;
}

// Specialization for 3D Coulomb friction -- hybrid solver
template< typename Scalar, local_soc_solver::Strategy Strat  >
struct LocalSOCSolver< 3, Scalar, true, Strat >
{
  enum { Dimension = 3 } ;

  typedef MatrixTraits< Dimension, Scalar > Traits ;
  typedef typename Traits::Vector Vector ;
  typedef typename Traits::Matrix Matrix ;

  static Scalar solve(
		  const typename Traits::Matrix &A,
		  const typename Traits::Vector &b,
		  typename Traits::Vector &x,
		  const Scalar mu, const Scalar tol, const Scalar scaling
		  )
  {
	  // see [Daviet et al 2011], Appendix B.2

	  // Newton solver
	  typedef FischerBurmeister< Dimension, Scalar, true > FBFunc ;
	  FBFunc fb( mu, A, b, scaling ) ;
	  NonSmoothNewton< FBFunc > nsNewton( fb, tol )  ;

	  if( Strat == local_soc_solver::PureNewton )
	  {
		  return nsNewton.solve( x ) ;
	  }

	  if( Traits::np(b) >= 0. )
	  {
		  // Take-off case
		  x.setZero() ;
		  return 0. ;
	  }
	  if( NumTraits< Scalar >::isZero( mu ) )
	  {
		  //Frictionless case
		  if( A(0,0) < NumTraits< Scalar >::epsilon() )
		  {
			  // Degenerate problem
			  x.setZero() ;
			  return b[0]*b[0] ;
		  } else {
			  Traits::tp( x ).setZero() ;
			  Traits::np( x ) = - Traits::np( b ) / A(0,0);
			  return 0. ;
		  }
	  }

	  double res = 0. ;

	  if( Strat == local_soc_solver::Hybrid )
	  {
		  res = nsNewton.solve( x ) ;
		  if( res < tol ) return res ;
	  }

	  // Continuing enumerative fallback

	  Vector x0 = x ;
	  x = DenseLU< Scalar, Dimension >( A ).solve( -b ) ;
	  if( mu * Traits::np( x ) >= Traits::tp( x ).norm() )
	  {
		  // Sticking case
		  return 0. ;
	  }

	 // Sliding case
	 if( !solveSliding( A, b, x, mu ) )
		 x = x0 ;

	 // Refinement of final solution
	 if( Strat == local_soc_solver::RevHybrid ) {
		 res = nsNewton.solve( x ) ;
	 } else if( Strat == local_soc_solver::Hybrid  ) {
		 const double refinedRes = nsNewton.solve( x ) ;
		 if( refinedRes <= res )
			 return refinedRes ;

		 //This can happen if the quartic solver returned a very bad value, like an
		 // unreastically huge alpha
		 x = x0 ;
	 }

	 return res ;
  }

  static bool solveSliding(
		  const typename Traits::Matrix &W,
		  const typename Traits::Vector &b,
		  typename Traits::Vector &r,
		  const Scalar mu
		  )
  {
	  // see [Daviet et al 2011], Appendix B.1

	 typedef Eigen::Matrix< Scalar, 2, 1 > Vec2 ;
	 typedef Eigen::Matrix< Scalar, 2, 2 > Mat2 ;

	 const Scalar wN = W(0,0) ;
	 if( wN < NumTraits< Scalar >::epsilon() )
		 return false ; // Could we do something better ?

	 const Vec2 wT = W.template block< 2, 1 >( 1, 0 ) ;
	 const Mat2 WT = W.template block< 2, 2 >( 1, 1 ) ;

	 const Scalar bN = Traits::np( b );
	 const Vec2 bT = Traits::tp( b ) ;

	 const Vec2 wT_wN = wT/wN ;
	 const Mat2 Wbar = WT - wT_wN * wT.transpose() ;
	 const Vec2 bbar = bT/bN - wT_wN ;

	 const Scalar A = Wbar.trace() - wT.dot( bbar ) ;
	 const Vec2   B ( Wbar(1,1)*bbar[0] - Wbar(1,0)*bbar[1],
					  Wbar(0,0)*bbar[1] - Wbar(0,1)*bbar[0] ) ;
	 const Scalar C = Wbar.determinant() - wT.dot( B ) ;
	 const Scalar D = wN*wN / ( mu*mu ) ;

	 const Scalar coeffs[4] = {
		 C*C - D * B.squaredNorm(),
		 2*( C*A - D * bbar.dot( B ) ),
		 2*C + A*A - D * bbar.squaredNorm(),
		 2*A
	 } ;

	 Scalar roots[4] ;
	 const unsigned nRoots =
			 polynomial::getRealRoots( coeffs, roots, polynomial::StrictlyPositiveRoots ) ;

	 if( 0 == nRoots ) return false ;

	 Scalar alpha = roots[0] ;
	 //Get the minimal one, this is as good an heuristic as any
	 for ( unsigned i = 1 ; i != nRoots ; ++ i )
	 {
		 if( roots[i] < alpha ) alpha = roots[i] ;
	 }

//	 std::cout << "Found " << alpha << std::endl ;

	 const Mat2 M = Wbar + alpha * Mat2::Identity() ;
	 Traits::tp( r ) = - bN * M.fullPivLu().solve( bbar ) ;
	 Traits::np( r ) = Traits::tp( r ).norm() / mu ;

	 return true ;
  }

} ;

// Specialization for 3D SOC complementarity -- hybrid solver
template< typename Scalar, local_soc_solver::Strategy Strat  >
struct LocalSOCSolver< 3, Scalar, false, Strat >
{
  enum { Dimension = 3 } ;

  typedef MatrixTraits< Dimension, Scalar > Traits ;
  typedef typename Traits::Vector Vector ;
  typedef typename Traits::Matrix Matrix ;

  static Scalar solve(
		  const typename Traits::Matrix &A,
		  const typename Traits::Vector &b,
		  typename Traits::Vector &x,
		  const Scalar mu, const Scalar tol, const Scalar scaling
		  )
  {
	  // Newton solver
	  typedef FischerBurmeister< Dimension, Scalar, false > FBFunc ;
	  FBFunc fb( mu, A, b, scaling ) ;
	  NonSmoothNewton< FBFunc > nsNewton( fb, tol )  ;

	  if( Strat == local_soc_solver::PureNewton )
	  {
		  return nsNewton.solve( x ) ;
	  }

	  if( Traits::np(b) >= mu * Traits::tp(b).norm() )
	  {
		  // b in dual cone
		  x.setZero() ;
		  return 0. ;
	  }

	  if( NumTraits< Scalar >::isZero( mu ) )
	  {
		  //Frictionless case
		  if( A(0,0) < NumTraits< Scalar >::epsilon() )
		  {
			  // Degenerate problem
			  x.setZero() ;
			  return b[0]*b[0] ;
		  } else {
			  Traits::tp( x ).setZero() ;
			  Traits::np( x ) = - Traits::np( b ) / A(0,0);
			  return 0. ;
		  }
	  }

	  double res = 0. ;

	  if( Strat == local_soc_solver::Hybrid )
	  {
		  res = nsNewton.solve( x ) ;
		  if( res < tol ) return res ;
	  }

	  // Continuing enumerative fallback

	  Vector x0 = x ;
	  x = DenseLU< Scalar, Dimension >( A ).solve( -b ) ;
	  if( mu * Traits::np( x ) >= Traits::tp( x ).norm() )
	  {
		  // Sticking case
		  return 0. ;
	  }

	 // Sliding case
	 if( !solveOrthogonality( A, b, x, mu ) )
		 x = x0 ;

	 // Refinement of final solution
	 if( Strat == local_soc_solver::RevHybrid ) {
		 res = nsNewton.solve( x ) ;
	 } else if( Strat == local_soc_solver::Hybrid  ) {
		 const double refinedRes = nsNewton.solve( x ) ;
		 if( refinedRes <= res )
			 return refinedRes ;

		 //This can happen if the quartic solver returned a very bad value
		 x = x0 ;
	 }

	 return res ;
  }

  static bool solveOrthogonality(
		  const typename Traits::Matrix &W,
		  const typename Traits::Vector &b,
		  typename Traits::Vector &r,
		  const Scalar mu
		  )
  {
	  // see doc/sage/poySOC.sage

	 typedef Eigen::Matrix< Scalar, 2, 1 > Vec2 ;
	 typedef Eigen::Matrix< Scalar, 2, 2 > Mat2 ;

	 const Scalar wN = W(0,0) ;
	 if( wN < NumTraits< Scalar >::epsilon() )
		 return false ; // Could we do something better ?

	 const Scalar A = (b[1]*W(1,2) - b[2]*W(1,1))*mu*mu + b[0]*W(0,2) - b[2]*W(0,0) ;
	 const Scalar B = (b[0]*W(1,2) + b[1]*W(0,2) - 2*b[2]*W(0,1))*mu ;
	 const Scalar C = 2*( (b[1]*W(2,2) - b[2]*W(1,2))*mu*mu - b[0]*W(0,1) + b[1]*W(0,0) );
	 const Scalar D = 2*( b[0]*(W(1,1) - W(2,2)) - b[1]*W(0,1) + b[2]*W(0,2))*mu ;
	 const Scalar E = -6*(b[0]*W(1,2) - b[1]*W(0,2))*mu ;

	 Scalar coeffs[5] ;
	 coeffs[0] =  A + B ;
	 coeffs[1] =  C - D ;
	 coeffs[2] =  E ;
	 coeffs[3] =  C + D ;
	 coeffs[4] =  B - A ;

	 Scalar roots[4] ;
	 const unsigned nRoots =
			 bogus::polynomial::getRealRoots( coeffs, roots, bogus::polynomial::AllRoots ) ;

	 bool found = false ;

	 for ( unsigned i = 0 ; i != nRoots ; ++ i )
	 {
		 Scalar t = roots[i] ;

		 const Scalar CT = ( 1 - t*t ) / ( 1 + t*t ) ;
		 const Scalar ST = 2*t / ( 1 + t*t ) ;

		 const Eigen::Vector3d dir ( 1, mu*CT, mu*ST ) ;

		 const Scalar den = ( mu * W.col(1) + CT * W.col( 0 )).dot( dir ) ;
		 if( bogus::NumTraits< Scalar >::isZero( den ) )
			 continue ;

		 const Scalar rN = -(CT*b[0] + b[1]*mu)/den ;
		 if( rN <= 0 )
			 continue ;

		 r = rN * dir ;
		 const Scalar uN = W.col(0).dot(r) + b[0] ;

		 if( uN > 0 )
		 {
			 found = true ;
			 break ;
		 }
	 }

	 return found ;
  }

} ;

}

#endif
