#ifndef BOGUS_LOCAL_SOC_SOLVER_IMPL_HPP
#define BOGUS_LOCAL_SOC_SOLVER_IMPL_HPP

#include "LocalSOCSolver.hpp"
#include "FischerBurmeister.hpp"

#include "../Utils/NonSmoothNewton.hpp"
#include "../Utils/NonSmoothNewton.impl.hpp"

#include "../Utils/Polynomial.hpp"
#include "../Utils/Polynomial.impl.hpp"

#define BOGUS_PURE_ENUMERATIVE

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

// Specialization for 3D Coulomb friction -- hybrid solver
template< typename Scalar >
struct LocalSOCSolver< 3, Scalar, true >
{
  enum { Dimension = 3 } ;

  typedef MatrixTraits< Dimension, Scalar > Traits ;
  typedef typename Traits::Vector Vector ;
  typedef typename Traits::Matrix Matrix ;

  static Scalar solve(
		  const typename Traits::Matrix &A,
		  const typename Traits::Vector &b,
		  typename Traits::Vector &x,
		  const Scalar mu, const Scalar tol
		  )
  {
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

	 // Newton solver
	 typedef FischerBurmeister< Dimension, Scalar, true > FBFunc ;
	 FBFunc fb( mu, A, b ) ;
	 NonSmoothNewton< FBFunc > nsNewton( fb, tol )  ;

#ifndef BOGUS_PURE_ENUMERATIVE
	 const double newtonRes = nsNewton.solve( x ) ;
	 if( newtonRes < tol ) return newtonRes ;
#endif

	 // Continuing enumerative fallback

	 Vector u = -b, x0 = x ;
	 x = LinearSolver< Dimension, Scalar >::solve( A, u ) ;
	 if( mu * Traits::np( x ) >= Traits::tp( x ).norm() )
	 {
		 // Sticking case
		 return 0. ;
	 }

	 // Sliding case
	 if( !solveSliding( A, b, x, mu ) )
		 x = x0 ;

#ifdef BOGUS_PURE_ENUMERATIVE
	 return 0. ;
#else
	 // Refinement of final solution
	 const double refinedRes = nsNewton.solve( x ) ;
	 if( refinedRes > newtonRes ) {
		 //This can happen if the quartic solver returned a very bad value, like an
		 // unreastically huge alpha
		 x = x0 ;
		 return newtonRes ;
	 }

	 return refinedRes ;
#endif
  }

  static bool solveSliding(
		  const typename Traits::Matrix &W,
		  const typename Traits::Vector &b,
		  typename Traits::Vector &r,
		  const Scalar mu
		  )
  {
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

}

#endif
