#ifndef BOGUS_FISCHER_BURMEISTER_HPP
#define BOGUS_FISCHER_BURMEISTER_HPP

#include "../Utils/NumTraits.hpp"

namespace bogus {

template< unsigned Dimension, typename Scalar >
struct FBBaseFunction
{
	typedef MatrixTraits< Dimension, Scalar > Traits ;
	typedef typename Traits::Vector Vector ;
	typedef typename Traits::Matrix Matrix ;

	static void compute( const Scalar mu, const Vector& x, const Vector& y, Vector& fb ) ;

	static void computeJacobian(
                const Scalar mu, const Vector& x, const Vector& y, 
	        Vector& fb, Matrix& dFb_dx, Matrix& dFb_dy ) ;

private:
	template <bool JacobianAsWell >
	static void compute(
			const Vector& x, const Vector& y, Vector& fb,
			Matrix& dFb_dx, Matrix& dFb_dy ) ;

} ;

template< unsigned Dimension, typename Scalar, bool DeSaxceCOV >
struct FischerBurmeister
{
  typedef MatrixTraits< Dimension, Scalar > Traits ;
  typedef FBBaseFunction< Dimension, Scalar > BaseFunction ;
  typedef typename Traits::Vector Vector ;
  typedef typename Traits::Matrix Matrix ;

  static void compute( const Scalar mu, const Vector& x, const Vector& y, Vector& fb ) ;
	
  
  static void computeJacobian(
      const Scalar mu, const Vector& x, const Matrix& A, const Vector& b, 
      Vector& fb, Matrix& dFb_dx ) ;

} ;



}

#endif
