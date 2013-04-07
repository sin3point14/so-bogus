#ifndef BOGUS_FISCHER_BURMEISTER_IMPL_HPP
#define BOGUS_FISCHER_BURMEISTER_IMPL_HPP

#include "FischerBurmeister.hpp"

namespace bogus {

template< unsigned Dimension, typename Scalar >
void FBBaseFunction< Dimension, Scalar >::compute(
	 const Scalar mu, const Vector& x, const Vector& y, Vector& fb )
{
	static Matrix a,b ; //unused

	if( NumTraits< Scalar >::isZero( mu ) )
	{
		Traits::np( fb ) = x[0] + y[0] - std::sqrt( x[0]*x[0] + y[0]*y[0] ) ;
		Traits::tp( fb ) = Traits::tp( x ) ;
	} else {
		Vector xh = x, yh = y ;
		Traits::np( xh ) *= mu ;
		Traits::tp( yh ) *= mu ;
		compute< false >( xh, yh, fb, a, b ) ;
	}
}

template< unsigned Dimension, typename Scalar >
void FBBaseFunction< Dimension, Scalar >::computeJacobian(
		const Scalar mu, const Vector& x, const Vector& y,
		Vector& fb, Matrix& dFb_dx, Matrix& dFb_dy )
{
	if( NumTraits< Scalar >::isZero( mu ) )
    {
        const Scalar z = std::sqrt( x[0]*x[0] + y[0]*y[0] ) ;
        Traits::np( fb ) = x[0] + y[0] - z ;
		Traits::tp( fb ) = Traits::tp( x ) ;

        dFb_dx.setIdentity() ;
        dFb_dy.setZero() ;

        if( !NumTraits< Scalar >::isZero( z ) )
        {
            dFb_dx(0,0) = 1. - x[0] / z ;
            dFb_dy(0,0) = 1. - y[0] / z ;
        }
      } else {
		Vector xh = x, yh = y ;
		Traits::np( xh ) *= mu ;
		Traits::tp( yh ) *= mu ;
		compute< false >( xh, yh, fb, dFb_dx, dFb_dy ) ;
                dFb_dx.col( 0 ) *= mu ;
                dFb_dy.template block< Dimension, Dimension-1 > ( 0, 1 ) *= mu ;
	}

}

template< unsigned Dimension, typename Scalar >
template< bool JacobianAsWell >
void FBBaseFunction< Dimension, Scalar >::compute(
			const Vector& x, const Vector& y, Vector& fb,
			Matrix& dFb_dx, Matrix& dFb_dy )
{
	(void ) x ;
	(void ) y ;
	(void ) fb ;
	(void ) dFb_dx ;
	(void ) dFb_dy ;
}

template< unsigned Dimension, typename Scalar, bool DeSaxceCOV >
void FischerBurmeister< Dimension, Scalar, DeSaxceCOV >::
compute( const Scalar mu, const Vector& x, const Vector& y, Vector& fb ) 
{
  if ( DeSaxceCOV )
  {
    Vector yt ( y );
    Traits::np( yt ) += mu * Traits::tp( y ).norm() ;
    BaseFunction::compute( mu, x, yt, fb ) ;
  } else {
    BaseFunction::compute( mu, x, y, fb ) ;
  }
}

template< unsigned Dimension, typename Scalar, bool DeSaxceCOV >
void FischerBurmeister< Dimension, Scalar, DeSaxceCOV >::compute(
    const Vector& x, Vector& fb ) const 
{
  const Vector y = m_A * x + m_b ;
  compute( m_mu, x, y, fb ) ;
}
	
template< unsigned Dimension, typename Scalar, bool DeSaxceCOV >
void FischerBurmeister< Dimension, Scalar, DeSaxceCOV >::computeJacobian(
      const Vector& x, Vector& fb, Matrix& dFb_dx ) const 
{
  Vector y ( m_A * x + m_b ) ;
  Scalar s = 0. ;

  if ( DeSaxceCOV )
  {
    s = Traits::tp( y ).norm() ;
    Traits::np( y ) += m_mu * s ;
  }

  Matrix dFb_dy ;
  BaseFunction::computeJacobian( m_mu, x, y, fb, dFb_dx, dFb_dy ) ;
  
  if ( DeSaxceCOV && !NumTraits< Scalar >::isZero( s ) )
  {
    dFb_dy.template block< Dimension, Dimension - 1 >( 0, 1 ).noalias() +=
      dFb_dy.col( 0 ) *  ( m_mu / s ) * Traits::tp( y ).transpose() ;

  }

  dFb_dx.noalias() += dFb_dy * m_A ;

}


}

#endif
