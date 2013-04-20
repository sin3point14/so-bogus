/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver          
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>                       
 *
 * This Source Code Form is subject to the terms of the Mozilla Public 
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_FISCHER_BURMEISTER_IMPL_HPP
#define BOGUS_FISCHER_BURMEISTER_IMPL_HPP

#include "FischerBurmeister.hpp"
#include <iostream>
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
		compute< true >( xh, yh, fb, dFb_dx, dFb_dy ) ;
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
	// see [Daviet et al 2011], Appendix A.1

	Vector z2 ;

	z2[0] = x.squaredNorm() + y.squaredNorm() ;
	Traits::tp( z2 ) = x[0] * Traits::tp( x ) + y[0] * Traits::tp( y ) ;
	const Scalar nz2t = Traits::tp( z2 ).norm() ;

	Vector omega1, omega2 ;
	omega1[0] = omega2[0] = .5 ;

	if( NumTraits< Scalar >::isZero( nz2t ) )
	{
		Traits::tp( omega1 ).setZero() ;
		omega1[1] = -.5 ;
		Traits::tp( omega2 ).setZero() ;
		omega2[1] =  .5 ;
	} else {
		Traits::tp( omega1 ) = -.5* Traits::tp( z2 ) / nz2t ;
		Traits::tp( omega2 ) =  .5* Traits::tp( z2 ) / nz2t ;
	}

	const Scalar rlambda1 = std::sqrt( std::max( (Scalar) 0, z2[0] - 2*nz2t ) ) ;
	const Scalar rlambda2 = std::sqrt( z2[0] + 2*nz2t ) ;

	const Vector z = rlambda1 * omega1 + rlambda2 * omega2 ;
	fb = x + y - z ;

	if( JacobianAsWell )
	{

		if( NumTraits< Scalar >::isZero( rlambda2 ) )
		{
			// x = y = 0
			dFb_dx.setZero() ;
			dFb_dy.setZero() ;
		} else {
			if( NumTraits< Scalar >::isZero( rlambda1 ) )
			{
				const Scalar izn = 1. / ( x[0]*x[0] + y[0]*y[0] ) ;
				dFb_dx = ( 1. - x[0]*izn ) * Matrix::Identity() ;
				dFb_dy = ( 1. - y[0]*izn ) * Matrix::Identity() ;
			} else {
				const Scalar det = rlambda1 * rlambda2 ;

				Matrix L, invLz ;

				invLz(0, 0) = z[0] ;
				invLz.template block< 1, Dimension - 1>( 0, 1 ) = -Traits::tp( z ).transpose() ;
				invLz.template block< Dimension - 1, 1>( 1, 0 ) = -Traits::tp( z ) ;
				invLz.template block< Dimension - 1, Dimension - 1>( 1, 1 ) =
					( det * MatrixTraits< Dimension - 1, Scalar >::Matrix::Identity() +
					  Traits::tp( z ) * Traits::tp( z ).transpose() ) / z[0] ;
				invLz /= det ;

				L = x[0] * Matrix::Identity() ;
				L.template block< 1, Dimension - 1>( 0, 1 ) = Traits::tp( x ).transpose() ;
				L.template block< Dimension - 1, 1>( 1, 0 ) = Traits::tp( x ) ;
				dFb_dx.setIdentity() ;
				dFb_dx.noalias() -= invLz * L ;

				L = y[0] * Matrix::Identity() ;
				L.template block< 1, Dimension - 1>( 0, 1 ) = Traits::tp( y ).transpose() ;
				L.template block< Dimension - 1, 1>( 1, 0 ) = Traits::tp( y ) ;
				dFb_dy.setIdentity() ;
				dFb_dy.noalias() -= invLz * L ;

			}
		}

	}

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
  const Vector y  = m_A * x + m_b ;
  const Vector xs = m_scaling * x ;
  compute( m_mu, xs, y, fb ) ;
}

template< unsigned Dimension, typename Scalar, bool DeSaxceCOV >
void FischerBurmeister< Dimension, Scalar, DeSaxceCOV >::computeJacobian(
	  const Vector& x, Vector& fb, Matrix& dFb_dx ) const
{
  const Vector xs = m_scaling * x ;
  Vector y ( m_A * x + m_b ) ;
  Scalar s = 0. ;
//  std::cout << m_A << std::endl ;
//  std::cout << m_b.transpose() << std::endl ;

  if ( DeSaxceCOV )
  {
	s = Traits::tp( y ).norm() ;
	Traits::np( y ) += m_mu * s ;
  }
//  std::cout << y.transpose() << std::endl ;

  Matrix dFb_dy ;
  BaseFunction::computeJacobian( m_mu, xs, y, fb, dFb_dx, dFb_dy ) ;
//  std::cout << dFb_dx.transpose() << std::endl ;
//  std::cout << dFb_dy.transpose() << std::endl ;

  if ( DeSaxceCOV && !NumTraits< Scalar >::isZero( s ) )
  {
	dFb_dy.template block< Dimension, Dimension - 1 >( 0, 1 ).noalias() +=
	  dFb_dy.col( 0 ) *  ( m_mu / s ) * Traits::tp( y ).transpose() ;

  }

  dFb_dx *= m_scaling ;
  dFb_dx.noalias() += dFb_dy * m_A ;

}


}

#endif
