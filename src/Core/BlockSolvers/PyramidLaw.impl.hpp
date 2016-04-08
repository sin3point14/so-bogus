/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_PYRAMIDLAW_IMPL_HPP
#define BOGUS_PYRAMIDLAW_IMPL_HPP

#include "PyramidLaw.hpp"

#include "../BlockSolvers.hpp"
#include "../Utils/NumTraits.hpp"

namespace bogus {

template < DenseIndexType Dimension, typename Scalar >
PyramidLaw< Dimension, Scalar >::PyramidLaw(const unsigned n, const double *mu )
	: m_mu(mu), m_n(n)
{
}

template < DenseIndexType Dimension, typename Scalar >
void PyramidLaw< Dimension, Scalar >::projectOnConstraint( const unsigned problemIndex,
																	   typename Traits::Vector &x ) const
{
	typedef typename LocalProblemTraits< Dimension-1, Scalar >::Array TgComp ;
	const Scalar mu   = m_mu[ problemIndex ] ;

	const Scalar mun  = mu * Traits::np(x) ;

	if( mun < Traits::tp(x).template lpNorm< Eigen::Infinity >() ) {

		TgComp xp = Traits::tp(x).array().max( 0 ) ;
		TgComp xm = xp - Traits::tp(x).array() ;

		// Displ vector to pyramid (at constant altitude)
		TgComp depl = ( (xm - xm.min( mun )) - (xp - xp.min( mun )) ) / (1 + mu*mu) ;

		// Projection
		Traits::tp(x) += depl.matrix() ;
		Traits::np(x) += mu * depl.matrix().template lpNorm< Eigen::Infinity >()  ;

		// Inside normal cone
		if( Traits::np(x) < NumTraits< Scalar >::epsilon() ) {
			x.setZero() ;
		}
	}
}

template < DenseIndexType Dimension, typename Scalar >
bool PyramidLaw< Dimension, Scalar >::solveLocal(const unsigned problemIndex,
			const typename Traits::Matrix &A,
			const typename Traits::Vector &b,
			typename Traits::Vector &x , const Scalar ) const
{

	const Scalar mu  = m_mu[ problemIndex ] ;

	typename Traits::Vector xp = x.array().max( 0 ) ;
	typename Traits::Vector xm = xp - x ;

	// TODO solve analytically or with fb, deal with not positive diag coeff
	for( unsigned it = 0; it < 2*Dimension ; ++it ) {

		xm[0] = 0 ;
		Scalar l = A.col(0).dot( xp - xm ) - A(0,0)*xp[0] + b[0] ;
		xp[0] = std::max( 0., -l / A(0,0) ) ;

		Scalar mun = mu * xp[0] ;

		for( unsigned k = 1; k < Dimension ; ++k )
		{
			Scalar l = A.col(k).dot( xp - xm ) - A(k,k)*(xp[k]-xm[k]) + b[k] ;
			xp[k] = std::min( mun, -std::min(0., l) / A(k,k) ) ;
			xm[k] = std::min( mun,  std::max(0., l) / A(k,k) ) ;
		}
	}

	x = xp - xm ;

	return true ;
}


template < DenseIndexType Dimension, typename Scalar >
Scalar PyramidLaw< Dimension, Scalar >::eval( const unsigned problemIndex,
							   const typename Traits::Vector &x,
							   const typename Traits::Vector &y ) const
{
	typedef typename LocalProblemTraits< Dimension-1, Scalar >::Array TgComp ;
	const Scalar mu  = m_mu[ problemIndex ] ;

	TgComp xp = Traits::tp(x).array().max( 0 ) ;
	TgComp yp = Traits::tp(y).array().max( 0 ) ;
	TgComp xm = xp - Traits::tp(x).array() ;
	TgComp ym = yp - Traits::tp(y).array() ;

	const Scalar xn = Traits::np(x) ;
	const Scalar yn = Traits::np(y) ;

	const Scalar fb = xn + yn - std::sqrt( xn*xn + yn*yn ) ;
	Scalar err = fb * fb ;

	xm = mu*Traits::np(x) - xm;
	xp = mu*Traits::np(x) - xp;

	TgComp fbm = xm + yp - ( xm*xm + yp*yp ).sqrt() ;
	TgComp fbp = xp + ym - ( xp*xp + ym*ym ).sqrt() ;

	err += fbm.matrix().squaredNorm() + fbp.matrix().squaredNorm() ;
	return err ;
}

} // bogus


#endif
