/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef BOGUS_CONJUGATE_GRADIENT_INPL_HPP
#define BOGUS_CONJUGATE_GRADIENT_INPL_HPP

#include "ConjugateGradient.hpp"
#include "../Utils/NumTraits.hpp"

namespace bogus {

template < typename BlockMatrixType, template< typename BlockMatrixT > class PreconditionerType >
ConjugateGradient< BlockMatrixType, PreconditionerType >::ConjugateGradient(
        const BlockMatrixBase< BlockMatrixType > & matrix )
    : Base( matrix, matrix.rows(), NumTraits< Scalar >::epsilon() ),
      m_preconditioner( matrix )
{
}

template < typename BlockMatrixType, template< typename BlockMatrixT > class PreconditionerType >
template < typename RhsT, typename ResT >
typename ConjugateGradient< BlockMatrixType, PreconditionerType >::Scalar
ConjugateGradient< BlockMatrixType, PreconditionerType >::solve( const RhsT &b, ResT &x ) const
{
	typedef typename GlobalProblemTraits::DynVector Vector ;

	const Scalar scale = 1. / ( 1 + m_matrix.rows() ) ;

	Vector r = b - m_matrix*x ;
	Scalar res = r.squaredNorm() * scale;
	const Scalar res0 = b.squaredNorm() * scale ;

	if( res > res0 ) {
		x.setZero() ;
		r = b;
		res = res0 ;
	}

	if( res < m_tol ) return res ;

	Vector z ;
	m_preconditioner.template apply< false >( r, z ) ;
	Vector p = z;

	Scalar zr0 = r.dot( z ) ;
	Scalar zr1 ;

	Vector Mp( m_matrix.rows() ) ;

	unsigned k ;
	for( k = 0 ; k < m_maxIters ; ++k )
	{
		Mp.setZero() ;
		m_matrix.template multiply< false >( p, Mp ) ;
		const Scalar alpha = zr0 / ( p.dot( Mp ) ) ;
		x += alpha * p ;
		r -= alpha * Mp ;

		res = r.squaredNorm() * scale ;
		if( res < m_tol ) break ;

		m_preconditioner.template apply< false >( r, z ) ;
		zr1 = z.dot( r ) ;

		p = z + ( zr1 / zr0 ) * p ;

		zr0 = zr1 ;
	}

	return res ;
}

template < typename BlockMatrixType, template< typename BlockMatrixT > class PreconditionerType >
template < typename RhsT, typename ResT >
typename ConjugateGradient< BlockMatrixType, PreconditionerType >::Scalar
ConjugateGradient< BlockMatrixType, PreconditionerType >::solve_BiCG( const RhsT &b, ResT &x ) const
{
	typedef typename GlobalProblemTraits::DynVector Vector ;

	const Scalar scale = 1. / ( 1 + m_matrix.rows() ) ;

	Vector r = b - m_matrix*x ;
	Scalar res = r.squaredNorm() * scale;
	const Scalar res0 = b.squaredNorm() * scale ;

	if( res > res0 ) {
		x.setZero() ;
		r = b;
		res = res0 ;
	}

	if( res < m_tol ) return res ;

	Vector b_ = b ;
	Vector x_ = x ;

	Vector r_ = b_ ;
	m_matrix.template multiply< true >( x_, r_ );

	Vector p, p_ ;
	m_preconditioner.template apply< false >( r,  p ) ;
	m_preconditioner.template apply< true >( r_, p_ ) ;

	Vector z = p, z_ = p_ ;
	Scalar zr0 = r_.dot( z ) ;
	Scalar zr1 ;

	Vector Mp( m_matrix.rows() ) ;

	unsigned k ;
	for( k = 0 ; k < m_maxIters ; ++k )
	{
		Mp.setZero() ;
		m_matrix.template multiply< false >( p, Mp ) ;
		const Scalar alpha = zr0 / ( p_.dot( Mp ) ) ;
		x  += alpha * p  ;
		x_ += alpha * p_ ;

		r  -= alpha * Mp ;
		r_ -= alpha * ( m_matrix.transpose() * p_ );

		res = r.squaredNorm() * scale ;
		if( res < m_tol ) break ;

		m_preconditioner.template apply< false >( r, z ) ;
		zr1 = z.dot( r_ ) ;

		const Scalar beta = ( zr1 / zr0 ) ;

		p = z + beta * p ;
		m_preconditioner.template apply< true >( r_, z_ ) ;
		p_ = z_ + beta * p_ ;

		zr0 = zr1 ;
	}


	return res ;
}


template < typename BlockMatrixType, template< typename BlockMatrixT > class PreconditionerType >
template < typename RhsT, typename ResT >
typename ConjugateGradient< BlockMatrixType, PreconditionerType >::Scalar
ConjugateGradient< BlockMatrixType, PreconditionerType >::solve_BiCGSTAB( const RhsT &b, ResT &x ) const
{
	typedef typename GlobalProblemTraits::DynVector Vector ;

	const Scalar scale = 1. / ( 1 + m_matrix.rows() ) ;

	Vector r = b - m_matrix*x ;
	Scalar res = r.squaredNorm() * scale;

/*	const Scalar res0 = b.squaredNorm() * scale ;
	if( res > res0 ) {
		x.setZero() ;
		r = b;
		res = res0 ;
	}*/

	if( res < m_tol ) return res ;

	Vector r0h = r ;

	Scalar rho0 = 1, rho1 ;
	Scalar alpha = 1, w = 1 ;

	Vector nu = Vector::Zero( r.rows() );
	Vector p = nu ;
	Vector s, t ( m_matrix.rows() ) , u, y, z ;

	unsigned k  ;
	for( k = 0 ; k < m_maxIters ; ++k )
	{
		rho1 = r0h.dot( r ) ;
		const Scalar beta = ( rho1 / rho0 ) * ( alpha / w ) ;
		p = r + beta * ( p - w * nu ) ;
		m_preconditioner.template apply< false >( p, y ) ;

		nu.setZero() ;
		m_matrix.template multiply< false >( y, nu ) ;
		alpha = rho1 / r0h.dot( nu ) ;
		s = r - alpha * nu ;
		m_preconditioner.template apply< false >( s, z ) ;

		t.setZero() ;
		m_matrix.template multiply< false >( z, t ) ;
		m_preconditioner.template apply< false >( t, u ) ;

		if ( NumTraits< Scalar >::isZero( u.squaredNorm() ) )
		{
			w = 1 ;
		} else {
			w = z.dot( u ) / u.squaredNorm() ;
		}

		x += alpha*y + w*z ;
		r = s - w*t ;

		res = r.squaredNorm() * scale;
		if( res < m_tol ) break ;

	}

	return res ;
}

} // namespace bogus

#endif
