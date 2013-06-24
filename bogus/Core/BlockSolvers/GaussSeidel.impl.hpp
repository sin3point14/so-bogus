/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_BLOCK_GAUSS_SEIDEL_IMPL_HPP
#define BOGUS_BLOCK_GAUSS_SEIDEL_IMPL_HPP

#include "GaussSeidel.hpp"
#include "BlockSolverBase.impl.hpp"
#include "../Block/BlockMatrix.hpp"

namespace bogus
{

template < typename BlockMatrixType >
GaussSeidel< BlockMatrixType >::GaussSeidel( )
	: Base( NULL, 250, 1.e-6 ), m_deterministic( true ),
	  m_evalEvery( 25 ), m_skipTol( m_tol * m_tol ), m_skipIters( 10 ),
	  m_autoRegularization( 0. )
{
}

template < typename BlockMatrixType >
GaussSeidel< BlockMatrixType >::GaussSeidel( const BlockMatrixBase< BlockMatrixType > & M )
	: Base( &M, 250, 1.e-6 ), m_deterministic( true ),
	  m_evalEvery( 25 ), m_skipTol( m_tol * m_tol ), m_skipIters( 10 ),
	  m_autoRegularization( 0. )
{
	setMatrix( M ) ;
}

template < typename BlockMatrixType >
void GaussSeidel< BlockMatrixType >::setMatrix( const BlockMatrixBase< BlockMatrixType > & M )
{
	m_matrix = &M ;

	const unsigned n = M.rowsOfBlocks() ;
	m_localMatrices.resize( n ) ;
	m_scaling.resize( n ) ;
	m_regularization.resize( n ) ;

	for( unsigned i = 0 ; i < n ; ++i )
	{
		m_localMatrices[i] = M.diagonal( i ) ;

		if( m_autoRegularization > 0. )
		{
			m_regularization(i) = std::max( 0., m_autoRegularization - m_localMatrices[i].eigenvalues().real().minCoeff() ) ;
			m_localMatrices[i].diagonal().array() += m_regularization(i) ;
		} else m_regularization(i) = 0. ;

		m_scaling[i] = std::max( 1., m_localMatrices[i].trace() ) ;
	}

}

template < typename BlockMatrixType >
template < typename NSLaw, typename RhsT, typename ResT >
typename GaussSeidel< BlockMatrixType >::Scalar GaussSeidel< BlockMatrixType >::eval( const NSLaw &law,
							const RhsT &x, const ResT &y ) const
{
	const int n = (int) m_localMatrices.size() ;
	Scalar sum = 0., res ;

	typename NSLaw::Traits::Vector lx, ly ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for private( lx, ly, res ) reduction ( + : sum )
#endif
	for( int i = 0 ; i < n ; ++ i )
	{
		lx = m_matrix->rowSegment( x, i ) * m_scaling[i] ;
		ly = m_matrix->rowSegment( y, i ) ;
		res = law.eval( i, lx, ly ) ;
		sum += res ;
	}

	return sum / ( 1 + n );

}

template < typename BlockMatrixType >
template < typename NSLaw, typename RhsT, typename ResT >
typename GaussSeidel< BlockMatrixType >::Scalar GaussSeidel< BlockMatrixType >::solve( const NSLaw &law,
							const RhsT &b, ResT &x ) const
{
	assert( m_matrix ) ;

	typedef LocalProblemTraits< GlobalProblemTraits::dimension, typename GlobalProblemTraits::Scalar > LocalProblemTraits ;

	const unsigned n = m_localMatrices.size() ;

	typename GlobalProblemTraits::DynVector y ( b + (*m_matrix)*x ) ;
	typename GlobalProblemTraits::DynVector x_best( GlobalProblemTraits::DynVector::Zero( x.rows() ) ) ;

	const Scalar err_init = eval( law, x, y ) ;
	const Scalar err_zero = eval( law, x_best, b ) ;

	Scalar err_best ;
	if( err_zero < err_init )
	{
		err_best = err_zero ;
		x.setZero() ;
	} else {
		err_best = err_init ;
		x_best = x ;
	}

	this->m_callback.trigger( 0, err_best ) ;
//	std::cout << err_init << " /// " << err_zero << std::endl ;

	std::vector< unsigned > skip( n, 0 ) ;

	// see [Daviet et al 2011], Algorithm 1
	for( unsigned GSIter = 1 ; GSIter <= m_maxIters ; ++GSIter )
	{
		typename LocalProblemTraits::Vector lb, lx, ldx ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for private( lb, lx, ldx ) if( !m_deterministic )
#endif
		for( int i = 0 ; i < (int) n ; ++ i )
		{
			if( skip[i] ) {
				--skip[i] ;
				continue ;
			}
			lb = m_matrix->rowSegment( b, i ) ;
			m_matrix->splitRowMultiply( i, x, lb ) ;
			lx = m_matrix->rowSegment( x, i ) ;
			ldx = -lx ;
			lb -= m_regularization(i) * lx ;

			const bool ok = law.solveLocal( i, m_localMatrices[i], lb, lx, m_scaling[ i ] ) ;
			ldx += lx ;

			if( !ok ) ldx *= .7 ;
			m_matrix->rowSegment( x, i ) += ldx ;

			const Scalar scaledSkipTol = m_scaling[ i ] * m_scaling[ i ] * m_skipTol ;
			if( ldx.squaredNorm() < scaledSkipTol || lx.squaredNorm() < scaledSkipTol )
			{
				skip[i] = m_skipIters ;
			}
		}

		if( 0 == ( GSIter % m_evalEvery ) )
		{
			y = b + (*m_matrix) * x ;
			const double err = eval( law, x, y ) ;

			this->m_callback.trigger( GSIter, err ) ;

			if( err < m_tol )
			{
				return err ;
			}

			if( err < err_best )
			{
				x_best = x ;
				err_best = err ;
			}
		}

	}

	x = x_best ;
	return err_best ;

}


} //namespace bogus


#endif
