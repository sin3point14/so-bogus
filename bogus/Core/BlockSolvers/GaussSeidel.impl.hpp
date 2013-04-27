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

#include <iostream>

namespace bogus
{


template < typename BlockMatrixType >
GaussSeidel< BlockMatrixType >::GaussSeidel( const BlockMatrixBase< BlockMatrixType > & M )
	: Base( M, 250, 1.e-6 ), m_deterministic( true ),
	  m_evalEvery( 25 ), m_skipTol( m_tol * m_tol ), m_skipIters( 10 )
{
	const unsigned d = GlobalProblemTraits::dimension ;
	const unsigned n = M.rowsOfBlocks() ;
	m_localMatrices.resize( n ) ;
	m_scaling.resize( n*d ) ;

	for( unsigned i = 0 ; i < n ; ++i )
	{
		m_localMatrices[i] = M.diagonal( i ) ;
		GlobalProblemTraits::segment( i, m_scaling )
				.setConstant( std::max( 1., m_localMatrices[i].trace() ) );
	}

}

template < typename BlockMatrixType >
template < typename NSLaw, typename RhsT, typename ResT >
typename GaussSeidel< BlockMatrixType >::Scalar GaussSeidel< BlockMatrixType >::solve( const NSLaw &law,
							const RhsT &b, ResT &x ) const
{
	typedef LocalProblemTraits< GlobalProblemTraits::dimension, typename GlobalProblemTraits::Scalar > LocalProblemTraits ;

	const unsigned d = GlobalProblemTraits::dimension ;
	const unsigned n = m_localMatrices.size() ;
	assert( n*d == b.rows() ) ;
	assert( n*d == x.rows() ) ;

	typename GlobalProblemTraits::DynVector y ( b + m_matrix*x ) ;
	typename GlobalProblemTraits::DynVector x_best( GlobalProblemTraits::DynVector::Zero( x.rows() ) ) ;
	typename GlobalProblemTraits::DynVector x_scaled( x.array() * m_scaling.array() ) ;

	const Scalar err_init = law.eval( x_scaled, y ) ;
	const Scalar err_zero = law.eval( x_best, b ) ;

	Scalar err_best ;
	if( err_zero < err_init )
	{
		err_best = err_zero ;
		x.setZero() ;
	} else {
		err_best = err_init ;
		x_best = x ;
	}

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
			lb = GlobalProblemTraits::segment( i, b ) ;
			m_matrix.splitRowMultiply( i, x, lb ) ;
			lx = GlobalProblemTraits::segment( i, x ) ;
			ldx = -lx ;

			const bool ok = law.solveLocal( i, m_localMatrices[i], lb, lx, m_scaling[ d*i ] ) ;
			ldx += lx ;

			if( !ok ) ldx *= .7 ;
			GlobalProblemTraits::segment( i, x ) += ldx ;

			const Scalar scaledSkipTol = m_scaling[ d*i ] * m_scaling[ d*i ] * m_skipTol ;
			if( ldx.squaredNorm() < scaledSkipTol || lx.squaredNorm() < scaledSkipTol )
			{
				skip[i] = m_skipIters ;
			}
		}

		if( 0 == ( GSIter % m_evalEvery ) )
		{
			y = b + m_matrix * x ;
			x_scaled = x.array() * m_scaling.array() ;
			const double err = law.eval( x_scaled, y ) ;

			std::cout << "Finished iteration " << GSIter
					  << " with residual " << err
					  << " (target: " << m_tol << " )" << std::endl ;

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
