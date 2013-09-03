/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_BLOCK_GAUSS_SEIDEL_IMPL_HPP
#define BOGUS_BLOCK_GAUSS_SEIDEL_IMPL_HPP

#include "GaussSeidel.hpp"
#include "BlockSolverBase.impl.hpp"

#include "../Block/BlockMatrixBase.hpp"

namespace bogus
{

template < typename BlockMatrixType >
void GaussSeidel< BlockMatrixType >::init( )
{
	m_tol = 1.e-6 ;
	m_maxIters = 250 ;
	m_enableColoring =  true ;
	m_useInfinityNorm = false ;
	m_evalEvery = 25  ;
	m_skipTol = m_tol * m_tol ;
	m_skipIters = 10 ;
	m_autoRegularization = 0. ;

}

template < typename BlockMatrixType >
void GaussSeidel< BlockMatrixType >::setMatrix( const BlockMatrixBase< BlockMatrixType > & M )
{
	m_matrix = &M ;

	const unsigned n = M.rowsOfBlocks() ;
	m_localMatrices.resize( n ) ;
	m_scaling.resize( n ) ;
	m_regularization.resize( n ) ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( std::ptrdiff_t i = 0 ; i < (std::ptrdiff_t) n ; ++i )
	{
		m_localMatrices[i] = M.diagonal( i ) ;

		if( m_autoRegularization > 0. )
		{
			m_regularization(i) = std::max( 0., m_autoRegularization - m_localMatrices[i].eigenvalues().real().minCoeff() ) ;
			m_localMatrices[i].diagonal().array() += m_regularization(i) ;
		} else m_regularization(i) = 0. ;

		m_scaling[i] = std::max( 1., m_localMatrices[i].trace() ) ;
	}

	updateColors() ;

}

template < typename BlockMatrixType >
void GaussSeidel< BlockMatrixType >::resetColors( )
{
	assert( m_matrix ) ;

	const std::ptrdiff_t n = static_cast< std::ptrdiff_t >( m_localMatrices.size() ) ;

	m_permutation.resize( n ) ;
	m_colors.clear() ;
	m_colors.push_back( 0 ) ;
	m_colors.push_back( n ) ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( std::ptrdiff_t i = 0 ; i < n ; ++ i ) { m_permutation[i] = i ; }
	m_colorsUpToDate = true ;
}

template < typename BlockMatrixType >
void GaussSeidel< BlockMatrixType >::computeColors()
{

	resetColors() ;
}

template < typename BlockMatrixType >
template < typename NSLaw, typename RhsT, typename ResT >
typename GaussSeidel< BlockMatrixType >::Scalar GaussSeidel< BlockMatrixType >::eval( const NSLaw &law,
							const RhsT &x, const ResT &y ) const
{
	const Segmenter< NSLaw::dimension, const RhsT, typename BlockMatrixType::Index >
			xSegmenter( x, m_matrix->rowOffsets() ) ;
	const Segmenter< NSLaw::dimension, const ResT, typename BlockMatrixType::Index >
			ySegmenter( y, m_matrix->rowOffsets() ) ;

	const std::ptrdiff_t n = (std::ptrdiff_t) m_localMatrices.size() ;

	Scalar err = 0., lres ;
	typename NSLaw::Traits::Vector lx, ly ;

	if( m_useInfinityNorm )
	{

	#ifndef BOGUS_DONT_PARALLELIZE
	#pragma omp parallel private( lx, ly, lres )
	#endif
		{
			lres = 0. ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp for
#endif
			for( int i = 0 ; i < n ; ++ i )
			{
				lx = xSegmenter[ i ] * m_scaling[i] ;
				ly = ySegmenter[ i ] ;
				lres = std::max( law.eval( i, lx, ly ), lres ) ;
			}

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp critical
#endif
			err = std::max( err, lres ) ;
		}

		return err ;

	} else {

	#ifndef BOGUS_DONT_PARALLELIZE
	#pragma omp parallel for private( lx, ly, lres ) reduction ( + : err )
	#endif
		for( int i = 0 ; i < n ; ++ i )
		{
			lx = xSegmenter[ i ] * m_scaling[i] ;
			ly = ySegmenter[ i ] ;
			lres = law.eval( i, lx, ly ) ;
			err += lres ;
		}

		return err / ( 1 + n );

	}
}

template < typename BlockMatrixType >
template < typename NSLaw, typename RhsT, typename ResT >
typename GaussSeidel< BlockMatrixType >::Scalar GaussSeidel< BlockMatrixType >::solve( const NSLaw &law,
							const RhsT &b, ResT &x ) const
{
	assert( m_matrix ) ;

	const Segmenter< NSLaw::dimension, const RhsT, typename BlockMatrixType::Index >
			bSegmenter( b, m_matrix->rowOffsets() ) ;
	Segmenter< NSLaw::dimension, ResT, typename BlockMatrixType::Index >
			xSegmenter( x, m_matrix->rowOffsets() ) ;

	typedef typename NSLaw::Traits LocalProblemTraits ;

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

	const std::ptrdiff_t n = static_cast< std::ptrdiff_t >( m_localMatrices.size() ) ;
	std::vector< unsigned char > skip( n, 0 ) ;

	// see [Daviet et al 2011], Algorithm 1
	for( unsigned GSIter = 1 ; GSIter <= m_maxIters ; ++GSIter )
	{
		typename LocalProblemTraits::Vector lb, lx, ldx ;

		for( unsigned c = 0 ; c+1 < m_colors.size() ; ++ c )
		{

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for private( lb, lx, ldx )
#endif
			for( std::ptrdiff_t pi = m_colors[c] ; pi < m_colors[c+1] ; ++ pi )
			{
				const std::size_t i = m_permutation[pi] ;

				if( skip[i] ) {
					--skip[i] ;
					continue ;
				}

				lx = xSegmenter[ i ] ;
				lb = bSegmenter[ i ] - m_regularization(i) * lx ;
				m_matrix->splitRowMultiply( i, x, lb ) ;
				ldx = -lx ;

				const bool ok = law.solveLocal( i, m_localMatrices[i], lb, lx, m_scaling[ i ] ) ;
				ldx += lx ;

				if( !ok ) ldx *= .7 ;
				xSegmenter[ i ] += ldx ;

				const Scalar scaledSkipTol = m_scaling[ i ] * m_scaling[ i ] * m_skipTol ;
				if( ldx.squaredNorm() < scaledSkipTol || lx.squaredNorm() < scaledSkipTol )
				{
					skip[i] = m_skipIters ;
				}
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
