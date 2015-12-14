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
#include "Coloring.impl.hpp"
#include "GaussSeidelBase.impl.hpp"

#ifndef BOGUS_DONT_PARALLELIZE
#include <omp.h>
#endif

namespace bogus
{

template < typename BlockMatrixType >
GaussSeidel< BlockMatrixType >& GaussSeidel< BlockMatrixType >::setMatrix( const BlockObjectBase< BlockMatrixType > & M )
{
	if( m_matrix != &M && ( m_matrix != BOGUS_NULL_PTR( const BlockObjectBase< BlockMatrixType >) ||
							m_coloring.size() != (std::size_t) M.rowsOfBlocks() )) {
		m_coloring.update( false, M.derived() );
	}

	m_matrix = &M ;

	Base::updateLocalMatrices() ;

	return *this ;
}

template < typename BlockMatrixType >
template < typename NSLaw, typename RhsT, typename ResT >
typename GaussSeidel< BlockMatrixType >::Scalar GaussSeidel< BlockMatrixType >::solve( const NSLaw &law,
																						 const RhsT &b, ResT &x, bool tryZeroAsWell ) const
{
	assert( m_matrix ) ;

	const Segmenter< NSLaw::dimension, const RhsT, typename BlockMatrixType::Index >
			bSegmenter( b, m_matrix->rowOffsets() ) ;
	Segmenter< NSLaw::dimension, ResT, typename BlockMatrixType::Index >
			xSegmenter( x, m_matrix->rowOffsets() ) ;

	typedef typename NSLaw::Traits LocalProblemTraits ;

	typename GlobalProblemTraits::DynVector y ( b + (*m_matrix)*x ) ;
	typename GlobalProblemTraits::DynVector x_best( GlobalProblemTraits::DynVector::Zero( x.rows() ) ) ;

	const Scalar err_init = Base::eval( law, y, x      ) ;
	Scalar err_zero, err_best = err_init ;

	bool useZero = false ;

	if( tryZeroAsWell )
	{
		err_zero = Base::eval( law, b, x_best ) ;

		if( err_zero < err_init )
		{
			err_best = err_zero ;
			x.setZero() ;
			useZero = true ;
		}
	}

	if( !useZero ) {
		x_best = x ;
	}
	this->m_callback.trigger( 0, err_best ) ;

	const Index n = m_matrix->rowsOfBlocks() ;
	std::vector< unsigned char > skip( n, 0 ) ;

#ifndef BOGUS_DONT_PARALLELIZE
	const int curMaxThreads = omp_get_max_threads() ;
	const int newMaxThreads = m_maxThreads == 0 ? curMaxThreads : m_maxThreads ;
	omp_set_num_threads( newMaxThreads ) ;
#endif

	Scalar ndxRef = 0 ; //Reference step size
	const Scalar absSkipTol = std::min( m_skipTol, m_tol ) ;
	const Scalar absSkipIters = std::min( m_skipIters, (unsigned) std::sqrt( (Scalar) n ) ) ;

	unsigned GSIter ;
	for( GSIter = 1 ; GSIter <= m_maxIters ; ++GSIter )
	{

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel if (m_maxThreads != 1 && n > newMaxThreads*newMaxThreads )
		{
#endif

			typename LocalProblemTraits::Vector lb, lx, ldx ;
			for( unsigned c = 0 ; c+1 < m_coloring.colors.size() ; ++ c )
			{

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp for
#endif
				for( std::ptrdiff_t pi = m_coloring.colors[c] ; pi < m_coloring.colors[c+1] ; ++ pi )
				{

					const std::size_t i = m_coloring.permutation[ pi ] ;

					if( skip[i] ) {
						--skip[i] ;
						continue ;
					}

					lx = xSegmenter[ i ] ;
					lb = bSegmenter[ i ] - m_regularization(i) * lx ;
					Base::explicitMatrix().splitRowMultiply( i, x, lb ) ;
					ldx = -lx ;

					const bool ok = law.solveLocal( i, m_localMatrices[i], lb, lx, m_scaling[ i ] ) ;
					ldx += lx ;

					if( !ok ) { ldx *= .5 ; }
					xSegmenter[ i ] += ldx ;

					const Scalar nx2 = m_scaling[ i ] * m_scaling[ i ] * lx.squaredNorm() ;
					const Scalar ndx2 = m_scaling[ i ] * m_scaling[ i ] * ldx.squaredNorm() ;
					if( ndx2 > ndxRef ) ndxRef = ndx2 ;

					if(  std::min(nx2, ndx2) < absSkipTol ||
						 ndx2 < m_skipTol * std::min( nx2, ndxRef ) )
					{
						skip[i] = absSkipIters ;
					}
				}

			}
#ifndef BOGUS_DONT_PARALLELIZE
		}
#endif

		if( 0 == ( GSIter % m_evalEvery ) )
		{

			y = b ;
			m_matrix->template multiply< false >( x, y, 1, 1 ) ;
			const double err = Base::eval( law, y, x ) ;

			this->m_callback.trigger( GSIter, err ) ;

			if( err < m_tol )
			{
				err_best = err ;
				break ;
			}

			if( err < err_best )
			{
				x_best = x ;
				err_best = err ;
			}

			ndxRef /= m_evalEvery ;
		}

	}

#ifndef BOGUS_DONT_PARALLELIZE
	omp_set_num_threads( curMaxThreads ) ;
#endif

	if( GSIter > m_maxIters ) x = x_best ;

	return err_best ;

}


} //namespace bogus


#endif
