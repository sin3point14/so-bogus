/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2016 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_PRODUCT_GAUSS_SEIDEL_IMPL_HPP
#define BOGUS_PRODUCT_GAUSS_SEIDEL_IMPL_HPP

#include "ProductGaussSeidel.hpp"
#include "GaussSeidelBase.impl.hpp"

#ifndef BOGUS_DONT_PARALLELIZE
#include <omp.h>
#endif

namespace bogus
{

template < typename BlockMatrixType, typename DiagonalType >
ProductGaussSeidel< BlockMatrixType, DiagonalType >&
ProductGaussSeidel< BlockMatrixType, DiagonalType >::setMatrix( const BlockObjectBase< BlockMatrixType > & M )
{
	m_matrix = &M ;

	updateLocalMatrices() ;

	return *this ;
}

namespace block_solvers_impl {

template <typename T>
struct SelfProductAccumulator {

	T& acc ;
	SelfProductAccumulator( T& mat ) : acc(mat)
	{}

	template <typename Index, typename Matrix >
	void operator() (const Index, const Matrix &block )
	{
		acc += block * block.transpose() ;
	}
};

} //block_solvers_impl

template < typename BlockMatrixType, typename DiagonalType >
void ProductGaussSeidel< BlockMatrixType, DiagonalType >::updateLocalMatrices( )
{

	if( !m_matrix )
		return ;

	const Index n = m_matrix->rowsOfBlocks() ;
	m_localMatrices.resize( n ) ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( Index i = 0 ; i <  n ; ++i )
	{
		m_localMatrices[i].setZero() ;

		block_solvers_impl::SelfProductAccumulator< typename Base::DiagonalMatrixType >
				acc( m_localMatrices[i] ) ;
		m_matrix->derived().eachBlockOfRow( i, acc ) ;
	}

	Base::processLocalMatrices() ;
}

template < typename BlockMatrixType, typename DiagonalType >
template < typename NSLaw,  typename VecT, typename ResT >
void ProductGaussSeidel< BlockMatrixType, DiagonalType >::innerLoop(
		bool parallelize, const NSLaw &law, const VecT& b,
		std::vector< unsigned char > &skip, Scalar &ndxRef,
		VecT& Mx, ResT &x	) const
{
	typedef typename NSLaw::Traits LocalProblemTraits ;

	Segmenter< NSLaw::dimension, ResT, typename BlockMatrixType::Index >
			xSegmenter( x, m_matrix->rowOffsets() ) ;
	const Segmenter< NSLaw::dimension, const VecT, typename BlockMatrixType::Index >
			bSegmenter( b, m_matrix->rowOffsets() ) ;

	const Scalar absSkipTol = std::min( m_skipTol, m_tol ) ;
	const Scalar absSkipIters = std::min( m_skipIters, (unsigned) std::sqrt( (Scalar) skip.size() ) ) ;

	const std::ptrdiff_t n = m_matrix->rowsOfBlocks() ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel if ( parallelize )
	{
#endif
		typename LocalProblemTraits::Vector lb, lx, ldx ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp for
#endif
		for( std::ptrdiff_t i = 0 ; i < n ; ++ i )
		{

			if( skip[i] ) {
				--skip[i] ;
				continue ;
			}

			lx = xSegmenter[ i ] ;
			lb = bSegmenter[ i ] - m_localMatrices[i] * lx ;

			m_matrix->derived().template rowMultiplyPrecompose<false>( i, Mx, lb, m_diagonal.get() ) ;

			ldx = -lx ;

			const bool ok = law.solveLocal( i, m_localMatrices[i], lb, lx, m_scaling[ i ] ) ;
			ldx += lx ;

			if( !ok ) { ldx *= .5 ; }
			xSegmenter[ i ] += ldx ;

			m_matrix->derived().template colMultiply<true >( i, ldx, Mx ) ;

			const Scalar nx2 = m_scaling[ i ] * m_scaling[ i ] * lx.squaredNorm() ;
			const Scalar ndx2 = m_scaling[ i ] * m_scaling[ i ] * ldx.squaredNorm() ;
			// Thread-safety:
			// we do not care if we lose an update of ndxRef to a data race,
			// but we need its bytes to be consistent.
			// Ok on x86(_64) if ndxRef is aligned as assignment is atomic
			if( ndx2 > ndxRef ) ndxRef = ndx2 ;

			if(  std::min(nx2, ndx2) < absSkipTol ||
				 ndx2 < m_skipTol * std::min( nx2, ndxRef ) )
			{
				skip[i] = absSkipIters ;
			}
		}

#ifndef BOGUS_DONT_PARALLELIZE
		}
#endif

}



template < typename BlockMatrixType, typename DiagonalType >
template < typename NSLaw, typename RhsT, typename ResT >
typename ProductGaussSeidel< BlockMatrixType, DiagonalType >::Scalar
ProductGaussSeidel< BlockMatrixType, DiagonalType >::solve( const NSLaw &law,
									   const RhsT &b, ResT &x, bool tryZeroAsWell ) const
{
	const bogus::Zero< Scalar > zero( x.rows(), x.rows() ) ;
	return solveWithLinearConstraints( law, zero, b, x, tryZeroAsWell, 0 ) ;
}

template < typename BlockMatrixType, typename DiagonalType >
template < typename NSLaw, typename RhsT, typename ResT, typename LSDerived, typename HDerived >
typename ProductGaussSeidel< BlockMatrixType, DiagonalType >::Scalar
ProductGaussSeidel< BlockMatrixType, DiagonalType >::solveWithLinearConstraints( const NSLaw &law,
								   const BlockObjectBase< LSDerived >& Cinv,
								   const BlockObjectBase<  HDerived >& H,
								   const Scalar alpha,
								   const RhsT &b, const RhsT &c,
								   ResT &x,
								   bool tryZeroAsWell, unsigned solveEvery ) const
{
	const typename GlobalProblemTraits::DynVector k = b + H * Cinv * c ;

	return solveWithLinearConstraints( law, -alpha * H * Cinv * H.transpose(), k, x, tryZeroAsWell, solveEvery ) ;
}


template < typename BlockMatrixType, typename DiagonalType >
template < typename NSLaw, typename RhsT, typename ResT, typename WDerived >
typename ProductGaussSeidel< BlockMatrixType, DiagonalType >::Scalar
ProductGaussSeidel< BlockMatrixType, DiagonalType >::solveWithLinearConstraints( const NSLaw &law,
								   const BlockObjectBase < WDerived >& W,
								   const RhsT &b, ResT &x,
								   bool tryZeroAsWell, unsigned solveEvery) const
{
	assert( m_matrix ) ;
	assert( m_diagonal.valid() ) ;
	assert( solveEvery == 0 || 0 == m_evalEvery % solveEvery ) ;

	typename GlobalProblemTraits::DynVector Mx( m_matrix->cols() ), y, x_best ;

	typename GlobalProblemTraits::DynVector w = b;
	W.template multiply< false >(x, w, 1, 1) ;

	Scalar err_best = std::numeric_limits< Scalar >::max() ;

	y = w ;
	m_matrix->template multiply< true  >( x, Mx, 1, 0 ) ;
	m_matrix->template multiply< false >( m_diagonal.get()*Mx, y, 1, 1 ) ;
	Base::evalAndKeepBest( law, x, y, x_best, err_best ) ;

	if( tryZeroAsWell && Base::tryZero( law, b, x, x_best, err_best ) ) {
		w = b ;
		Mx.setZero() ;
	}

	this->m_callback.trigger( 0, err_best ) ;

	const Index n = m_matrix->rowsOfBlocks() ;

	WithMaxThreads wmt( m_maxThreads ) ;
	const int newMaxThreads = wmt.nThreads() ;
	const bool parallelize = (newMaxThreads != 1 && n > newMaxThreads*newMaxThreads ) ;

	std::vector< unsigned char > skip( n, 0 ) ;
	Scalar ndxRef = 0 ; //Reference step size

	unsigned GSIter ;
	for( GSIter = 1 ; GSIter <= m_maxIters ; ++GSIter )
	{

		innerLoop( parallelize, law, w, skip, ndxRef, Mx, x ) ;

		if( solveEvery > 0 && 0 == ( GSIter % solveEvery ) )
		{
			w = b ;
			W.template multiply< false >(x, w, 1, 1) ;
		}

		if( 0 == ( GSIter % m_evalEvery ) )
		{
			y = w ;
			m_matrix->template multiply< true  >( x, Mx, 1, 0 ) ;
			m_matrix->template multiply< false >( m_diagonal.get()*Mx, y, 1, 1 ) ;
			const Scalar err = Base::evalAndKeepBest( law, x, y, x_best, err_best ) ;

			this->m_callback.trigger( GSIter, err ) ;

			if( err < m_tol )
			{
				break ;
			}

			ndxRef /= m_evalEvery ;
		}

	}

	if( GSIter > m_maxIters ) x = x_best ;

	return err_best ;

}


} //namespace bogus


#endif
