/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_CONSTRAINED_SOLVER_BASE_IMPL_HPP
#define BOGUS_CONSTRAINED_SOLVER_BASE_IMPL_HPP

#include "ConstrainedSolverBase.hpp"

#include "../Block/Access.hpp"
#include "../Block/BlockMatrixBase.hpp"

namespace bogus {


template < typename Derived, typename BlockMatrixType >
template < typename NSLaw, typename RhsT, typename ResT >
typename ConstrainedSolverBase< Derived, BlockMatrixType >::Scalar
ConstrainedSolverBase< Derived, BlockMatrixType >::eval( const NSLaw &law,
							const ResT &y, const RhsT &x ) const
{
	const Segmenter< NSLaw::dimension, const RhsT, typename BlockMatrixType::Index >
			xSegmenter( x, m_matrix->rowOffsets() ) ;
	const Segmenter< NSLaw::dimension, const ResT, typename BlockMatrixType::Index >
			ySegmenter( y, m_matrix->rowOffsets() ) ;

	typedef typename BlockMatrixTraits< BlockMatrixType >::Index Index ;

	const Index n = m_matrix->rowsOfBlocks() ;

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
			for( Index i = 0 ; i < n ; ++ i )
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
		for( Index i = 0 ; i < n ; ++ i )
		{
			lx = xSegmenter[ i ] * m_scaling[i] ;
			ly = ySegmenter[ i ] ;
			lres = law.eval( i, lx, ly ) ;
			err += lres ;
		}

		return err / ( 1 + n );

	}
}

template< typename Derived >
typename BlockObjectBase< Derived >::Scalar estimate_row_scaling( const BlockObjectBase< Derived >& ,
					 typename BlockObjectBase< Derived >::Index )
{
	return 1. ;
}

template< typename Derived >
typename Derived::BlockMatrixType::Scalar estimate_row_scaling( const typename Derived::BlockMatrixType& mat,
					 typename Derived::BlockMatrixType::Index row )
{
	typedef typename BlockMatrixTraits< Derived >::BlockType LocalMatrixType ;
	typedef ProblemTraits< LocalMatrixType > GlobalProblemTraits ;

	if ( mat.rows() != mat.cols() || mat.rowsOfBlocks() != mat.colsOfBlocks() ) {
		// Non-square matrix
		return 1. ;
	}

	const typename Derived::BlockPtr ptr = mat.diagonalBlockPtr( row ) ;
	if( ptr == Derived::InvalidBlockPtr ) {
		// Empty diagonal block
		return 1. ;
	} else {
		return std::max( 1., GlobalProblemTraits::asConstMatrix( mat.block( ptr ) ).trace() ) ;
	}
}

template < typename Derived, typename BlockMatrixType >
void ConstrainedSolverBase< Derived,BlockMatrixType >::updateScalings()
{
	if( !m_matrix )
	{
		return ;
	}

	const Index n = m_matrix->rowsOfBlocks() ;
	m_scaling.resize( n ) ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( Index i = 0 ; i <  n ; ++i )
	{
		m_scaling[i] = estimate_row_scaling( m_matrix->derived(), i ) ;
	}

}

template < typename Derived, typename BlockMatrixType >
template < typename NSLaw, typename VectorT >
void ConstrainedSolverBase<  Derived,BlockMatrixType >::projectOnConstraints(
		const NSLaw &law, VectorT &x ) const
{
	Segmenter< NSLaw::dimension, VectorT, typename BlockMatrixType::Index >
			xSegmenter( x, m_matrix->rowOffsets() ) ;

	const Index n = m_matrix->rowsOfBlocks() ;
	typename NSLaw::Traits::Vector lx ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for private( lx )
#endif
	for( Index i = 0 ; i < n ; ++ i )
	{
		lx = xSegmenter[ i ] ;
		law.projectOnConstraint( i, lx ) ;
		xSegmenter[ i ] = lx ;
	}

}

template < typename Derived, typename BlockMatrixType >
template < typename NSLaw, typename RhsT, typename ResT >
void ConstrainedSolverBase<  Derived,BlockMatrixType >::dualityCOV(
		const NSLaw &law, const RhsT &u, ResT &s ) const
{
	const Segmenter< NSLaw::dimension, const RhsT, typename BlockMatrixType::Index >
			uSegmenter( u, m_matrix->rowOffsets() ) ;
	Segmenter< NSLaw::dimension, ResT, typename BlockMatrixType::Index >
			sSegmenter( s, m_matrix->rowOffsets() ) ;

	const Index n = m_matrix->rowsOfBlocks() ;
	typename NSLaw::Traits::Vector ls ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for private( ls )
#endif
	for( Index i = 0 ; i < n ; ++ i )
	{
		law.dualityCOV( i, uSegmenter[i], ls ) ;
		sSegmenter[ i ] = ls ;
	}

}


} // namespace bogus

#endif
