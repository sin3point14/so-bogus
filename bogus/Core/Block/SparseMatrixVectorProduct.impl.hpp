/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_SPARSE_MATRIXVECTOR_PRODUCT_HPP
#define BOGUS_SPARSE_MATRIXVECTOR_PRODUCT_HPP

#include "SparseBlockMatrix.hpp"
#include "Expressions.hpp"
#include "Access.hpp"

namespace bogus {



template < bool Transpose, typename BlockT, typename IndexT, typename RhsT, typename ResT, typename ScalarT >
static void innerRowMultiply( const BlockT* blocks, const IndexT &index,
							const typename IndexT::Index outerIdx, const RhsT& rhs, ResT& res, ScalarT alpha )
{
	const Segmenter< BlockDims< BlockT, Transpose >::Cols, const RhsT, typename IndexT::Index >
			segmenter( rhs, index.innerOffsetsData() ) ;

	for( typename IndexT::InnerIterator it( index, outerIdx ) ; it ; ++ it )
	{
		res += alpha * ( BlockGetter< Transpose >::get( blocks[ it.ptr() ] ) * segmenter[ it.inner() ] )  ;
	}
}

template < bool Transpose, typename BlockT, typename IndexT, typename RhsT, typename ResT, typename ScalarT >
static void innerColMultiply( const BlockT* blocks, const IndexT &index,
							  const typename IndexT::Index outerIdx, const RhsT& rhs, ResT& res, ScalarT alpha )
{
	Segmenter< BlockDims< BlockT, Transpose >::Rows, ResT, typename IndexT::Index >
			segmenter( res, index.innerOffsetsData() ) ;

	for( typename IndexT::InnerIterator it( index, outerIdx ) ; it ; ++ it )
	{
		segmenter[ it.inner() ] += alpha * ( BlockGetter< Transpose >::get( blocks[ it.ptr() ] ) * rhs ) ;
	}
}

template < typename Derived >
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::splitRowMultiply( const Index row, const RhsT& rhs, ResT& res ) const
{
	assert( Traits::is_symmetric || !Traits::is_col_major ) ;

	const Segmenter< BlockDims< BlockType, Traits::is_col_major >::Cols, const RhsT, Index >
			segmenter( rhs, this->colMajorIndex().innerOffsetsData() ) ;

	for( typename SparseBlockMatrixBase< Derived >::MajorIndexType::InnerIterator it( m_majorIndex, row ) ;
		 it ; ++ it )
	{
		if( it.inner() != row )
			res += BlockGetter< Traits::is_col_major >::get( block( it.ptr() ) ) * segmenter[ it.inner() ] ;
	}
	if( Traits::is_symmetric )
	{
		if( Traits::transpose_can_be_cached && m_transposeIndex.valid )
		{
			const bool is_col_major = ((bool) Traits::is_col_major) ^ ( !Traits::transpose_can_be_cached )  ;
			innerRowMultiply< is_col_major >( this->data(), m_transposeIndex, row, rhs, res, 1 ) ;
		} else {
			assert( m_minorIndex.valid ) ;
			innerRowMultiply< !Traits::is_col_major >( this->data(), m_minorIndex, row, rhs, res, 1 ) ;
		}
	}
}

template < bool Symmetric, bool NativeOrder, bool Transpose >
struct SparseBlockMatrixVectorMultiplier
{
	template < typename Derived, typename RhsT, typename ResT, typename ScalarT >
	static void multiply( const SparseBlockMatrixBase< Derived >& matrix,  const RhsT& rhs, ResT& res, ScalarT alpha )
	{
		const int ResSegDim = BlockDims< typename Derived::BlockType, Transpose >::Rows ;
		typedef Segmenter< ResSegDim, ResT, typename Derived::Index > ResSegmenter ;
		ResSegmenter resSegmenter( res, matrix.minorIndex().innerOffsetsData() ) ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
		for( int i = 0 ; i < (int) matrix.majorIndex().outerSize() ; ++i )
		{
			typename ResSegmenter::ReturnType seg = resSegmenter[ i ] ;
			innerRowMultiply< Transpose >( matrix.data(), matrix.majorIndex(), i, rhs, seg, alpha ) ;
		}
	}
} ;

template < bool NativeOrder, bool Transpose >
struct SparseBlockMatrixVectorMultiplier< true, NativeOrder, Transpose >
{
#ifndef BOGUS_DONT_PARALLELIZE
	template < typename Derived, typename RhsT, typename ResT, typename LocalResT, typename ScalarT >
	static void multiplyAndReduct( const SparseBlockMatrixBase< Derived >& matrix,  const RhsT& rhs, ResT& res, const LocalResT&, ScalarT alpha )
	{
		const int SegDim = BlockDims< typename Derived::BlockType, Transpose >::Rows ;

		typedef Segmenter< SegDim, LocalResT, typename Derived::Index > ResSegmenter ;
		typedef Segmenter< SegDim, const RhsT, typename Derived::Index > RhsSegmenter ;
		const RhsSegmenter rhsSegmenter( rhs, matrix.majorIndex().innerOffsetsData() ) ;

		typedef typename SparseBlockMatrixBase< Derived >::MajorIndexType MajorIndexType ;
#pragma omp parallel
		{
			LocalResT locRes( res.rows() ) ;
			locRes.setZero() ;

			ResSegmenter resSegmenter( locRes, matrix.minorIndex().innerOffsetsData() ) ;

#pragma omp for
			for( int i = 0 ; i < (int) matrix.majorIndex().outerSize() ; ++i )
			{
				typename RhsSegmenter::ConstReturnType rhs_seg( rhsSegmenter[ i ] ) ;
				typename ResSegmenter::ReturnType res_seg(resSegmenter[ i ] ) ;
				for( typename MajorIndexType::InnerIterator it( matrix.majorIndex(), i ) ;
					 it ; ++ it )
				{
					const typename Derived::BlockType &b = matrix.block( it.ptr() ) ;
					res_seg += alpha * ( b *  rhsSegmenter[ it.inner() ] ) ;
					if( it.inner() != i )
						resSegmenter[ it.inner() ] += alpha * ( transpose_block( b ) * rhs_seg ) ;
				}
			}

#pragma omp critical
			res += locRes ;
		}
	}
#endif

	template < typename Derived, typename RhsT, typename ResT, typename ScalarT >
	static void multiply( const SparseBlockMatrixBase< Derived >& matrix,  const RhsT& rhs, ResT& res, ScalarT alpha )
	{
		const int SegDim = BlockDims< typename Derived::BlockType, Transpose >::Rows ;

		typedef Segmenter< SegDim, ResT, typename Derived::Index > ResSegmenter ;
		ResSegmenter resSegmenter( res, matrix.minorIndex().innerOffsetsData() ) ;
		typedef Segmenter< SegDim, const RhsT, typename Derived::Index > RhsSegmenter ;
		const RhsSegmenter rhsSegmenter( rhs, matrix.majorIndex().innerOffsetsData() ) ;

		typedef typename SparseBlockMatrixBase< Derived >::Index Index ;
		if( matrix.transposeIndex().valid )
		{
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
			for( int i = 0 ; i < (int) matrix.majorIndex().outerSize() ; ++i )
			{
				typename ResSegmenter::ReturnType seg( resSegmenter[ i ] ) ;
				innerRowMultiply< false >( matrix.data(), matrix.majorIndex()    , i, rhs, seg, alpha ) ;
				innerRowMultiply< false >( matrix.data(), matrix.transposeIndex(), i, rhs, seg, alpha ) ;
			}
		} else {
#ifdef BOGUS_DONT_PARALLELIZE
			for( Index i = 0 ; i < matrix.majorIndex().outerSize() ; ++i )
			{
				typename RhsSegmenter::ConstReturnType rhs_seg( rhsSegmenter[ i ] ) ;
				typename ResSegmenter::ReturnType      res_seg( resSegmenter[ i ] ) ;
				for( typename SparseBlockMatrixBase< Derived >::MajorIndexType::InnerIterator it( matrix.majorIndex(), i ) ;
					 it ; ++ it )
				{
					const typename Derived::BlockType &b = matrix.block( it.ptr() ) ;
					res_seg += alpha * ( b *  rhsSegmenter[ it.inner() ] ) ;
					if( it.inner() != i )
						resSegmenter[ it.inner() ] += alpha * ( transpose_block( b ) * rhs_seg ) ;
				}
			}
#else
			multiplyAndReduct( matrix, rhs, res, getBlockProductResVec( res ), alpha ) ;
#endif
		}
	}
} ;

template < bool Transpose >
struct OutOfOrderSparseBlockMatrixVectorMultiplier
{

#ifndef BOGUS_DONT_PARALLELIZE
	template < typename Derived, typename RhsT, typename ResT, typename LocalResT, typename ScalarT >
	static void multiplyAndReduct( const SparseBlockMatrixBase< Derived >& matrix,  const RhsT& rhs, ResT& res, const LocalResT&, ScalarT alpha )
	{
		typedef typename SparseBlockMatrixBase< Derived >::Index Index ;

		const int RhsSegDim = BlockDims< typename Derived::BlockType, Transpose >::Cols ;
		typedef Segmenter< RhsSegDim, const RhsT, Index > RhsSegmenter ;
		const RhsSegmenter rhsSegmenter( rhs, matrix.minorIndex().innerOffsetsData() ) ;

#pragma omp parallel
		{
			LocalResT locRes( res.rows() ) ;
			locRes.setZero() ;

#pragma omp for
			for( int i = 0 ; i < (int) matrix.majorIndex().outerSize() ; ++i )
			{
				innerColMultiply< Transpose >( matrix.data(), matrix.majorIndex(), i, rhsSegmenter[i], locRes, alpha ) ;
			}

#pragma omp critical
			res += locRes ;
		}
	}
#endif

	template < typename Derived, typename RhsT, typename ResT, typename ScalarT >
	static void multiply( const SparseBlockMatrixBase< Derived >& matrix,  const RhsT& rhs, ResT& res, ScalarT alpha )
	{
		typedef typename SparseBlockMatrixBase< Derived >::Index Index ;

		const int RhsSegDim = BlockDims< typename Derived::BlockType, Transpose >::Cols ;
		typedef Segmenter< RhsSegDim, const RhsT, Index > RhsSegmenter ;
		const RhsSegmenter rhsSegmenter( rhs, matrix.minorIndex().innerOffsetsData() ) ;

#ifdef BOGUS_DONT_PARALLELIZE
		for( Index i = 0 ; i < matrix.majorIndex().outerSize() ; ++i )
		{
			innerColMultiply< Transpose >( matrix.data(), matrix.majorIndex(), i, rhsSegmenter[i], res, alpha ) ;
		}
#else
		multiplyAndReduct( matrix, rhs, res, getBlockProductResVec( res ), alpha ) ;
#endif
	}

} ;

template < bool Transpose >
struct SparseBlockMatrixVectorMultiplier< false, false, Transpose >
		: public OutOfOrderSparseBlockMatrixVectorMultiplier< Transpose >
{} ;

template < >
struct SparseBlockMatrixVectorMultiplier< false, false, true >
		: public OutOfOrderSparseBlockMatrixVectorMultiplier< true >
{
	typedef OutOfOrderSparseBlockMatrixVectorMultiplier< true > Base;

	template < typename Derived, typename RhsT, typename ResT, typename ScalarT >
	static void multiply( const SparseBlockMatrixBase< Derived >& matrix,  const RhsT& rhs, ResT& res, ScalarT alpha )
	{
		typedef BlockMatrixTraits< Derived > Traits ;

		if( Traits::transpose_can_be_cached && matrix.transposeIndex().valid )
		{
			const int ResSegDim = BlockDims< typename Derived::BlockType, true >::Rows ;
			typedef Segmenter< ResSegDim, ResT, typename Derived::Index > ResSegmenter ;
			ResSegmenter resSegmenter( res, matrix.majorIndex().innerOffsetsData() ) ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
			for( int i = 0 ; i < (int) matrix.transposeIndex().outerSize() ; ++i )
			{
				typename ResSegmenter::ReturnType seg( resSegmenter[ i ] ) ;
				innerRowMultiply< !Traits::transpose_can_be_cached >( matrix.data(), matrix.transposeIndex(), i, rhs, seg, alpha ) ;
			}
		} else {
			Base::multiply( matrix, rhs, res, alpha ) ;
		}
	}
} ;


template < typename Derived >
template < bool Transpose, typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::multiply( const RhsT& rhs, ResT& res, typename RhsT::Scalar alpha, typename ResT::Scalar beta  ) const
{
	if( ( typename ResT::Scalar ) 0 == beta )
		res.setZero() ;
	if( ( typename ResT::Scalar ) 1 != beta )
		res *= beta ;

	SparseBlockMatrixVectorMultiplier
			< Traits::is_symmetric, bool(Transpose) == bool(Traits::is_col_major), Transpose >
			::multiply( *this, rhs, res, alpha )  ;

}


}

#endif
