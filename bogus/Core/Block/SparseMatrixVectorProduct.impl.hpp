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

//! Implementation for non-symmetric, in order matrix/vector products
template < bool UseBSRLibrary, bool Symmetric, bool NativeOrder, bool Transpose >
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

//! Implementation for symmetric products
template < bool NativeOrder, bool Transpose >
struct SparseBlockMatrixVectorMultiplier< false, true, NativeOrder, Transpose >
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

//! Implemntation for non-symmetric, out-of-order products
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

//! Implemntation for non-symmetric, out-of-order products
template < bool Transpose >
struct SparseBlockMatrixVectorMultiplier< false, false, false, Transpose >
		: public OutOfOrderSparseBlockMatrixVectorMultiplier< Transpose >
{} ;

//! Implemntation for non-symmetric, out-of-order transposed products using cached transpose
template < >
struct SparseBlockMatrixVectorMultiplier< false, false, false, true >
		: public OutOfOrderSparseBlockMatrixVectorMultiplier< true >
{
	typedef OutOfOrderSparseBlockMatrixVectorMultiplier< true > Base;

	template < typename Derived, typename RhsT, typename ResT, typename ScalarT >
	static void multiply( const SparseBlockMatrixBase< Derived >& matrix,  const RhsT& rhs, ResT& res, ScalarT alpha )
	{
		if( Derived::has_square_or_dynamic_blocks && matrix.transposeIndex().valid )
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
				innerRowMultiply< !Derived::has_square_or_dynamic_blocks >( matrix.data(), matrix.transposeIndex(), i, rhs, seg, alpha ) ;
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

#ifdef BOGUS_WITH_MKL
	const bool use_bsr_library = is_bsr_compatible  ;
#else
	const bool use_bsr_library = false  ;
#endif

	SparseBlockMatrixVectorMultiplier
			< use_bsr_library, Traits::is_symmetric, bool(Transpose) == bool(Traits::is_col_major), Transpose >
			::multiply( *this, rhs, res, alpha )  ;

}

// Split-row-multiply


//! Implementation for for non-symmetric matrices
template < bool Symmetric, bool NativeOrder >
struct SparseBlockSplitRowMultiplier
{
	template < typename Derived,  typename RhsT, typename ResT >
	static void splitRowMultiply( const SparseBlockMatrixBase< Derived >& matrix, typename Derived::Index row, const RhsT& rhs, ResT& res )
	{
		const Segmenter< BlockDims< typename Derived::BlockType, !NativeOrder >::Cols, const RhsT, typename Derived::Index >
				segmenter( rhs, matrix.colOffsets() ) ;

		assert( matrix.rowMajorIndex().valid ) ;

		for( typename SparseBlockMatrixBase< Derived >::RowIndexType::InnerIterator it( matrix.rowMajorIndex(), row ) ;
			 it ; ++ it )
		{
			if( it.inner() != row )
				res += BlockGetter< !NativeOrder >::get( matrix.block( it.ptr() ) ) * segmenter[ it.inner() ] ;
		}
	}
} ;

//! Implementation for symmetric matrices
template < bool NativeOrder >
struct SparseBlockSplitRowMultiplier< true, NativeOrder >
{
	template < typename Derived, typename RhsT, typename ResT >
	static void splitRowMultiply( const SparseBlockMatrixBase< Derived >& matrix, typename Derived::Index row, const RhsT& rhs, ResT& res )
	{
		const Segmenter< BlockDims< typename Derived::BlockType, !NativeOrder >::Cols, const RhsT, typename Derived::Index >
				segmenter( rhs, matrix.colOffsets() ) ;

		for( typename SparseBlockMatrixBase< Derived >::MajorIndexType::InnerIterator it( matrix.majorIndex(), row ) ;
			 it && it.inner() < row ; ++ it )
		{
			res += BlockGetter< !NativeOrder >::get( matrix.block( it.ptr() ) ) * segmenter[ it.inner() ] ;
		}

		if( Derived::has_square_or_dynamic_blocks && matrix.transposeIndex().valid )
		{
			innerRowMultiply< NativeOrder ^ ( (bool) Derived::has_square_or_dynamic_blocks ) >
					( matrix.data(), matrix.transposeIndex(), row, rhs, res, 1 ) ;
		} else {
			assert( matrix.minorIndex().valid ) ;
			innerRowMultiply< NativeOrder >( matrix.data(), matrix.minorIndex(), row, rhs, res, 1 ) ;
		}

	}
} ;


template < typename Derived >
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::splitRowMultiply( const Index row, const RhsT& rhs, ResT& res ) const
{
	SparseBlockSplitRowMultiplier< Traits::is_symmetric, !Traits::is_col_major >
			::splitRowMultiply( *this, row, rhs, res )  ;
}


}

#endif
