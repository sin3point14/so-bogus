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
#include "BlockTranspose.hpp"

namespace bogus {

template < typename BlockT, typename IndexT, typename GetterT, typename RhsT, typename ResT, typename ScalarT >
static void innerRowMultiply( const BlockT* blocks, const IndexT &index, const GetterT& getter,
							const typename IndexT::Index outerIdx, const RhsT& rhs, ResT& res, ScalarT alpha )
{
	for( typename IndexT::InnerIterator it( index, outerIdx ) ; it ; ++ it )
	{
		res += alpha * ( getter.get( blocks[ it.ptr() ] ) * index.innerSegment( rhs, it.inner() ) ) ;
	}
}

template < typename BlockT, typename IndexT, typename GetterT, typename RhsT, typename ResT, typename ScalarT >
static void innerColMultiply( const BlockT* blocks, const IndexT &index,  const GetterT& getter,
							const typename IndexT::Index outerIdx, const RhsT& rhs, ResT& res, ScalarT alpha )
{
	for( typename IndexT::InnerIterator it( index, outerIdx ) ; it ; ++ it )
	{
		index.innerSegment( res, it.inner() ) += alpha * ( getter.get( blocks[ it.ptr() ] ) * rhs ) ;
	}
}

template < typename Derived >
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::splitRowMultiply( const Index row, const RhsT& rhs, ResT& res ) const
{
	assert( Traits::is_symmetric || !Traits::is_col_major ) ;
	BlockGetter< Traits::is_col_major > getter ;
	BlockGetter< !Traits::is_col_major > t_getter ;

	for( typename SparseBlockMatrixBase< Derived >::MajorIndexType::InnerIterator it( m_majorIndex, row ) ;
		 it ; ++ it )
	{
		if( it.inner() != row )
			res += getter.get( block( it.ptr() ) ) * colSegment( rhs, it.inner() ) ;
	}
	if( Traits::is_symmetric )
	{
		if( m_transposeIndex.valid )
		{
			innerRowMultiply( this->data(), m_transposeIndex, getter, row, rhs, res, 1 ) ;
		} else {
			assert( m_minorIndex.valid ) ;
			innerRowMultiply( this->data(), m_minorIndex, t_getter, row, rhs, res, 1 ) ;
		}
	}
}

template < bool Symmetric, bool NativeOrder, bool Transpose >
struct SparseBlockMatrixVectorMultiplier
{
	template < typename Derived, typename RhsT, typename ResT, typename ScalarT >
	static void multiply( const SparseBlockMatrixBase< Derived >& matrix,  const RhsT& rhs, ResT& res, ScalarT alpha )
	{
		BlockGetter< Transpose > getter ;
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
		for( int i = 0 ; i < (int) matrix.majorIndex().outerSize() ; ++i )
		{
			typename ResT::SegmentReturnType seg( matrix.minorIndex().innerSegment( res, i ) ) ;
			innerRowMultiply( matrix.data(), matrix.majorIndex(), getter, i, rhs, seg, alpha ) ;
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
		typedef typename SparseBlockMatrixBase< Derived >::MajorIndexType MajorIndexType ;
#pragma omp parallel
		{
			LocalResT locRes( res.rows() ) ;
			locRes.setZero() ;

#pragma omp for
			for( int i = 0 ; i < (int) matrix.majorIndex().outerSize() ; ++i )
			{
				typename RhsT::ConstSegmentReturnType rhs_seg( matrix.majorIndex().innerSegment( rhs, i ) ) ;
				typename ResT::SegmentReturnType locRes_seg( matrix.majorIndex().innerSegment( locRes, i ) ) ;
				for( typename MajorIndexType::InnerIterator it( matrix.majorIndex(), i ) ;
					 it ; ++ it )
				{
					const typename Derived::BlockType &b = matrix.block( it.ptr() ) ;
					locRes_seg += alpha * ( b * matrix.minorIndex().innerSegment( rhs, it.inner() ) ) ;
					if( it.inner() != (typename MajorIndexType::Index) i )
						matrix.minorIndex().innerSegment( locRes, it.inner() ) += alpha * ( transpose_block( b ) * rhs_seg ) ;
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
		typedef typename SparseBlockMatrixBase< Derived >::Index Index ;
		if( matrix.transposeIndex().valid )
		{
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
			for( int i = 0 ; i < (int) matrix.majorIndex().outerSize() ; ++i )
			{
				typename ResT::SegmentReturnType seg( matrix.majorIndex().innerSegment( res, i ) ) ;
				innerRowMultiply( matrix.data(), matrix.majorIndex()    , BlockGetter< false >(), i, rhs, seg, alpha ) ;
				innerRowMultiply( matrix.data(), matrix.transposeIndex(), BlockGetter< false >(), i, rhs, seg, alpha ) ;
			}
		} else {
#ifdef BOGUS_DONT_PARALLELIZE
			for( Index i = 0 ; i < matrix.majorIndex().outerSize() ; ++i )
			{
				typename RhsT::ConstSegmentReturnType rhs_seg( matrix.majorIndex().innerSegment( rhs, i ) ) ;
				typename ResT::SegmentReturnType res_seg( matrix.majorIndex().innerSegment( res, i ) ) ;
				for( typename SparseBlockMatrixBase< Derived >::MajorIndexType::InnerIterator it( matrix.majorIndex(), i ) ;
					 it ; ++ it )
				{
					const typename Derived::BlockType &b = matrix.block( it.ptr() ) ;
					res_seg += alpha * ( b * matrix.minorIndex().innerSegment( rhs, it.inner() ) ) ;
					if( it.inner() != i )
						matrix.minorIndex().innerSegment( res, it.inner() ) += alpha * ( transpose_block( b ) * rhs_seg ) ;
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
		BlockGetter< Transpose > getter ;

#pragma omp parallel
		{
			LocalResT locRes( res.rows() ) ;
			locRes.setZero() ;

#pragma omp for
			for( int i = 0 ; i < (int) matrix.majorIndex().outerSize() ; ++i )
			{
				typename RhsT::ConstSegmentReturnType seg( matrix.minorIndex().innerSegment( rhs, i ) ) ;
				innerColMultiply( matrix.data(), matrix.majorIndex(), getter, i, seg, locRes, alpha ) ;
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
#ifdef BOGUS_DONT_PARALLELIZE
		for( Index i = 0 ; i < matrix.majorIndex().outerSize() ; ++i )
		{
			typename RhsT::ConstSegmentReturnType seg( matrix.minorIndex().innerSegment( rhs, i ) ) ;
			innerColMultiply( matrix.data(), matrix.majorIndex(), BlockGetter< Transpose >(), i, seg, res, alpha ) ;
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
		if( matrix.transposeIndex().valid )
		{
			BlockGetter< false > getter ;
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
			for( int i = 0 ; i < (int) matrix.transposeIndex().outerSize() ; ++i )
			{
				typename ResT::SegmentReturnType seg( matrix.majorIndex().innerSegment( res, i ) ) ;
				innerRowMultiply( matrix.data(), matrix.transposeIndex(), getter, i, rhs, seg, alpha ) ;
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
