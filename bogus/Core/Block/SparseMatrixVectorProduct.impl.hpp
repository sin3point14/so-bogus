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

template < typename BlockT, typename IndexT, typename GetterT, typename RhsT, typename ResT >
static void innerRowMultiply( const BlockT* blocks, const IndexT &index, const GetterT& getter,
							const typename IndexT::Index outerIdx, const RhsT& rhs, ResT& res ) 
{
	for( typename IndexT::InnerIterator it( index, outerIdx ) ; it ; ++ it )
	{
		res += getter.get( blocks[ it.ptr() ] ) * index.innerSegment( rhs, it.inner() ) ;
	}
}

template < typename BlockT, typename IndexT, typename GetterT, typename RhsT, typename ResT >
static void innerColMultiply( const BlockT* blocks, const IndexT &index,  const GetterT& getter,
							const typename IndexT::Index outerIdx, const RhsT& rhs, ResT& res ) 
{
	for( typename IndexT::InnerIterator it( index, outerIdx ) ; it ; ++ it )
	{
		index.innerSegment( res, it.inner() ) += getter.get( blocks[ it.ptr() ] ) * rhs  ;
	}
}

template < typename Derived >
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::splitRowMultiply( const Index row, const RhsT& rhs, ResT& res ) const
{
	assert( Traits::is_symmetric || !Traits::is_col_major ) ;

	for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( m_majorIndex, row ) ;
		 it ; ++ it )
	{
		if( it.inner() != row )
			res += block( it.ptr() ) * colSegment( rhs, it.inner() ) ;
	}
	if( Traits::is_symmetric )
	{
		if( m_transposeIndex.valid )
		{
			BlockGetter< false > getter ;
			innerRowMultiply( this->data(), m_transposeIndex, getter, row, rhs, res ) ;
		} else {
			BlockGetter< true > getter ;
			assert( m_minorIndex.valid ) ;
			innerRowMultiply( this->data(), m_minorIndex, getter, row, rhs, res ) ;
		}
	}
}

template < bool Symmetric, bool NativeOrder, bool Transpose >
struct SparseBlockMatrixVectorMultiplier
{
	template < typename Derived, typename RhsT, typename ResT >
	static void multiply( const SparseBlockMatrixBase< Derived >& matrix,  const RhsT& rhs, ResT& res )
	{
		BlockGetter< Transpose > getter ;
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
		for( int i = 0 ; i < (int) matrix.majorIndex().outerSize() ; ++i )
		{
			typename ResT::SegmentReturnType seg( matrix.minorIndex().innerSegment( res, i ) ) ;
			innerRowMultiply( matrix.data(), matrix.majorIndex(), getter, i, rhs, seg ) ;
		}
	}
} ;

template < bool NativeOrder, bool Transpose >
struct SparseBlockMatrixVectorMultiplier< true, NativeOrder, Transpose >
{
#ifndef BOGUS_DONT_PARALLELIZE
	template < typename Derived, typename RhsT, typename ResT, typename LocalResT >
	static void multiplyAndReduct( const SparseBlockMatrixBase< Derived >& matrix,  const RhsT& rhs, ResT& res, const LocalResT& )
	{
		typedef typename SparseBlockMatrixBase< Derived >::SparseIndexType SparseIndexType ;
#pragma omp parallel
		{
			LocalResT locRes( res.rows() ) ;
			locRes.setZero() ;

#pragma omp for
			for( int i = 0 ; i < (int) matrix.majorIndex().outerSize() ; ++i )
			{
				typename RhsT::ConstSegmentReturnType rhs_seg( matrix.majorIndex().innerSegment( rhs, i ) ) ;
				typename ResT::SegmentReturnType locRes_seg( matrix.majorIndex().innerSegment( locRes, i ) ) ;
				for( typename SparseIndexType::InnerIterator it( matrix.majorIndex(), i ) ;
					 it ; ++ it )
				{
					const typename Derived::BlockType &b = matrix.block( it.ptr() ) ;
					locRes_seg += b * matrix.minorIndex().innerSegment( rhs, it.inner() ) ;
					if( it.inner() != (typename SparseIndexType::Index) i )
						matrix.minorIndex().innerSegment( locRes, it.inner() ) += transpose_block( b ) * rhs_seg  ;
				}
			}

#pragma omp critical
			res += locRes ;
		}
	}
#endif

	template < typename Derived, typename RhsT, typename ResT >
	static void multiply( const SparseBlockMatrixBase< Derived >& matrix,  const RhsT& rhs, ResT& res )
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
				innerRowMultiply( matrix.data(), matrix.majorIndex()    , BlockGetter< false >(), i, rhs, seg ) ;
				innerRowMultiply( matrix.data(), matrix.transposeIndex(), BlockGetter< false >(), i, rhs, seg ) ;
			}
		} else {
#ifdef BOGUS_DONT_PARALLELIZE
			for( Index i = 0 ; i < matrix.majorIndex().outerSize() ; ++i )
			{
				typename RhsT::ConstSegmentReturnType rhs_seg( matrix.majorIndex().innerSegment( rhs, i ) ) ;
				typename ResT::SegmentReturnType res_seg( matrix.majorIndex().innerSegment( res, i ) ) ;
				for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( matrix.majorIndex(), i ) ;
					 it ; ++ it )
				{
					const typename Derived::BlockType &b = matrix.block( it.ptr() ) ;
					res_seg += b * matrix.minorIndex().innerSegment( rhs, it.inner() ) ;
					if( it.inner() != i )
						matrix.minorIndex().innerSegment( res, it.inner() ) += transpose_block( b ) * rhs_seg  ;
				}
			}
#else
			multiplyAndReduct( matrix, rhs, res, getBlockProductResVec( res ) ) ;
#endif
		}
	}
} ;

template < bool Transpose >
struct OutOfOrderSparseBlockMatrixVectorMultiplier
{

#ifndef BOGUS_DONT_PARALLELIZE
	template < typename Derived, typename RhsT, typename ResT, typename LocalResT >
	static void multiplyAndReduct( const SparseBlockMatrixBase< Derived >& matrix,  const RhsT& rhs, ResT& res, const LocalResT& )
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
				innerColMultiply( matrix.data(), matrix.majorIndex(), getter, i, seg, locRes ) ;
			}

#pragma omp critical
			res += locRes ;
		}
	}
#endif

	template < typename Derived, typename RhsT, typename ResT >
	static void multiply( const SparseBlockMatrixBase< Derived >& matrix,  const RhsT& rhs, ResT& res )
	{
		typedef typename SparseBlockMatrixBase< Derived >::Index Index ;
#ifdef BOGUS_DONT_PARALLELIZE
		for( Index i = 0 ; i < matrix.majorIndex().outerSize() ; ++i )
		{
			typename RhsT::ConstSegmentReturnType seg( matrix.minorIndex().innerSegment( rhs, i ) ) ;
			innerColMultiply( matrix.data(), matrix.majorIndex(), BlockGetter< Transpose >(), i, seg, res ) ;
		}
#else
		multiplyAndReduct( matrix, rhs, res, getBlockProductResVec( res ) ) ;
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

	template < typename Derived, typename RhsT, typename ResT >
	static void multiply( const SparseBlockMatrixBase< Derived >& matrix,  const RhsT& rhs, ResT& res )
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
				innerRowMultiply( matrix.data(), matrix.transposeIndex(), getter, i, rhs, seg ) ;
			}
		} else {
			Base::multiply( matrix, rhs, res ) ;
		}
	}
} ;


template < typename Derived >
template < bool Transpose, typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::multiply( const RhsT& rhs, ResT& res ) const
{

	SparseBlockMatrixVectorMultiplier
			< Traits::is_symmetric, bool(Transpose) == bool(Traits::is_col_major), Transpose >
			::multiply( *this, rhs, res )  ;

}


}

#endif
