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
#include "BlockUtils.hpp"

namespace bogus {

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
			innerRowMultiply( m_transposeIndex, getter, row, rhs, res ) ;
		} else {
			BlockGetter< true > getter ;
			assert( m_minorIndex.valid ) ;
			innerRowMultiply( m_minorIndex, getter, row, rhs, res ) ;
		}
	}
}

template < typename Derived >
template < typename IndexT, typename GetterT, typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::innerRowMultiply( const IndexT &index, const GetterT& getter,
														 const Index outerIdx, const RhsT& rhs, ResT& res ) const
{
	for( typename IndexT::InnerIterator it( index, outerIdx ) ;
		 it ; ++ it )
	{
		res += getter.get( block( it.ptr() ) ) * index.innerSegment( rhs, it.inner() ) ;
	}
}

template < typename Derived >
template < typename IndexT, typename GetterT, typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::innerColMultiply( const IndexT &index,  const GetterT& getter,
														 const Index outerIdx, const RhsT& rhs, ResT& res ) const
{
	for( typename IndexT::InnerIterator it( index, outerIdx ) ;
		 it ; ++ it )
	{
		index.innerSegment( res, it.inner() ) += getter.get( block( it.ptr() ) ) * rhs  ;
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
			matrix.innerRowMultiply( matrix.majorIndex(), getter, i, rhs, seg ) ;
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
#pragma omp parallel
		{
			LocalResT locRes( res.rows() ) ;
			locRes.setZero() ;

#pragma omp for
			for( int i = 0 ; i < (int) matrix.majorIndex().outerSize() ; ++i )
			{
				typename RhsT::ConstSegmentReturnType rhs_seg( matrix.majorIndex().innerSegment( rhs, i ) ) ;
				typename ResT::SegmentReturnType locRes_seg( matrix.majorIndex().innerSegment( locRes, i ) ) ;
				for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( matrix.majorIndex(), i ) ;
					 it ; ++ it )
				{
					const typename Derived::BlockType &b = matrix.block( it.ptr() ) ;
					locRes_seg += b * matrix.minorIndex().innerSegment( rhs, it.inner() ) ;
					if( it.inner() != (Index) i )
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
		if( matrix.transposeIndex().valid )
		{
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
			for( int i = 0 ; i < (int) matrix.majorIndex().outerSize() ; ++i )
			{
				typename ResT::SegmentReturnType seg( matrix.rowSegment( res, i ) ) ;
				matrix.innerRowMultiply( matrix.majorIndex()    , BlockGetter< false >(), i, rhs, seg ) ;
				matrix.innerRowMultiply( matrix.transposeIndex(), BlockGetter< false >(), i, rhs, seg ) ;
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
				matrix.innerColMultiply( matrix.majorIndex(), getter, i, seg, locRes ) ;
			}

#pragma omp critical
			res += locRes ;
		}
	}
#endif

	template < typename Derived, typename RhsT, typename ResT >
	static void multiply( const SparseBlockMatrixBase< Derived >& matrix,  const RhsT& rhs, ResT& res )
	{
#ifdef BOGUS_DONT_PARALLELIZE
		for( Index i = 0 ; i < matrix.majorIndex().outerSize() ; ++i )
		{
			typename RhsT::ConstSegmentReturnType seg( matrix.minorIndex().innerSegment( rhs, i ) ) ;
			matrix.innerColMultiply( matrix.majorIndex(), BlockGetter< Transpose >(), i, seg, res ) ;
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
				matrix.innerRowMultiply( matrix.transposeIndex(), getter, i, rhs, seg ) ;
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
