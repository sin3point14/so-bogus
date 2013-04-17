#ifndef BOGUS_SPARSE_MATRIXVECTOR_PRODUCT_HPP
#define BOGUS_SPARSE_MATRIXVECTOR_PRODUCT_HPP

#include "SparseBlockMatrix.hpp"
#include "Expressions.hpp"

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
			innerRowMultiply( m_transposeIndex, row, rhs, res ) ;
		} else {
			assert( m_minorIndex.valid ) ;
			innerRowTransposedMultiply( m_minorIndex, row, rhs, res ) ;
		}
	}
}

template < typename Derived >
template < typename IndexT, typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::innerRowMultiply( const IndexT &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const
{
	for( typename IndexT::InnerIterator it( index, outerIdx ) ;
		 it ; ++ it )
	{
		res += block( it.ptr() ) * index.innerSegment( rhs, it.inner() ) ;
	}
}

template < typename Derived >
template < typename IndexT, typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::innerRowTransposedMultiply( const IndexT &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const
{
	for( typename IndexT::InnerIterator it( index, outerIdx ) ;
		 it ; ++ it )
	{
		res += transpose_block( block( it.ptr() ) ) * index.innerSegment( rhs, it.inner() ) ;
	}
}

template < typename Derived >
template < typename IndexT, typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::innerColMultiply( const IndexT &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const
{
	for( typename IndexT::InnerIterator it( index, outerIdx ) ;
		 it ; ++ it )
	{
		index.innerSegment( res, it.inner() ) += block( it.ptr() ) * rhs  ;
	}
}

template < typename Derived >
template < typename IndexT, typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::innerColTransposedMultiply( const IndexT &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const
{
	for( typename IndexT::InnerIterator it( index, outerIdx ) ;
		 it ; ++ it )
	{
		index.innerSegment( res, it.inner() ) += transpose_block( block( it.ptr() ) ) * rhs  ;
	}
}

template < typename Derived >
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::multiply( const RhsT& rhs, ResT& res, bool transposed ) const
{

	if( Traits::is_symmetric )
	{
		if( m_transposeIndex.valid )
		{
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
			for( Index i = 0 ; i < majorIndex().outerSize() ; ++i )
			{
				typename ResT::SegmentReturnType seg( rowSegment( res, i ) ) ;
				innerRowMultiply( majorIndex(), i, rhs, seg ) ;
				innerRowMultiply( m_transposeIndex, i, rhs, seg ) ;
			}
		} else {
#ifdef BOGUS_DONT_PARALLELIZE
			for( Index i = 0 ; i < majorIndex().outerSize() ; ++i )
			{
				typename RhsT::ConstSegmentReturnType rhs_seg( majorIndex().innerSegment( rhs, i ) ) ;
				typename ResT::SegmentReturnType res_seg( majorIndex().innerSegment( res, i ) ) ;
				for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( m_majorIndex, i ) ;
					 it ; ++ it )
				{
					const BlockType &b = block( it.ptr() ) ;
					res_seg += b * m_minorIndex.innerSegment( rhs, it.inner() ) ;
					if( it.inner() != i )
						m_minorIndex.innerSegment( res, it.inner() ) += transpose_block( b ) * rhs_seg  ;
				}
			}
#else
		multiplyAndReduct( rhs, res, transposed, getBlockProductResVec( res ) ) ;
#endif
		}

	} else if( transposed )
	{
		if ( Traits::is_col_major )  {
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
			for( Index i = 0 ; i < colsOfBlocks() ; ++i )
			{
				typename ResT::SegmentReturnType seg( colSegment( res, i ) ) ;
				innerRowTransposedMultiply( colMajorIndex(), i, rhs, seg ) ;
			}
		} else {
			if( m_transposeIndex.valid )
			{
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
				for( Index i = 0 ; i < colsOfBlocks() ; ++i )
				{
					typename ResT::SegmentReturnType seg( colSegment( res, i ) ) ;
					innerRowMultiply( m_transposeIndex, i, rhs, seg ) ;
				}
			} else {
#ifdef BOGUS_DONT_PARALLELIZE
				for( Index i = 0 ; i < rowsOfBlocks() ; ++i )
				{
					typename RhsT::ConstSegmentReturnType seg( rowSegment( rhs, i ) ) ;
					innerColTransposedMultiply( rowMajorIndex(), i, seg, res ) ;
				}
#else
				multiplyAndReduct( rhs, res, transposed, getBlockProductResVec( res ) ) ;
#endif
			}
		}
	} else if ( Traits::is_col_major )  {
#ifdef BOGUS_DONT_PARALLELIZE
		for( Index i = 0 ; i < colsOfBlocks() ; ++i )
		{
			typename RhsT::ConstSegmentReturnType seg( colSegment( rhs, i ) ) ;
			innerColMultiply( colMajorIndex(), i, seg, res ) ;
		}
#else
		multiplyAndReduct( rhs, res, transposed, getBlockProductResVec( res ) ) ;
#endif
	} else {
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
		for( Index i = 0 ; i < rowsOfBlocks() ; ++i )
		{
			typename ResT::SegmentReturnType seg( rowSegment( res, i ) ) ;
			innerRowMultiply( rowMajorIndex(), i, rhs, seg ) ;
		}
	}
}


template < typename Derived >
template < typename RhsT, typename ResT, typename LocalResT >
void SparseBlockMatrixBase< Derived >::multiplyAndReduct( const RhsT& rhs, ResT& res, bool transpose, const LocalResT& ) const
{
#ifdef BOGUS_DONT_PARALLELIZE
	(void) rhs ; (void) res ; (void) transpose ;
	assert( false && "multiplyAndReduct should never be called if BOGUS_DONT_PARALLELIZE" ) ;
#else
#pragma omp parallel
{
	LocalResT locRes( res.rows() ) ;
	locRes.setZero() ;

	if( Traits::is_symmetric )
	{
#pragma omp for
		for( Index i = 0 ; i < majorIndex().outerSize() ; ++i )
		{
			typename RhsT::ConstSegmentReturnType rhs_seg( majorIndex().innerSegment( rhs, i ) ) ;
			typename ResT::SegmentReturnType res_seg( majorIndex().innerSegment( res, i ) ) ;
			for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( m_majorIndex, i ) ;
				 it ; ++ it )
			{
				const BlockType &b = block( it.ptr() ) ;
				res_seg += b * m_minorIndex.innerSegment( rhs, it.inner() ) ;
				if( it.inner() != i )
					m_minorIndex.innerSegment( locRes, it.inner() ) += transpose_block( b ) * rhs_seg  ;
			}
		}
	} else if ( transpose ) {
#pragma omp for
		for( Index i = 0 ; i < rowsOfBlocks() ; ++i )
		{
			typename RhsT::ConstSegmentReturnType seg( rowSegment( rhs, i ) ) ;
			innerColTransposedMultiply( rowMajorIndex(), i, seg, locRes ) ;
		}
	} else if ( Traits::is_col_major ) {
#pragma omp for
		for( Index i = 0 ; i < colsOfBlocks() ; ++i )
		{
			typename RhsT::ConstSegmentReturnType seg( colSegment( rhs, i ) ) ;
			innerColMultiply( colMajorIndex(), i, seg, locRes ) ;
		}
	}


#pragma omp critical
	res += locRes ;
}
#endif
}


}

#endif
