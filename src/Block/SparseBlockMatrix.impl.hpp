#ifndef SPARSE_BLOCK_INDEX_IMPL_HPP
#define SPARSE_BLOCK_INDEX_IMPL_HPP


#include "SparseBlockMatrix.hpp"

namespace bogus {

template < typename Derived >
void SparseBlockMatrixBase< Derived >::setRows(
		const std::vector< Index > &rowsPerBlocks )
{
	setInnerOffets( colMajorIndex(), rowsPerBlocks );
	this->m_rows = colMajorIndex().innerOffsets.back() ;
	rowMajorIndex().resizeOuter( rowsPerBlocks.size() ) ;

	if( Traits::is_symmetric ) setCols( rowsPerBlocks ) ;
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::setCols(
		const std::vector< Index > &colsPerBlocks )
{
	setInnerOffets( rowMajorIndex(), colsPerBlocks );
	this->m_cols = rowMajorIndex().innerOffsets.back() ;

	colMajorIndex().resizeOuter( colsPerBlocks.size() ) ;
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::setInnerOffets(
		SparseIndexType &index, const std::vector<Index> &blockSizes)
{
	index.innerOffsets.resize( blockSizes.size() + 1 ) ;
	index.innerOffsets[0] = 0 ;
	for ( unsigned i = 0 ; i < blockSizes.size() ; ++ i )
	{
		index.innerOffsets[ i+1 ] = index.innerOffsets[ i ] + blockSizes[i] ;
	}
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::finalize()
{
	m_majorIndex.finalize() ;
}

template < typename Derived >
bool SparseBlockMatrixBase< Derived >::computeMinorIndex()
{
	// Warning if Traits::is_compressed, the col major index will not be
	// valid unless the transpose Matrix is cached as well

	if ( m_minorIndexValid ) return true ;

	SparseBlockIndex< > cmIndex ;
	computeMinorIndex( cmIndex ) ;

	m_minorIndex.assign( cmIndex ) ;
	m_minorIndexValid = ! Traits::is_compressed ;

	return m_minorIndexValid ;
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::computeMinorIndex(SparseBlockIndex< > &cmIndex) const
{

	cmIndex.resizeOuter( colsOfBlocks() );

	for ( unsigned i = 0 ; i < rowsOfBlocks() ; ++ i )
	{
		// For a symmetric matrix, do not store diagonal block in col-major index
		for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( rowMajorIndex(), i ) ;
			 it && !( Traits::is_symmetric && it.inner() == i ) ; ++ it )
		{
			cmIndex.insertBack( it.inner(), i, it.ptr() );
		}
	}

}

template < typename Derived >
const SparseBlockIndex< >& SparseBlockMatrixBase< Derived >::getUncompressedMinorIndex(SparseBlockIndex< > &cmIndex) const
{
	if( m_minorIndexValid && !m_transposeCached ) return m_minorIndex.asUncompressed() ;
	computeMinorIndex( cmIndex ) ;
	return cmIndex ;
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::cacheTranspose()
{
	if ( m_transposeCached ) return ;

	SparseBlockIndex< > cmIndex ;
	if ( !m_minorIndexValid )
	{
		computeMinorIndex( cmIndex ) ;

		m_minorIndex.assign( cmIndex ) ;
		m_minorIndexValid = ! Traits::is_compressed ;
	}

	this->m_blocks.reserve( 2*m_nBlocks ) ;

	for ( unsigned i = 0 ; i < colsOfBlocks() ; ++ i )
	{
		typename SparseBlockIndex< >::InnerIterator uncompressed_it
			 ( getUncompressedMinorIndex( cmIndex ), i ) ;
		for( typename SparseIndexType::InnerIterator it( m_minorIndex, i ) ;
			 it ; ++ it, ++uncompressed_it )
		{
			allocateBlock() = block( uncompressed_it.ptr() ).transpose() ;
			m_minorIndex.setPtr( it, ( BlockPtr ) ( this->m_blocks.size() - 1 )  ) ;
		}
	}

	m_minorIndexValid   = true ;
	m_transposeCached = true ;
}

template < typename Derived >
SparseBlockIndexBase& SparseBlockMatrixBase< Derived >::getIndex(
		const bool transpose, const bool colWise,
		TransposeMode &indexTransposeMode,
		SparseBlockIndex< >& aux ) const
{
	if ( (colWise ^ transpose) == Traits::is_col_major )
	{
		// Native order
		indexTransposeMode = transpose ? TransposeAll : NoTranspose ;
		return m_majorIndex ;
	} else {
		// Minor order
		if( m_transposeCached )
		{
			indexTransposeMode = transpose ? NoTranspose : TransposeAll ;
			return m_minorIndex ;
		}

		indexTransposeMode = transpose ? TransposeAll : NoTranspose ;
		return getUncompressedMinorIndex( aux ) ;
	}
}

template < typename Derived >
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::multiply( const RhsT& rhs, ResT& res, bool transposed ) const
{
	if( Traits::is_symmetric )
	{
		if( m_transposeCached )
		{
			for( Index i = 0 ; i < rowsOfBlocks() ; ++i )
			{
				typename ResT::SegmentReturnType seg( rowSegment( res, i ) ) ;
				innerRowMultiply( rowMajorIndex(), i, rhs, seg ) ;
				innerRowMultiply( colMajorIndex(), i, rhs, seg ) ;
			}
		} else {
			for( Index i = 0 ; i < rowsOfBlocks() ; ++i )
			{
				typename RhsT::ConstSegmentReturnType rhs_seg( rowSegment( rhs, i ) ) ;
				typename ResT::SegmentReturnType res_seg( rowSegment( res, i ) ) ;
				for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( m_majorIndex, i ) ;
					 it ; ++ it )
				{
					const BlockType &b = block( it.ptr() ) ;
					res_seg += b * colSegment( rhs, it.inner() ) ;
					if( it.inner() != i )
						colSegment( res, it.inner() ) += b.transpose() * rhs_seg  ;
				}
			}
		}

	} else if( transposed )
	{
		if ( Traits::is_col_major )  {
			for( Index i = 0 ; i < colsOfBlocks() ; ++i )
			{
				typename ResT::SegmentReturnType seg( colSegment( res, i ) ) ;
				innerRowTransposedMultiply( colMajorIndex(), i, rhs, seg ) ;
			}
		} else {
			if( m_transposeCached )
			{
				for( Index i = 0 ; i < colsOfBlocks() ; ++i )
				{
					typename ResT::SegmentReturnType seg( colSegment( res, i ) ) ;
					innerRowMultiply( colMajorIndex(), i, rhs, seg ) ;
				}
			} else {
				for( Index i = 0 ; i < rowsOfBlocks() ; ++i )
				{
					typename RhsT::ConstSegmentReturnType seg( rowSegment( rhs, i ) ) ;
					innerColTransposedMultiply( rowMajorIndex(), i, seg, res ) ;
				}
			}
		}
	} else if ( Traits::is_col_major )  {
		for( Index i = 0 ; i < colsOfBlocks() ; ++i )
		{
			typename RhsT::ConstSegmentReturnType seg( colSegment( rhs, i ) ) ;
			innerColMultiply( colMajorIndex(), i, seg, res ) ;
		}
	} else {
		for( Index i = 0 ; i < rowsOfBlocks() ; ++i )
		{
			typename ResT::SegmentReturnType seg( rowSegment( res, i ) ) ;
			innerRowMultiply( rowMajorIndex(), i, rhs, seg ) ;
		}
	}
}

template < typename Derived >
const typename SparseBlockMatrixBase< Derived >::BlockType& SparseBlockMatrixBase< Derived >::diagonal( const Index row ) const
{
	for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( rowMajorIndex(), row ) ;
		 it ; ++ it )
	{
		if( it.inner() == row )
			return this->m_blocks( it.ptr() ) ;
	}
	assert( 0 && "No diagonal block" ) ;
	return this->m_blocks[ (BlockPtr)-1 ] ;
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
		if( m_transposeCached )
		{
			innerRowMultiply( m_minorIndex, row, rhs, res ) ;
		} else {
			assert( m_minorIndexValid ) ;
			innerRowTransposedMultiply( m_minorIndex, row, rhs, res ) ;
		}
	}
}

template < typename Derived >
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::innerRowMultiply( const SparseIndexType &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const
{
	for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( index, outerIdx ) ;
		 it ; ++ it )
	{
		res += block( it.ptr() ) * index.innerSegment( rhs, it.inner() ) ;
	}
}

template < typename Derived >
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::innerRowTransposedMultiply( const SparseIndexType &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const
{
	for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( index, outerIdx ) ;
		 it ; ++ it )
	{
		res += block( it.ptr() ).transpose() * index.innerSegment( rhs, it.inner() ) ;
	}
}

template < typename Derived >
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::innerColMultiply( const SparseIndexType &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const
{
	for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( index, outerIdx ) ;
		 it ; ++ it )
	{
		index.innerSegment( res, it.inner() ) += block( it.ptr() ) * rhs  ;
	}
}

template < typename Derived >
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::innerColTransposedMultiply( const SparseIndexType &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const
{
	for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( index, outerIdx ) ;
		 it ; ++ it )
	{
		index.innerSegment( res, it.inner() ) += block( it.ptr() ).transpose() * rhs  ;
	}
}

template < typename Derived >
std::ostream& operator<<( std::ostream &out, const SparseBlockMatrixBase< Derived > &sbm )
{
	out << " Total rows: " << sbm.rows() << " / cols: " << sbm.cols() << std::endl ;
	for ( unsigned i = 0 ; i < sbm.rowsOfBlocks() ; ++ i )
	{
		out << i << ": " ;
		for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( sbm.rowMajorIndex(), i ) ;
			 it ; ++ it )
		{
			out << "(" << it.inner() << ";" << it.ptr() << ")" ;
		}
		out << std::endl ;
	}
	out << " Blocks (" << sbm.nBlocks() << ")" << std::endl ;
	for ( unsigned i = 0 ; i < sbm.nBlocks() ; ++ i )
	{
		out << sbm.block(i) << std::endl ;
		out << " --- " << std::endl ;
	}
	return out ;
}

template < typename Derived >
template < typename BlockT2 >
void SparseBlockMatrixBase<Derived>::cloneStructure( const SparseBlockMatrix< BlockT2, SparseBlockMatrixBase<Derived>::Traits::flags > &source )
{
	m_nBlocks = source.nBlocks() ;
	rowMajorIndex() = source.rowMajorIndex() ;
	colMajorIndex() = source.colMajorIndex() ;
	m_minorIndexValid = false ;
	m_transposeCached = false ;

	this->m_cols = source.cols() ;
	this->m_rows = source.cols() ;
	this->m_blocks.resize( m_nBlocks ) ;
}

template < typename Derived >
template < typename LhsDerived, typename RhsDerived >
void SparseBlockMatrixBase<Derived>::setFromProduct( const SparseBlockMatrixBase< LhsDerived >& lhs,
					 const SparseBlockMatrixBase< RhsDerived >& rhs,
					 bool lhsTransposed, bool rhsTransposed )
{
	typedef BlockMatrixTraits< LhsDerived > LhsTraits ;
	typedef BlockMatrixTraits< RhsDerived > RhsTraits ;

	TransposeMode transposeModeLhs ; TransposeMode transposeModeRhs ;
	SparseBlockIndex< > auxIndexLhs, auxIndexRhs ;

	assert( ! LhsTraits::is_symmetric ) ;
	assert( ! RhsTraits::is_symmetric ) ;

	const SparseBlockIndexBase &lhsIdx = getIndex( lhsTransposed, true, transposeModeLhs, auxIndexLhs ) ;
	const SparseBlockIndexBase &rhsIdx = getIndex( rhsTransposed, true, transposeModeRhs, auxIndexRhs ) ;

	if( lhsIdx.isCompressed() )
	{
		if( rhsIdx.isCompressed() )
		{
			setFromProduct( lhsIdx.asCompressed(), rhsIdx.asCompressed(),
							lhs.blocks(), rhs.blocks(), transposeModeLhs, transposeModeRhs ) ;
		} else {
			setFromProduct( lhsIdx.asCompressed(), rhsIdx.asUncompressed(),
							lhs.blocks(), rhs.blocks(), transposeModeLhs, transposeModeRhs ) ;
		}
	} else {
		if( rhsIdx.isCompressed() )
		{
			setFromProduct( lhsIdx.asUncompressed(), rhsIdx.asCompressed(),
							lhs.blocks(), rhs.blocks(), transposeModeLhs, transposeModeRhs ) ;
		} else {
			setFromProduct( lhsIdx.asUncompressed(), rhsIdx.asUncompressed(),
							lhs.blocks(), rhs.blocks(), transposeModeLhs, transposeModeRhs ) ;
		}
	}

}

template < typename Derived >
template < typename LhsIndex, typename RhsIndex, typename LhsBlock, typename RhsBlock  >
void SparseBlockMatrixBase<Derived>::setFromProduct(
		const LhsIndex &lhsIdx,
		const RhsIndex &rhsIdx,
		const std::vector< LhsBlock > &lhsData,
		const std::vector< RhsBlock > &rhsData,
		TransposeMode transposeLhs, TransposeMode transposeRhs
		)
{
	(void) lhsIdx ;
	(void) rhsIdx ;
	(void) lhsData ;
	(void) rhsData;
	(void) transposeLhs;
	(void) transposeRhs;
}

template < typename BlockT, int Flags >
class SparseBlockMatrix : public SparseBlockMatrixBase< SparseBlockMatrix< BlockT, Flags > >
{
} ;

}

#endif
