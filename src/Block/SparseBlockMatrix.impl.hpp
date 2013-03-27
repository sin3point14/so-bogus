#ifndef SPARSE_BLOCK_INDEX_IMPL_HPP
#define SPARSE_BLOCK_INDEX_IMPL_HPP


#include "SparseBlockMatrix.hpp"

namespace bogus {

template < typename Derived >
void SparseBlockMatrixBase< Derived >::setRows(
		const std::vector< Index > &rowsPerBlocks )
{
	setInnerOffets( m_colMajorIndex, rowsPerBlocks );
	this->m_rows = m_colMajorIndex.innerOffsets.back() ;
	m_rowMajorIndex.resizeOuter( rowsPerBlocks.size() ) ;

	if( Traits::is_symmetric ) setCols( rowsPerBlocks ) ;
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::setCols(
		const std::vector< Index > &colsPerBlocks )
{
	setInnerOffets( m_rowMajorIndex, colsPerBlocks );
	this->m_cols = m_rowMajorIndex.innerOffsets.back() ;

	m_colMajorIndex.resizeOuter( colsPerBlocks.size() ) ;
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
	m_rowMajorIndex.finalize() ;
}

template < typename Derived >
bool SparseBlockMatrixBase< Derived >::computeColMajorIndex()
{
	// Warning if Traits::is_compressed, the col major index will not be
	// valid unless the transpose Matrix is cached as well

	if ( m_colMajorValid ) return true ;

	SparseBlockIndex< > cmIndex ;
	computeColMajorIndex( cmIndex ) ;

	m_colMajorIndex.assign( cmIndex ) ;
	m_colMajorValid = ! Traits::is_compressed ;

	return m_colMajorValid ;
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::computeColMajorIndex(SparseBlockIndex< > &cmIndex) const
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
const SparseBlockIndex< >& SparseBlockMatrixBase< Derived >::getUncompressedColMajorIndex(SparseBlockIndex< > &cmIndex) const
{
	if( m_colMajorValid && !m_transposeCached ) return m_colMajorIndex.asUncompressed() ;
	computeColMajorIndex( cmIndex ) ;
	return cmIndex ;
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::cacheTranspose()
{
	if ( m_transposeCached ) return ;

	SparseBlockIndex< > cmIndex ;
	if ( !m_colMajorValid )
	{
		computeColMajorIndex( cmIndex ) ;

		m_colMajorIndex.assign( cmIndex ) ;
		m_colMajorValid = ! Traits::is_compressed ;
	}

	this->m_blocks.reserve( 2*m_nBlocks ) ;

	for ( unsigned i = 0 ; i < colsOfBlocks() ; ++ i )
	{
		typename SparseBlockIndex< >::InnerIterator uncompressed_it
			 ( getUncompressedColMajorIndex( cmIndex ), i ) ;
		for( typename SparseIndexType::InnerIterator it( m_colMajorIndex, i ) ;
			 it ; ++ it, ++uncompressed_it )
		{
			allocateBlock() = block( uncompressed_it.ptr() ).transpose() ;
			m_colMajorIndex.setPtr( it, ( BlockPtr ) ( this->m_blocks.size() - 1 )  ) ;
		}
	}

	m_colMajorValid   = true ;
	m_transposeCached = true ;
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
				innerMultiply( rowMajorIndex(), i, rhs, seg ) ;
				innerMultiply( m_colMajorIndex, i, rhs, seg ) ;
			}
		} else {
			for( Index i = 0 ; i < rowsOfBlocks() ; ++i )
			{
				typename RhsT::ConstSegmentReturnType rhs_seg( rowSegment( rhs, i ) ) ;
				typename ResT::SegmentReturnType res_seg( rowSegment( res, i ) ) ;
				for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( rowMajorIndex(), i ) ;
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
		if( m_transposeCached )
		{
			for( Index i = 0 ; i < colsOfBlocks() ; ++i )
			{
				typename ResT::SegmentReturnType seg( colSegment( res, i ) ) ;
				innerMultiply( m_colMajorIndex, i, rhs, seg ) ;
			}
		} else {
			for( Index i = 0 ; i < rowsOfBlocks() ; ++i )
			{
				typename RhsT::ConstSegmentReturnType seg( rowSegment( rhs, i ) ) ;
				innerTransposedMultiply( rowMajorIndex(), i, seg, res ) ;
			}
		}
	} else {
		for( Index i = 0 ; i < rowsOfBlocks() ; ++i )
		{
			typename ResT::SegmentReturnType seg( rowSegment( res, i ) ) ;
			innerMultiply( rowMajorIndex(), i, rhs, seg ) ;
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
	for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( rowMajorIndex(), row ) ;
		 it ; ++ it )
	{
		if( it.inner() != row )
			res += block( it.ptr() ) * colSegment( rhs, it.inner() ) ;
	}
	if( Traits::is_symmetric )
	{
		if( m_transposeCached )
		{
			innerMultiply( m_colMajorIndex, row, rhs, res ) ;
		} else {
			assert( m_colMajorValid ) ;
			for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( m_colMajorIndex, row ) ;
				 it ; ++ it )
			{
				res  += block( it.ptr() ).transpose() * m_colMajorIndex.innerSegment( rhs, it.inner() )  ;
			}
		}
	}
}

template < typename Derived >
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::innerMultiply( const SparseIndexType &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const
{
	for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( index, outerIdx ) ;
		 it ; ++ it )
	{
		res += block( it.ptr() ) * index.innerSegment( rhs, it.inner() ) ;
	}
}

template < typename Derived >
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::innerTransposedMultiply( const SparseIndexType &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const
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
	m_rowMajorIndex = source.rowMajorIndex() ;
	m_colMajorIndex = source.colMajorIndex() ;
	m_colMajorValid = false ;
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
	SparseBlockIndex< > cmIndex ;

	assert( ! LhsTraits::is_symmetric ) ;
	assert( ! RhsTraits::is_symmetric ) ;

	if( !lhsTransposed )
	{
		if( rhsTransposed )
		{
			setFromProduct( lhs.rowMajorIndex(), rhs.rowMajorIndex(),
							lhs.blocks(), rhs.blocks(), false, true ) ;
		} else {
			if ( rhs.m_transposeCached )
			{
				setFromProduct( lhs.rowMajorIndex(), rhs.colMajorIndex(),
								lhs.blocks(), rhs.blocks(), false, true ) ;
			} else {
				setFromProduct( lhs.rowMajorIndex(), getUncompressedColMajorIndex(cmIndex),
								lhs.blocks(), rhs.blocks(), false, false ) ;
			}
		}
	} else if( lhs.m_transposeCached ) {
		if( rhsTransposed )
		{
			setFromProduct( lhs.colMajorIndex(), rhs.rowMajorIndex(),
							lhs.blocks(), rhs.blocks(), false, true ) ;
		} else {
			if ( rhs.m_transposeCached )
			{
				setFromProduct( lhs.colMajorIndex(), rhs.colMajorIndex(),
								lhs.blocks(), rhs.blocks(), false, true ) ;
			} else {
				setFromProduct( lhs.colMajorIndex(), getUncompressedColMajorIndex(cmIndex),
								lhs.blocks(), rhs.blocks(), false, false ) ;
			}
		}
	} else {
		SparseBlockIndex< > cmIndexL ;
		SparseBlockIndex< > &lhsIndex = getUncompressedColMajorIndex( cmIndexL ) ;

		if( rhsTransposed )
		{
			setFromProduct( lhsIndex, rhs.rowMajorIndex(),
							lhs.blocks(), rhs.blocks(), true, true ) ;
		} else {
			if ( rhs.m_transposeCached )
			{
				setFromProduct( lhsIndex, rhs.colMajorIndex(),
								lhs.blocks(), rhs.blocks(), true, true ) ;
			} else {
				setFromProduct( lhsIndex, getUncompressedColMajorIndex(cmIndex),
								lhs.blocks(), rhs.blocks(), true, false ) ;
			}
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
		bool transposeLhs, bool transposeRhs
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
