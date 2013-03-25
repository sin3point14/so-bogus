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
	SparseBlockIndex< BlockType > cmIndex ;
	computeColMajorIndex( cmIndex ) ;

	return m_colMajorValid ;
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::computeColMajorIndex(SparseBlockIndex< BlockType > &cmIndex)
{
	if ( m_colMajorValid ) return ;

	// Warning if Traits::is_compressed, the col major index will not be
	// valid unless the transpose Matrix is cached as well

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

	m_colMajorIndex.assign( cmIndex ) ;
	m_colMajorValid = ! Traits::is_compressed ;
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::cacheTranspose()
{
	if ( m_transposeCached ) return ;

	SparseBlockIndex< BlockType > cmIndex ;
	computeColMajorIndex( cmIndex ) ;

	this->m_blocks.reserve( 2*m_nBlocks ) ;

	for ( unsigned i = 0 ; i < colsOfBlocks() ; ++ i )
	{
		typename SparseBlockIndex< BlockType >::InnerIterator uncompressed_it
			 ( m_colMajorValid ? m_colMajorIndex.asUncompressed() : cmIndex, i ) ;
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

template < typename BlockT, int Flags >
class SparseBlockMatrix : public SparseBlockMatrixBase< SparseBlockMatrix< BlockT, Flags > >
{
} ;

}

#endif
