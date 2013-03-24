#ifndef SPARSE_BLOCK_INDEX_IMPL_HPP
#define SPARSE_BLOCK_INDEX_IMPL_HPP


#include "SparseBlockMatrix.hpp"

namespace bogus {

template < typename Derived >
void SparseBlockMatrixBase< Derived >::setRows(
		const std::vector< Index > &rowsPerBlocks )
{
	m_rowOffsets.resize( rowsPerBlocks.size() + 1 ) ;
	m_rowOffsets[0] = 0 ;
	for ( unsigned i = 0 ; i < rowsPerBlocks.size() ; ++ i )
	{
		m_rowOffsets[ i+1 ] = m_rowOffsets[ i ] + rowsPerBlocks[i] ;
	}

	this->m_rows = m_rowOffsets.back() ;
	m_rowMajorIndex.resizeOuter( rowsPerBlocks.size() ) ;

}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::setCols(
		const std::vector< Index > &colsPerBlocks )
{
	m_colOffsets.resize( colsPerBlocks.size() + 1 ) ;
	m_colOffsets[0] = 0 ;
	for ( unsigned i = 0 ; i < colsPerBlocks.size() ; ++ i )
	{
		m_colOffsets[ i+1 ] = m_colOffsets[ i ] + colsPerBlocks[i] ;
	}
	this->m_cols = m_colOffsets.back() ;
}

template < typename Derived >
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::multiply( const RhsT& rhs, ResT& res, bool transposed ) const
{
	if( transposed )
	{
		for( Index i = 0 ; i < blockRows() ; ++i )
		{
			typename ResT::ConstSegmentReturnType seg( rowSegment( rhs, i ) ) ;
			colMultiply( i, seg, res ) ;
		}
	} else {
		for( Index i = 0 ; i < blockRows() ; ++i )
		{
			typename ResT::SegmentReturnType seg( rowSegment( res, i ) ) ;
			rowMultiply( i, rhs, seg ) ;
		}
	}
}

template < typename Derived >
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::rowMultiply( const Index row, const RhsT& rhs, ResT& res ) const
{
	for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( rowMajorIndex(), row ) ;
		 it ; ++ it )
	{
		res += block( it.ptr() ) * colSegment( rhs, it.inner() ) ;
	}
}

template < typename Derived >
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::colMultiply( const Index col, const RhsT& rhs, ResT& res ) const
{
	for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( rowMajorIndex(), col ) ;
		 it ; ++ it )
	{
		colSegment( res, it.inner() ) += block( it.ptr() ).transpose() * rhs  ;
	}
}

template < typename Derived >
std::ostream& operator<<( std::ostream &out, const SparseBlockMatrixBase< Derived > &sbm )
{
	out << " Total rows: " << sbm.rows() << " / cols: " << sbm.cols() << std::endl ;
	for ( unsigned i = 0 ; i < sbm.blockRows() ; ++ i )
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

template < typename BlockT, bool Compressed >
class SparseBlockMatrix : public SparseBlockMatrixBase< SparseBlockMatrix< BlockT, Compressed > >
{
	typedef typename BlockMatrixTraits< SparseBlockMatrix< BlockT, Compressed > >::SparseIndexType SparseIndexType ;
	typedef typename SparseIndexType::BlockPtr BlockPtr ;
public:
	void finalize() { }

} ;

template < typename BlockT >
class SparseBlockMatrix< BlockT, true > : public SparseBlockMatrixBase< SparseBlockMatrix< BlockT, true > >
{
	typedef typename BlockMatrixTraits< SparseBlockMatrix< BlockT, true > >::SparseIndexType SparseIndexType ;
	typedef typename SparseIndexType::BlockPtr BlockPtr ;
public:

	void finalize() { this->m_rowMajorIndex.finalize() ; }

} ;



}

#endif
