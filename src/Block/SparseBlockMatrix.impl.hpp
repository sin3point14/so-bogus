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


}

#endif
