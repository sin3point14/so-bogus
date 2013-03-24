#ifndef SPARSEBLOCKMATRIX_HH
#define SPARSEBLOCKMATRIX_HH

#include "BlockMatrix.hpp"
#include "SparseBlockIndex.hpp"

#include <iosfwd>

namespace bogus
{

template < typename Derived >
class SparseBlockMatrixBase : public BlockMatrixBase< Derived >
{

public:
	typedef typename BlockMatrixTraits< Derived >::SparseIndexType SparseIndexType ;
	typedef typename SparseIndexType::BlockPtr BlockPtr ;
	typedef typename BlockMatrixBase< Derived >::BlockType BlockType ;

	SparseBlockMatrixBase()
		: m_nBlocks(0), m_rowMajorIndex( m_colOffsets, m_rowOffsets )
	{}

	void setRows( const std::vector< Index > &rowsPerBlocks ) ;
	void setRows( Index n_blocks, Index rows_per_block )
	{
		setRows( std::vector< Index >( n_blocks, rows_per_block ) ) ;
	}

	void setCols( const std::vector< Index > &colsPerBlocks ) ;
	void setCols( Index n_blocks, Index cols_per_block )
	{
		setCols( std::vector< Index >( n_blocks, cols_per_block ) ) ;
	}

	Index rowsOfBlocks() const { return m_rowOffsets.size() - 1 ; }
	Index colsOfBlocks() const { return m_colOffsets.size() - 1 ; }
	Index blockRows( Index row ) const { return m_rowOffsets[ row + 1 ] - m_rowOffsets[ row ] ; }
	Index blockCols( Index col ) const { return m_colOffsets[ col + 1 ] - m_colOffsets[ col ] ; }

	void reserve( std::size_t nBlocks )
	{
		this->m_blocks.reserve( nBlocks ) ;
	}

	void prealloc( std::size_t nBlocks )
	{
		this->m_blocks.resize( nBlocks ) ;
		m_nBlocks = nBlocks ;
	}

	BlockType& insertBack( Index row, Index col )
	{
		m_rowMajorIndex.insertBack( row, col, m_nBlocks ) ;
		return allocateBlock() ;
	}

	void finalize()
	{
		this->derived().finalize() ;
	}

	std::size_t nBlocks() const { return m_nBlocks ; }

	const BlockType& block( BlockPtr ptr ) const
	{
		return this->m_blocks[ ptr ] ;
	}

	BlockType& block( BlockPtr ptr )
	{
		return this->m_blocks[ ptr ] ;
	}

	const SparseIndexType& rowMajorIndex() const
	{
		return m_rowMajorIndex ;
	}

	template < typename RhsT, typename ResT >
	void multiply( const RhsT& rhs, ResT& res, bool transposed = false ) const ;

	template < typename RhsT, typename ResT >
	void splitRowMultiply( const Index row, const RhsT& rhs, ResT& res, bool transposed = false ) const ;

protected:
	BlockType& allocateBlock()
	{
		++m_nBlocks ;
		this->m_blocks.push_back( BlockType() ) ;
		return this->m_blocks.back() ;
	}

	template < typename VecT >
	typename VecT::SegmentReturnType rowSegment( VecT& v, Index rowBlockIdx ) const
	{
		return v.segment( m_rowOffsets[ rowBlockIdx ],
						m_rowOffsets[ rowBlockIdx + 1 ] - m_rowOffsets[ rowBlockIdx ] ) ;
	}

	template < typename VecT >
	typename VecT::ConstSegmentReturnType rowSegment( const VecT& v, Index rowBlockIdx ) const
	{
		return v.segment( m_rowOffsets[ rowBlockIdx ],
						  blockRows( rowBlockIdx ) ) ;
	}

	template < typename VecT >
	typename VecT::SegmentReturnType colSegment( VecT& v, Index colBlockIdx ) const
	{
		return v.segment( m_colOffsets[ colBlockIdx ],
						 blockCols( colBlockIdx ) ) ;
	}

	template < typename VecT >
	typename VecT::ConstSegmentReturnType colSegment( const VecT& v, Index colBlockIdx ) const
	{
		return v.segment( m_colOffsets[ colBlockIdx ],
						m_colOffsets[ colBlockIdx + 1 ] - m_colOffsets[ colBlockIdx ] ) ;
	}

	template < typename RhsT, typename ResT >
	void rowMultiply( const Index row, const RhsT& rhs, ResT& res ) const ;

	template < typename RhsT, typename ResT >
	void colMultiply( const Index col, const RhsT& rhs, ResT& res ) const ;

	std::size_t m_nBlocks ;
	std::vector< Index > m_rowOffsets ;
	std::vector< Index > m_colOffsets ;

	SparseIndexType m_rowMajorIndex ;

} ;

template < typename Derived >
std::ostream& operator<<( std::ostream &out, const SparseBlockMatrixBase< Derived > &sbm ) ;

template < typename BlockT, bool Compressed = false >
class SparseBlockMatrix  ;

template < typename BlockT, bool Compressed >
struct BlockMatrixTraits< SparseBlockMatrix< BlockT, Compressed > >
{
	typedef BlockT BlockType ;
	typedef SparseBlockIndex< BlockT, Compressed > SparseIndexType ;
	typedef typename BlockMatrixTraits< SparseIndexType >::Index Index ;

	struct is_compressed {
		enum { value = Compressed } ;
	} ;
} ;

}

#endif // SPARSEBLOCKMATRIX_HH
