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
	typedef BlockMatrixTraits< Derived > Traits ;
	typedef typename Traits::SparseIndexType SparseIndexType ;
	typedef typename Traits::BlockType BlockType ;
	typedef typename SparseIndexType::BlockPtr BlockPtr ;

	SparseBlockMatrixBase()
		: m_nBlocks(0),
		  m_rowMajorIndex( m_colOffsets ),
		  m_colMajorComputed( false ), m_colMajorIndex( m_rowOffsets )
	{}

	void setRows( const std::vector< Index > &rowsPerBlocks ) ;
	void setRows( Index n_blocks, Index rows_per_block )
	{
		setRows( std::vector< Index >( n_blocks, rows_per_block ) ) ;
	}
	void setRows( Index n_blocks )
	{
		setRows( n_blocks, BlockType::RowsAtCompileTime ) ;
	}

	void setCols( const std::vector< Index > &colsPerBlocks ) ;
	void setCols( Index n_blocks, Index cols_per_block )
	{
		setCols( std::vector< Index >( n_blocks, cols_per_block ) ) ;
	}
	void setCols( Index n_blocks )
	{
		setCols( n_blocks, BlockType::ColsAtCompileTime ) ;
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
		m_rowMajorIndex.insertBack( row, col, m_nBlocks++ ) ;
		return allocateBlock() ;
	}

	void finalize() ;
	// Only relevant if Traits::is_symmetric
	void computeColMajorIndex() ;

	std::size_t nBlocks() const { return m_nBlocks ; }

	const BlockType& block( BlockPtr ptr ) const
	{
		return this->m_blocks[ ptr ] ;
	}

	BlockType& block( BlockPtr ptr )
	{
		return this->m_blocks[ ptr ] ;
	}

	const BlockType& diagonal( const Index row ) const ;

	const SparseIndexType& rowMajorIndex() const
	{
		return m_rowMajorIndex ;
	}

	template < typename RhsT, typename ResT >
	void multiply( const RhsT& rhs, ResT& res, bool transposed = false ) const ;

	template < typename RhsT, typename ResT >
	void splitRowMultiply( const Index row, const RhsT& rhs, ResT& res ) const ;

protected:
	BlockType& allocateBlock()
	{
		this->m_blocks.push_back( BlockType() ) ;
		return this->m_blocks.back() ;
	}

	template < typename VecT >
	typename VecT::SegmentReturnType rowSegment( VecT& v, Index rowBlockIdx ) const
	{
		return v.segment( m_rowOffsets[ rowBlockIdx ], blockRows( rowBlockIdx ) ) ;
	}

	template < typename VecT >
	typename VecT::ConstSegmentReturnType rowSegment( const VecT& v, Index rowBlockIdx ) const
	{
		return v.segment( m_rowOffsets[ rowBlockIdx ], blockRows( rowBlockIdx ) ) ;
	}

	template < typename VecT >
	typename VecT::SegmentReturnType colSegment( VecT& v, Index colBlockIdx ) const
	{
		return v.segment( m_colOffsets[ colBlockIdx ], blockCols( colBlockIdx ) ) ;
	}

	template < typename VecT >
	typename VecT::ConstSegmentReturnType colSegment( const VecT& v, Index colBlockIdx ) const
	{
		return v.segment( m_colOffsets[ colBlockIdx ], blockCols( colBlockIdx ) ) ;
	}

	template < typename RhsT, typename ResT >
	void innerMultiply( const SparseIndexType &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const ;

	template < typename RhsT, typename ResT >
	void innerTransposedMultiply( const SparseIndexType &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const ;

	std::size_t m_nBlocks ;
	std::vector< Index > m_rowOffsets ;
	std::vector< Index > m_colOffsets ;

	SparseIndexType m_rowMajorIndex ;

	bool m_colMajorComputed ;
	// For a symmetric matrix, do not store diagonal block in col-major index
	SparseIndexType m_colMajorIndex ;
} ;

template < typename Derived >
std::ostream& operator<<( std::ostream &out, const SparseBlockMatrixBase< Derived > &sbm ) ;

template < typename BlockT, int Flags = BlockMatrixFlags::NONE >
class SparseBlockMatrix  ;

template < typename BlockT, int Flags >
struct BlockMatrixTraits< SparseBlockMatrix< BlockT, Flags > >
{
	typedef BlockT BlockType ;

	enum {
		is_compressed = Flags & BlockMatrixFlags::COMPRESSED,
		is_symmetric  = Flags & BlockMatrixFlags::SYMMETRIC
	} ;

	typedef SparseBlockIndex< BlockT, is_compressed > SparseIndexType ;
	typedef typename BlockMatrixTraits< SparseIndexType >::Index Index ;
} ;

}

#endif // SPARSEBLOCKMATRIX_HH
