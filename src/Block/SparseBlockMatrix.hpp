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
		  m_colMajorValid( false ),
		  m_transposeCached( false )
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

	Index rowsOfBlocks() const { return rowOffsets().size() - 1 ; }
	Index colsOfBlocks() const { return colOffsets().size() - 1 ; }

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
	void cacheTranspose() ;
	bool computeColMajorIndex() ;

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

	void computeColMajorIndex(SparseBlockIndex< BlockType > &cmIndex) ;

	Index blockRows( Index row ) const { return rowOffsets()[ row + 1 ] - rowOffsets()[ row ] ; }
	Index blockCols( Index col ) const { return colOffsets()[ col + 1 ] - colOffsets()[ col ] ; }
	Index rowOffset( Index row ) const { return rowOffsets()[row] ; }
	Index colOffset( Index col ) const { return colOffsets()[col] ; }
	const std::vector< Index >& rowOffsets() const { return m_colMajorIndex.innerOffsets ; }
	const std::vector< Index >& colOffsets() const { return m_rowMajorIndex.innerOffsets ; }

	template < typename VecT >
	typename VecT::SegmentReturnType rowSegment( VecT& v, Index rowBlockIdx ) const
	{
		return v.segment( rowOffset( rowBlockIdx ), blockRows( rowBlockIdx ) ) ;
	}

	template < typename VecT >
	typename VecT::ConstSegmentReturnType rowSegment( const VecT& v, Index rowBlockIdx ) const
	{
		return v.segment( rowOffset( rowBlockIdx ), blockRows( rowBlockIdx ) ) ;
	}

	template < typename VecT >
	typename VecT::SegmentReturnType colSegment( VecT& v, Index colBlockIdx ) const
	{
		return v.segment( colOffset( colBlockIdx ), blockCols( colBlockIdx ) ) ;
	}

	template < typename VecT >
	typename VecT::ConstSegmentReturnType colSegment( const VecT& v, Index colBlockIdx ) const
	{
		return v.segment( colOffset( colBlockIdx ), blockCols( colBlockIdx ) ) ;
	}

	template < typename RhsT, typename ResT >
	void innerMultiply( const SparseIndexType &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const ;

	template < typename RhsT, typename ResT >
	void innerTransposedMultiply( const SparseIndexType &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const ;

	void setInnerOffets( SparseIndexType& index, const std::vector< Index > &blockSizes ) ;

	std::size_t m_nBlocks ;

	bool m_colMajorValid ;
	bool m_transposeCached ;

	SparseIndexType m_rowMajorIndex ;
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
