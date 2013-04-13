#ifndef BOGUS_SPARSEBLOCKMATRIX_HH
#define BOGUS_SPARSEBLOCKMATRIX_HH

#include "BlockMatrix.hpp"
#include "SparseBlockIndex.hpp"
#include "Expressions.hpp"

namespace bogus
{

template < typename BlockT, int Flags >
class SparseBlockMatrix  ;

template < bool Symmetric > struct SparseBlockMatrixFinalizer ;
template < typename Derived, bool Major > struct SparseBlockIndexGetter ;

template < typename Derived >
class SparseBlockMatrixBase : public BlockMatrixBase< Derived >
{

public:
	typedef BlockMatrixTraits< Derived > Traits ;

	typedef typename Traits::SparseIndexType SparseIndexType ;
	typedef typename Traits::RowIndexType RowIndexType ;
	typedef typename Traits::ColIndexType ColIndexType ;
	typedef typename Traits::UncompressedIndexType UncompressedIndexType ;

	typedef typename Traits::BlockType BlockType ;
	typedef typename SparseIndexType::BlockPtr BlockPtr ;

	SparseBlockMatrixBase() ;

	template < typename OtherDerived >
	Derived& operator= ( const SparseBlockMatrixBase< OtherDerived > &source ) ;

	template < typename OtherDerived >
	Derived& operator= ( const Transpose< SparseBlockMatrixBase< OtherDerived > > &source ) ;

	template < typename LhsT, typename RhsT >
	Derived& operator= ( const Product< LhsT, RhsT > &prod ) ;

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

	BlockType& insertBackOuterInner( Index outer, Index inner )
	{
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp atomic
#endif
		++m_nBlocks ;

		BlockPtr ptr ;
		allocateBlock( ptr ) ;

		m_majorIndex.insertBack( outer, inner, ptr ) ;
		m_minorIndex.valid = false ;

		return block(ptr) ;
	}
	BlockType& insertBack( Index row, Index col )
	{
		if( Traits::is_col_major )
			return insertBackOuterInner( col, row ) ;
		else
			return insertBackOuterInner( row, col ) ;
	}
	BlockType& insertBackAndResize( Index row, Index col )
	{
		BlockType& block = insertBack( row, col ) ;
		block.resize( blockRows( row ), blockCols( col ) ) ;
		return block ;
	}

	void finalize() ;
	void clear() ;

	bool computeMinorIndex() ;

	void cacheTranspose() ;
	bool transposeCached() const { return m_transposeIndex.valid ; }

	std::size_t nBlocks() const { return m_nBlocks ; }

	const BlockType& block( BlockPtr ptr ) const
	{
		return this->m_blocks[ ptr ] ;
	}

	BlockType& block( BlockPtr ptr )
	{
		return this->m_blocks[ ptr ] ;
	}

	// Warning: block has to exists
	const BlockType& diagonal( const Index row ) const ;
	// Warning: inefficient ; block has to exist
	const BlockType& block( Index row, Index col ) const ;

	const SparseIndexType& majorIndex() const
	{
		return m_majorIndex ;
	}
	const UncompressedIndexType& minorIndex() const
	{
		return m_minorIndex ;
	}
	const SparseIndexType& transposeIndex() const
	{
		return m_transposeIndex ;
	}
	const ColIndexType& colMajorIndex() const ;
	const RowIndexType& rowMajorIndex() const ;

	Transpose< SparseBlockMatrixBase > transpose() const { return Transpose< SparseBlockMatrixBase< Derived > >( *this ) ; }

	template < typename RhsT, typename ResT >
	void multiply( const RhsT& rhs, ResT& res, bool transposed = false ) const ;

	template < typename RhsT, typename ResT >
	void splitRowMultiply( const Index row, const RhsT& rhs, ResT& res ) const ;

	template < typename BlockT2 >
	void cloneStructure( const SparseBlockMatrix< BlockT2, Traits::flags > &source ) ;

	template < bool ColWise, typename LhsT, typename RhsT >
	void setFromProduct( const Product< LhsT, RhsT > &prod , double scale = 1. ) ;

protected:

	typedef SparseBlockMatrixFinalizer< Traits::is_symmetric > Finalizer ;
	friend struct SparseBlockIndexGetter< Derived, true > ;
	friend struct SparseBlockIndexGetter< Derived, false > ;

	void allocateBlock( BlockPtr &ptr )
	{
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp critical
#endif
		{
			ptr = this->m_blocks.size() ;
			this->m_blocks.push_back( BlockType() ) ;
		}
	}

	void computeMinorIndex( UncompressedIndexType &cmIndex) const ;

	const UncompressedIndexType& getOrComputeMinorIndex( UncompressedIndexType &tempIndex) const ;

	Index blockRows( Index row ) const { return rowOffsets()[ row + 1 ] - rowOffsets()[ row ] ; }
	Index blockCols( Index col ) const { return colOffsets()[ col + 1 ] - colOffsets()[ col ] ; }
	Index rowOffset( Index row ) const { return rowOffsets()[row] ; }
	Index colOffset( Index col ) const { return colOffsets()[col] ; }
	const std::vector< Index >& rowOffsets() const { return colMajorIndex().innerOffsets ; }
	const std::vector< Index >& colOffsets() const { return rowMajorIndex().innerOffsets ; }

	ColIndexType& colMajorIndex() ;
	RowIndexType& rowMajorIndex() ;

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

	template < typename IndexT, typename RhsT, typename ResT >
	void innerRowMultiply( const IndexT &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const ;
	template < typename IndexT, typename RhsT, typename ResT >
	void innerColMultiply( const IndexT &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const ;
	template < typename IndexT, typename RhsT, typename ResT >
	void innerRowTransposedMultiply( const IndexT &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const ;
	template < typename IndexT, typename RhsT, typename ResT >
	void innerColTransposedMultiply( const IndexT &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const ;

	template< typename IndexT >
	void setInnerOffets( IndexT& index, const std::vector< Index > &blockSizes ) const ;

	template < bool ColWise, typename LhsIndex, typename RhsIndex, typename LhsBlock, typename RhsBlock, typename LhsGetter, typename RhsGetter  >
	void setFromProduct( const LhsIndex &lhsIdx, const RhsIndex &rhsIdx,
						 const std::vector< LhsBlock > &lhsData,
						 const std::vector< RhsBlock > &rhsData,
						 const LhsGetter &lhsGetter,
						 const RhsGetter &rhsGetter,
						 double scale = 1.
						  ) ;

	template < typename RhsT, typename ResT, typename LocalResT >
	void multiplyAndReduct( const RhsT& rhs, ResT& res, bool transposed, const LocalResT& ) const ;

	std::size_t m_nBlocks ;

	SparseIndexType m_majorIndex ;
	// Minor index is always uncompressed, as the blocks cannot be contiguous
	// For a symmetric matrix, do not store diagonal block in the minor and transpose index
	UncompressedIndexType m_minorIndex ;
	SparseIndexType m_transposeIndex ;
} ;

template < typename BlockT, int Flags >
struct BlockMatrixTraits< SparseBlockMatrix< BlockT, Flags > >
{
	typedef BlockT BlockType ;

	enum {
		is_compressed  = Flags & flags::COMPRESSED,
		is_symmetric   = Flags & flags::SYMMETRIC,
		is_col_major   = Flags & flags::COL_MAJOR
	} ;
	enum {
		flags         = Flags
	} ;

	typedef SparseBlockIndex< is_compressed > SparseIndexType ;
	typedef typename SparseIndexType::Index Index ;

	typedef SparseBlockMatrix< BlockT, Flags > PlainObjectType ;

	typedef SparseBlockIndex< false > UncompressedIndexType ;
	typedef SparseBlockIndex< is_compressed && !is_col_major > RowIndexType ;
	typedef SparseBlockIndex< is_compressed && is_col_major > ColIndexType ;
} ;

template < typename BlockT, int Flags >
class SparseBlockMatrix : public  SparseBlockMatrixBase< SparseBlockMatrix< BlockT, Flags > >
{
	typedef SparseBlockMatrixBase< SparseBlockMatrix< BlockT, Flags > > Base ;
public:
	SparseBlockMatrix() : Base() {}

	template < typename RhsT >
	SparseBlockMatrix( const RhsT& rhs ) : Base()
	{
		Base::operator= ( rhs ) ;
	}

	template < typename RhsT >
	SparseBlockMatrix< BlockT, Flags>& operator=( const RhsT& rhs )
	{
		return Base::operator= ( rhs ) ;
	}
} ;

}

#endif // SPARSEBLOCKMATRIX_HH
