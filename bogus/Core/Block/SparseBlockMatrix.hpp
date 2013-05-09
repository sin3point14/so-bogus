/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_SPARSEBLOCKMATRIX_HH
#define BOGUS_SPARSEBLOCKMATRIX_HH

#include "BlockMatrix.hpp"
#include "SparseBlockIndex.hpp"
#include "Expressions.hpp"

namespace bogus
{

template < typename BlockT, int Flags > class SparseBlockMatrix  ;
template < bool Symmetric > struct SparseBlockMatrixFinalizer ;
template < typename Derived, bool Major > struct SparseBlockIndexGetter ;

template < typename Derived >
class SparseBlockMatrixBase : public BlockMatrixBase< Derived >
{

public:
	typedef BlockMatrixBase< Derived > Base ;
	typedef BlockMatrixTraits< Derived > Traits ;

	typedef typename Traits::Index Index ;
	typedef typename Traits::SparseIndexType SparseIndexType ;
	typedef typename Traits::RowIndexType RowIndexType ;
	typedef typename Traits::ColIndexType ColIndexType ;
	typedef typename Traits::UncompressedIndexType UncompressedIndexType ;

	typedef typename Traits::BlockType BlockType ;
	typedef typename Traits::BlockPtr BlockPtr ;
	static const BlockPtr InvalidBlockPtr ;

	typedef typename Base::ConstTransposeReturnType  ConstTransposeReturnType ;

	using Base::rows ;
	using Base::cols ;
	using Base::blocks ;
	using Base::derived ;

	SparseBlockMatrixBase() ;

	template < typename OtherDerived >
	Derived& operator= ( const SparseBlockMatrixBase< OtherDerived > &source ) ;

	template < typename OtherDerived >
	Derived& operator= ( const Transpose< SparseBlockMatrixBase< OtherDerived > > &source ) ;

	template < typename LhsT, typename RhsT >
	Derived& operator= ( const Product< LhsT, RhsT > &prod ) ;

	void setRows( const Index nBlocks, const Index* rowsPerBlock ) ;
	void setRows( const std::vector< Index > &rowsPerBlock )
	{
		setRows( rowsPerBlock.size(), &rowsPerBlock[0] ) ;
	}
	void setRows( const Index nBlocks, const Index rowsPerBlock )
	{
		setRows( std::vector< Index >( nBlocks, rowsPerBlock ) ) ;
	}
	void setRows( const Index nBlocks )
	{
		setRows( nBlocks, BlockType::RowsAtCompileTime ) ;
	}

	void setCols( const Index nBlocks, const Index* colsPerBlock ) ;
	void setCols( const std::vector< Index > &colsPerBlock )
	{
		setCols( colsPerBlock.size(), &colsPerBlock[0] ) ;
	}
	void setCols( const Index nBlocks, const Index colsPerBlock )
	{
		setCols( std::vector< Index >( nBlocks, colsPerBlock ) ) ;
	}
	void setCols( const Index nBlocks )
	{
		setCols( nBlocks, BlockType::ColsAtCompileTime ) ;
	}

	Index rowsOfBlocks() const { return rowOffsets().size() - 1 ; }
	Index colsOfBlocks() const { return colOffsets().size() - 1 ; }
	Index blockRows( Index row ) const { return rowOffsets()[ row + 1 ] - rowOffsets()[ row ] ; }
	Index blockCols( Index col ) const { return colOffsets()[ col + 1 ] - colOffsets()[ col ] ; }

	void reserve( std::size_t nBlocks )
	{
		m_blocks.reserve( nBlocks ) ;
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
	
	BlockType& insertBackOuterInner( Index outer, Index inner ) ;

	void finalize() ;
	void clear() ;

	template< typename PrecisionT >	
	void prune( const PrecisionT& precision ) ;

	bool computeMinorIndex() ;

	void cacheTranspose() ;
	bool transposeCached() const { return m_transposeIndex.valid ; }

	std::size_t nBlocks() const { return m_nBlocks ; }

	const BlockType& block( BlockPtr ptr ) const
	{
		return m_blocks[ ptr ] ;
	}

	BlockType& block( BlockPtr ptr )
	{
		return m_blocks[ ptr ] ;
	}

	BlockPtr blockPtr( Index row, Index col ) const ;

	// Warning: block has to exists
	BlockType& diagonal( const Index row ) ;
	const BlockType& diagonal( const Index row ) const ;

	// Warning: block has to exist
	BlockType& block( Index row, Index col )
	{ return block( blockPtr( row, col ) ) ; }
	const BlockType& block( Index row, Index col ) const
	{ return block( blockPtr( row, col ) ) ; }

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

	ConstTransposeReturnType transpose() const { return Transpose< SparseBlockMatrixBase< Derived > >( *this ) ; }

	template < bool Transpose, typename RhsT, typename ResT >
	void multiply( const RhsT& rhs, ResT& res ) const ;

	template < typename RhsT, typename ResT >
	void splitRowMultiply( const Index row, const RhsT& rhs, ResT& res ) const ;

	template< typename OtherDerived >
	void cloneDimensions( const BlockMatrixBase< OtherDerived > &source ) ;

	template < typename BlockT2 >
	void cloneStructure( const SparseBlockMatrix< BlockT2, Traits::flags > &source ) ;

	template < bool ColWise, typename LhsT, typename RhsT >
	void setFromProduct( const Product< LhsT, RhsT > &prod ) ;
	
#ifdef BOGUS_WITH_BOOST_SERIALIZATION
	template < typename Archive >
	void serialize( Archive & ar, const unsigned int file_version ) ;
#endif

protected:

	using Base::m_cols ;
	using Base::m_rows ;
	using Base::m_blocks ;

	typedef SparseBlockMatrixFinalizer< Traits::is_symmetric > Finalizer ;
	friend struct SparseBlockIndexGetter< Derived, true > ;
	friend struct SparseBlockIndexGetter< Derived, false > ;

	void allocateBlock( BlockPtr &ptr ) ;
	void prealloc( std::size_t nBlocks ) ;

	void computeMinorIndex( UncompressedIndexType &cmIndex) const ;

	const UncompressedIndexType& getOrComputeMinorIndex( UncompressedIndexType &tempIndex) const ;

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

	template< typename IndexT >
	void setInnerOffets( IndexT& index, const Index nBlocks, const Index *blockSizes ) const ;

	template < bool ColWise, typename LhsIndex, typename RhsIndex, typename LhsBlock, typename RhsBlock, typename LhsGetter, typename RhsGetter  >
	void setFromProduct( const LhsIndex &lhsIdx, const RhsIndex &rhsIdx,
						 const LhsBlock  *lhsData,
						 const RhsBlock  *rhsData,
						 const LhsGetter &lhsGetter,
						 const RhsGetter &rhsGetter
						  ) ;

	std::size_t m_nBlocks ;

	SparseIndexType m_majorIndex ;
	// Minor index is always uncompressed, as the blocks cannot be contiguous
	// For a symmetric matrix, do not store diagonal block in the minor and transpose index
	UncompressedIndexType m_minorIndex ;
	SparseIndexType m_transposeIndex ;
} ;

template < typename BlockT, int Flags >
struct BlockMatrixTraits< SparseBlockMatrix< BlockT, Flags > > : public BlockMatrixTraits< BlockObjectBase< SparseBlockMatrix< BlockT, Flags > > >
{
	typedef BlockMatrixTraits< BlockObjectBase< SparseBlockMatrix< BlockT, Flags > > > BaseTraits ;
	typedef typename BaseTraits::Index Index ;
	typedef typename BaseTraits::BlockPtr BlockPtr ;

	typedef BlockT BlockType ;

	enum {
		is_compressed  = Flags & flags::COMPRESSED,
		is_symmetric   = Flags & flags::SYMMETRIC,
		is_col_major   = Flags & flags::COL_MAJOR
	} ;
	enum {
		flags         = Flags
	} ;

	typedef SparseBlockIndex< is_compressed, Index, BlockPtr > SparseIndexType ;

	typedef SparseBlockIndex< false, Index, BlockPtr  > UncompressedIndexType ;
	typedef SparseBlockIndex< is_compressed && !is_col_major, Index, BlockPtr > RowIndexType ;
	typedef SparseBlockIndex< is_compressed && is_col_major , Index, BlockPtr > ColIndexType ;

	typedef Transpose< SparseBlockMatrixBase< SparseBlockMatrix< BlockT, Flags > > > ConstTransposeReturnType ;
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
