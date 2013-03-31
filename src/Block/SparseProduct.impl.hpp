#ifndef BOGUS_SPARSEPRODUCT_IMPL_HPP
#define BOGUS_SPARSEPRODUCT_IMPL_HPP

#include "Expressions.hpp"
#include "SparseBlockMatrix.hpp"


template < typename LhsT, typename RhsT >
bogus::Product< LhsT, RhsT > operator* ( const bogus::BlockObjectBase< LhsT >& lhs,
			 const bogus::BlockObjectBase< RhsT > &rhs )
{
	return bogus::Product<  LhsT, RhsT >( lhs, rhs ) ;
}

namespace bogus
{

template < typename Derived >
template < typename LhsT, typename RhsT >
Derived& SparseBlockMatrixBase<Derived>::operator=( const Product< LhsT, RhsT > &prod )
{
	setFromProduct( prod, 1. ) ;
	return this->derived() ;
}

template < typename Derived >
template < typename LhsT, typename RhsT >
void SparseBlockMatrixBase<Derived>::setFromProduct( const Product< LhsT, RhsT > &prod, double scale)
{
	typedef Product< LhsT, RhsT> Prod ;
	typedef BlockMatrixTraits< typename Prod::LhsTraits::MatrixType > LhsTraits ;
	typedef BlockMatrixTraits< typename Prod::RhsTraits::MatrixType > RhsTraits ;

	typedef BlockTranspose< LhsTraits::is_symmetric, Prod::transposeLhs > LhsGetter ;
	LhsGetter lhsGetter ;
	typedef BlockTranspose< RhsTraits::is_symmetric, Prod::transposeRhs > RhsGetter ;
	RhsGetter rhsGetter ;

	SparseBlockIndex< > auxIndexLhs, auxIndexRhs ;

	assert( ! LhsTraits::is_symmetric ) ;
	assert( ! RhsTraits::is_symmetric ) ;

	const SparseBlockIndexBase &lhsIdx = getIndex( Prod::transposeLhs, false, auxIndexLhs ) ;
	const SparseBlockIndexBase &rhsIdx = getIndex( Prod::transposeRhs, true, auxIndexRhs ) ;

	clear() ;
	if( Prod::transposeLhs )
	{
		this->m_rows = prod.lhs.cols() ;
		colMajorIndex().innerOffsets = prod.lhs.rowMajorIndex().innerOffsets;
	} else {
		this->m_rows = prod.lhs.rows() ;
		colMajorIndex().innerOffsets = prod.lhs.colMajorIndex().innerOffsets;
	}
	if( Prod::transposeRhs )
	{
		this->m_cols = prod.rhs.rows() ;
		rowMajorIndex().innerOffsets = prod.rhs.colMajorIndex().innerOffsets;
	} else {
		this->m_cols = prod.rhs.cols() ;
		rowMajorIndex().innerOffsets = prod.rhs.rowMajorIndex().innerOffsets;
	}

	if( lhsIdx.isCompressed() )
	{
		if( rhsIdx.isCompressed() )
		{
			setFromProduct( lhsIdx.asCompressed(), rhsIdx.asCompressed(),
							prod.lhs.blocks(), prod.rhs.blocks(),
							lhsGetter, rhsGetter, scale ) ;
		} else {
			setFromProduct( lhsIdx.asCompressed(), rhsIdx.asUncompressed(),
							prod.lhs.blocks(), prod.rhs.blocks(),
							lhsGetter, rhsGetter, scale ) ;
		}
	} else {
		if( rhsIdx.isCompressed() )
		{
			setFromProduct( lhsIdx.asUncompressed(), rhsIdx.asCompressed(),
							prod.lhs.blocks(), prod.rhs.blocks(),
							lhsGetter, rhsGetter, scale ) ;
		} else {
			setFromProduct( lhsIdx.asUncompressed(), rhsIdx.asUncompressed(),
							prod.lhs.blocks(), prod.rhs.blocks(),
							lhsGetter, rhsGetter, scale ) ;
		}
	}

}

template < typename Derived >
template < typename LhsIndex, typename RhsIndex, typename LhsBlock, typename RhsBlock, typename LhsGetter, typename RhsGetter >
void SparseBlockMatrixBase<Derived>::setFromProduct(const LhsIndex &lhsIdx,
		const RhsIndex &rhsIdx,
		const std::vector< LhsBlock > &lhsData,
		const std::vector< RhsBlock > &rhsData, const LhsGetter &lhsGetter, const RhsGetter &rhsGetter,
		double scale)
{
	typedef std::pair< std::vector< BlockPtr >, std::vector< BlockPtr > > BlockComputation ;

	assert( lhsIdx.innerSize() == rhsIdx.innerSize() ) ;
	rowMajorIndex().resizeOuter( colMajorIndex().innerSize() ) ;
	colMajorIndex().resizeOuter( rowMajorIndex().innerSize() ) ;

	std::vector< std::vector< BlockComputation > > to_compute ( majorIndex().outerSize() ) ;

	SparseBlockIndex< true > compressed ;
	compressed.resizeOuter( majorIndex().outerSize() ) ;

	BlockComputation currentBlock ;

	for( Index i = 0 ; i != majorIndex().outerSize() ; ++i )
	{
		const Index last = Traits::is_symmetric ? i+1 : minorIndex().outerSize() ;
		for( Index j = 0 ; j != last ; ++ j )
		{
			typename LhsIndex::InnerIterator lhs_it ( lhsIdx, Traits::is_col_major ? j : i ) ;
			typename RhsIndex::InnerIterator rhs_it ( rhsIdx, Traits::is_col_major ? i : j ) ;


			while( lhs_it && rhs_it )
			{
				if( lhs_it.inner() > rhs_it.inner() ) ++rhs_it ;
				else if( lhs_it.inner() < rhs_it.inner() ) ++lhs_it ;
				else {
					currentBlock.first.push_back( lhs_it.ptr() ) ;
					currentBlock.second.push_back( rhs_it.ptr() ) ;
				}
			}

			if( !currentBlock.first.empty() )
			{
				to_compute[i].push_back( currentBlock ) ;
				currentBlock.first.clear() ;
				currentBlock.second.clear() ;
				compressed.insertBack( i, j, 0 ) ;
			}
		}
	}

	compressed.finalize() ;

	prealloc( compressed.outer[ compressed.outerSize() ] ) ;

	std::vector< const BlockComputation* > flat_compute ( nBlocks() ) ;
	std::vector< Index > outerIndices ( nBlocks() ) ;

	for( Index i = 0 ; i != majorIndex().outerSize() ; ++i )
	{
		unsigned j = 0 ;
		for( typename SparseBlockIndex< true >::InnerIterator c_it( compressed, i ) ;
			 c_it ; ++c_it, ++j )
		{
			outerIndices[ c_it.ptr() ] = i  ;
			flat_compute[ c_it.ptr() ] = &to_compute[i][j] ;
		}
	}

	for( BlockPtr i = 0 ; i < nBlocks() ; ++ i )
	{
		BlockType& b = block( i ) ;
		const BlockComputation &bc = *flat_compute[i] ;
		b.setZero() ;
		for( unsigned j = 0 ; j != bc.first.size() ; ++ j)
		{
			b += lhsGetter.get( lhsData[ bc.first[j] ], false )
					* rhsGetter.get( rhsData[ bc.second[j] ], false ) ;
		}
		b *= scale ;
	}

	m_majorIndex = compressed ;
	assert( m_majorIndex.valid ) ;
}


}

#endif // SPARSEPRODUCT_IMPL_HPP
