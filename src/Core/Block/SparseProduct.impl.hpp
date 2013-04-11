#ifndef BOGUS_SPARSEPRODUCT_IMPL_HPP
#define BOGUS_SPARSEPRODUCT_IMPL_HPP

#include "Expressions.hpp"
#include "SparseBlockMatrix.hpp"

#include <map>

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
	setFromProduct< true >( prod, 1. ) ;
	return this->derived() ;
}

template < typename Derived >
template < bool ColWise, typename LhsT, typename RhsT >
void SparseBlockMatrixBase<Derived>::setFromProduct( const Product< LhsT, RhsT > &prod, double scale )
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

	const SparseBlockIndexBase &lhsIdx = prod.lhs.getIndex( Prod::transposeLhs, ColWise, auxIndexLhs ) ;
	const SparseBlockIndexBase &rhsIdx = prod.rhs.getIndex( Prod::transposeRhs, !ColWise, auxIndexRhs ) ;

	if( lhsIdx.isCompressed() )
	{
		if( rhsIdx.isCompressed() )
		{
			setFromProduct< ColWise >( lhsIdx.asCompressed(), rhsIdx.asCompressed(),
							prod.lhs.blocks(), prod.rhs.blocks(),
							lhsGetter, rhsGetter, scale ) ;
		} else {
			setFromProduct< ColWise >( lhsIdx.asCompressed(), rhsIdx.asUncompressed(),
							prod.lhs.blocks(), prod.rhs.blocks(),
							lhsGetter, rhsGetter, scale ) ;
		}
	} else {
		if( rhsIdx.isCompressed() )
		{
			setFromProduct< ColWise >( lhsIdx.asUncompressed(), rhsIdx.asCompressed(),
							prod.lhs.blocks(), prod.rhs.blocks(),
							lhsGetter, rhsGetter, scale ) ;
		} else {
			setFromProduct< ColWise >( lhsIdx.asUncompressed(), rhsIdx.asUncompressed(),
							prod.lhs.blocks(), prod.rhs.blocks(),
							lhsGetter, rhsGetter, scale ) ;
		}
	}


}

template < bool ColWise, typename Derived >
struct SparseBlockProductIndex
{
	typedef SparseBlockMatrixBase<Derived> SparseMatrixT ;
	typedef typename SparseMatrixT::Traits Traits ;
	typedef typename SparseMatrixT::BlockPtr BlockPtr ;
	typedef typename SparseMatrixT::Index Index ;
	typedef std::pair< std::vector< BlockPtr >, std::vector< BlockPtr > > BlockComputation ;

	typedef std::vector< BlockComputation > InnerType;
	typedef typename InnerType::const_iterator InnerIterator;

	static const BlockComputation* get( const InnerIterator& iter ) { return &(*iter) ; }

	SparseBlockIndex< true > compressed ;
	std::vector< InnerType > to_compute ;

	template< typename LhsIndex, typename RhsIndex >
	void compute(
		const SparseMatrixT &matrix,
		const LhsIndex &lhsIdx,
		const RhsIndex &rhsIdx )
	{
		assert( lhsIdx.innerSize() == rhsIdx.innerSize() ) ;

		const Index outerSize = matrix.majorIndex().outerSize() ;
		const Index innerSize = matrix.minorIndex().outerSize() ;
		to_compute.resize( outerSize ) ;

#ifndef BOGUS_DONT_PARALLELIZE
		SparseBlockIndex< false > uncompressed ;
		uncompressed.resizeOuter( outerSize ) ;
#else
		compressed.resizeOuter( outerSize ) ;
#endif

		BlockComputation currentBlock ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for private( currentBlock )
#endif
		for( Index i = 0 ; i < outerSize ; ++i )
		{
			const Index last = Traits::is_symmetric ? i+1 : innerSize ;
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
						++lhs_it ;
						++rhs_it ;
					}
				}

				if( !currentBlock.first.empty() )
				{
					to_compute[i].push_back( currentBlock ) ;
					currentBlock.first.clear() ;
					currentBlock.second.clear() ;
#ifndef BOGUS_DONT_PARALLELIZE
					uncompressed.insertBack( i, j, 0 ) ;
#else
					compressed.insertBack( i, j, 0 ) ;
#endif
				}
			}
		}

#ifndef BOGUS_DONT_PARALLELIZE
		compressed = uncompressed ;
#else
		compressed.finalize() ;
#endif

	}

} ;

template < typename Derived >
struct SparseBlockProductIndex< true, Derived >
{
	typedef SparseBlockMatrixBase<Derived> SparseMatrixT ;
	typedef typename SparseMatrixT::Traits Traits ;
	typedef typename SparseMatrixT::BlockPtr BlockPtr ;
	typedef typename SparseMatrixT::Index Index ;
	typedef std::pair< std::vector< BlockPtr >, std::vector< BlockPtr > > BlockComputation ;

	typedef std::map< Index, BlockComputation > InnerType;
	typedef typename InnerType::const_iterator InnerIterator;

	static const BlockComputation* get( const InnerIterator& iter ) { return &(iter->second) ; }

	SparseBlockIndex< true > compressed ;
	std::vector< InnerType > to_compute ;

	template< typename LhsIndex, typename RhsIndex >
	void compute(
		const SparseMatrixT &matrix,
		const LhsIndex &lhsIdx,
		const RhsIndex &rhsIdx )
	{
		assert( lhsIdx.outerSize() == rhsIdx.outerSize() ) ;

		const Index outerSize = matrix.majorIndex().outerSize() ;
		const Index productSize = lhsIdx.outerSize() ;

		to_compute.resize( outerSize ) ;
		compressed.resizeOuter( outerSize ) ;

#ifdef BOGUS_DONT_PARALLELIZE
		std::vector< InnerType > &loc_compute = to_compute ;
#else
#pragma omp parallel
		{
		std::vector< InnerType > loc_compute( outerSize ) ;
#pragma omp for
#endif
		for( Index i = 0 ; i < productSize ; ++i )
		{
			if( Traits::is_col_major )
			{
				for( typename RhsIndex::InnerIterator rhs_it ( rhsIdx, i ) ; rhs_it ; ++rhs_it )
				{
					for( typename LhsIndex::InnerIterator lhs_it ( lhsIdx, i ) ;
						 lhs_it && ( !Traits::is_symmetric || lhs_it.inner() <= rhs_it.inner() ) ;
						 ++lhs_it )
					{
						BlockComputation &bc = loc_compute[ rhs_it.inner() ][ lhs_it.inner() ] ;
						bc.first.push_back( lhs_it.ptr() ) ;
						bc.second.push_back( rhs_it.ptr() ) ;
					}
				}
			} else {
				for( typename LhsIndex::InnerIterator lhs_it ( lhsIdx, i ) ; lhs_it ; ++lhs_it )
				{
					for( typename RhsIndex::InnerIterator rhs_it ( rhsIdx, i ) ;
						 rhs_it && ( !Traits::is_symmetric || rhs_it.inner() <= lhs_it.inner() ) ;
						 ++rhs_it )
					{
						BlockComputation &bc = loc_compute[ lhs_it.inner() ][ rhs_it.inner() ] ;
						bc.first.push_back( lhs_it.ptr() ) ;
						bc.second.push_back( rhs_it.ptr() ) ;
					}
				}
			}
		}

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp critical
		{
		for( Index i = 0 ; i < outerSize ; ++i )
		{
			for( InnerIterator j = loc_compute[i].begin() ; j != loc_compute[i].end() ; ++j )
			{
				const BlockComputation &src_bc = j->second ;
				BlockComputation &dest_bc = to_compute[ i ][ j->first ] ;
				dest_bc.first.insert( dest_bc.first.end(), src_bc.first.begin(), src_bc.first.end() ) ;
				dest_bc.second.insert( dest_bc.second.end(), src_bc.second.begin(), src_bc.second.end() ) ;
			}
		}
		}
		}
#endif

		for( Index i = 0 ; i < outerSize ; ++i )
		{
			for( InnerIterator j = to_compute[i].begin() ; j != to_compute[i].end() ; ++j )
			{
				compressed.insertBack( i, j->first, 0 );
			}
		}

		compressed.finalize() ;

	}

} ;

template < typename Derived >
template < bool ColWise, typename LhsIndex, typename RhsIndex, typename LhsBlock, typename RhsBlock, typename LhsGetter, typename RhsGetter >
void SparseBlockMatrixBase<Derived>::setFromProduct(
		const LhsIndex &lhsIdx,
		const RhsIndex &rhsIdx,
		const std::vector< LhsBlock > &lhsData,
		const std::vector< RhsBlock > &rhsData,
		const LhsGetter &lhsGetter, const RhsGetter &rhsGetter,
		double scale)
{
	typedef SparseBlockProductIndex< ColWise, Derived > ProductIndex ;
	typedef typename ProductIndex::BlockComputation BlockComputation ;

	assert( lhsIdx.valid ) ;
	assert( rhsIdx.valid ) ;
	rowMajorIndex().resizeOuter( colMajorIndex().innerSize() ) ;
	colMajorIndex().resizeOuter( rowMajorIndex().innerSize() ) ;

	const unsigned outerSize = majorIndex().outerSize() ;

	ProductIndex productIndex  ;

	productIndex.compute( *this, lhsIdx, rhsIdx ) ;


	prealloc( productIndex.compressed.outer[ outerSize ] ) ;

	std::vector< const BlockComputation* > flat_compute ( nBlocks() ) ;
	std::vector< Index > outerIndices ( nBlocks() ) ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( Index i = 0 ; i < outerSize ; ++i )
	{
		typename ProductIndex::InnerIterator j = productIndex.to_compute[i].begin() ;
		for( typename SparseBlockIndex< true >::InnerIterator c_it( productIndex.compressed, i ) ;
			 c_it ; ++c_it, ++j )
		{
			outerIndices[ c_it.ptr() ] = i  ;
			flat_compute[ c_it.ptr() ] = ProductIndex::get( j ) ;
		}
	}

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( BlockPtr i = 0 ; i < nBlocks() ; ++ i )
	{
		BlockType& b = block( i ) ;
		const BlockComputation &bc = *flat_compute[i] ;
		b = lhsGetter.get( lhsData[ bc.first[0] ], false )
				* rhsGetter.get( rhsData[ bc.second[0] ], false ) ;
		for( unsigned j = 1 ; j != bc.first.size() ; ++ j)
		{
			b += lhsGetter.get( lhsData[ bc.first[j] ], false )
					* rhsGetter.get( rhsData[ bc.second[j] ], false ) ;
		}
		b *= scale ;
	}

	productIndex.compressed.valid  = true ;
	m_majorIndex = productIndex.compressed ;
	assert( m_majorIndex.valid ) ;

	m_minorIndex.valid = false ;
}

} //namespace bogus

#endif // SPARSEPRODUCT_IMPL_HPP
