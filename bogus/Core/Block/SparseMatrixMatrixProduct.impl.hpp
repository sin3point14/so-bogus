/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_SPARSE_MATRIXMATRIX_PRODUCT_IMPL_HPP
#define BOGUS_SPARSE_MATRIXMATRIX_PRODUCT_IMPL_HPP

#include "Expressions.hpp"
#include "Access.hpp"

#include "SparseBlockMatrixBase.hpp"
#include "SparseBlockIndexComputer.hpp"

#include <map>

template < typename LhsT, typename RhsT  >
bogus::Product< LhsT, RhsT > operator* ( const bogus::BlockObjectBase< LhsT >& lhs,
										 const bogus::BlockObjectBase< RhsT > &rhs )
{
	return bogus::Product< LhsT, RhsT >( lhs.derived(), rhs.derived() ) ;
}

namespace bogus
{

template < typename Derived >
template < typename LhsT, typename RhsT >
Derived& SparseBlockMatrixBase<Derived>::operator=( const Product< LhsT, RhsT > &prod )
{
	setFromProduct< true >( prod ) ;
	return derived() ;
}


template < typename Derived >
template < bool ColWise, typename LhsT, typename RhsT >
void SparseBlockMatrixBase<Derived>::setFromProduct( const Product< LhsT, RhsT > &prod )
{
	typedef Product< LhsT, RhsT> Prod ;

	typename Prod::Lhs::EvalType lhs = prod.lhs.object.eval() ;
	typename Prod::Rhs::EvalType rhs = prod.rhs.object.eval() ;
	typedef BlockMatrixTraits< typename Prod::PlainLhsMatrixType > LhsTraits ;
	typedef BlockMatrixTraits< typename Prod::PlainRhsMatrixType > RhsTraits ;

	BOGUS_STATIC_ASSERT( !Prod::transposeLhs || IsTransposable< typename LhsTraits::BlockType >::Value,
						 TRANSPOSE_IS_NOT_DEFINED_FOR_THIS_BLOCK_TYPE
	) ;
	BOGUS_STATIC_ASSERT( !Prod::transposeRhs || IsTransposable< typename RhsTraits::BlockType >::Value,
						 TRANSPOSE_IS_NOT_DEFINED_FOR_THIS_BLOCK_TYPE
	) ;

	typedef BlockTransposeOption< LhsTraits::is_symmetric, Prod::transposeLhs > LhsGetter ;
	LhsGetter lhsGetter ;
	typedef BlockTransposeOption< RhsTraits::is_symmetric, Prod::transposeRhs > RhsGetter ;
	RhsGetter rhsGetter ;

	clear() ;
	if( Prod::transposeLhs )
	{
		m_rows = lhs->cols() ;
		colMajorIndex().innerOffsets = lhs->rowMajorIndex().innerOffsets;
	} else {
		m_rows = lhs->rows() ;
		colMajorIndex().innerOffsets = lhs->colMajorIndex().innerOffsets;
	}
	if( Prod::transposeRhs )
	{
		m_cols = rhs->rows() ;
		rowMajorIndex().innerOffsets = rhs->colMajorIndex().innerOffsets;
	} else {
		m_cols = rhs->cols() ;
		rowMajorIndex().innerOffsets = rhs->rowMajorIndex().innerOffsets;
	}


	SparseBlockIndexComputer< typename Prod::PlainLhsMatrixType, LhsTraits::is_symmetric,
			ColWise, Prod::transposeLhs> lhsIndexComputer ( *lhs ) ;
	SparseBlockIndexComputer< typename Prod::PlainRhsMatrixType, RhsTraits::is_symmetric,
			!ColWise, Prod::transposeRhs> rhsIndexComputer ( *rhs ) ;


	setFromProduct< ColWise >( lhsIndexComputer.get(), rhsIndexComputer.get(),
							   lhs->data(), rhs->data(),
							   lhsGetter, rhsGetter, prod.lhs.scaling * prod.rhs.scaling ) ;

}

template < bool ColWise, typename Index, typename BlockPtr, bool is_symmetric, bool is_col_major >
struct SparseBlockProductIndex
{
	typedef std::vector< std::pair< bool, BlockPtr > > BlockComputationFactor ;
	typedef std::pair< BlockComputationFactor, BlockComputationFactor > BlockComputation ;

	typedef std::vector< BlockComputation > InnerType;
	typedef typename InnerType::const_iterator InnerIterator;

	static const BlockComputation* get( const InnerIterator& iter ) { return &(*iter) ; }

	SparseBlockIndex< true, Index, BlockPtr > compressed ;
	std::vector< InnerType > to_compute ;

	template< typename LhsIndex, typename RhsIndex >
	void compute(
			const Index outerSize,
			const Index innerSize,
			const LhsIndex &lhsIdx,
			const RhsIndex &rhsIdx )
	{
		assert( lhsIdx.innerSize() == rhsIdx.innerSize() ) ;

		to_compute.resize( outerSize ) ;

#ifndef BOGUS_DONT_PARALLELIZE
		SparseBlockIndex< false, Index, BlockPtr > uncompressed ;
		uncompressed.resizeOuter( outerSize ) ;
#else
		compressed.resizeOuter( outerSize ) ;
#endif

		BlockComputation currentBlock ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for private( currentBlock )
#endif
		for( int i = 0 ; i < (int) outerSize ; ++i )
		{
			const Index last = is_symmetric ? i+1 : innerSize ;
			for( Index j = 0 ; j != last ; ++ j )
			{
				const Index lhsOuter = is_col_major ? j : i ;
				const Index rhsOuter = is_col_major ? i : j ;
				typename LhsIndex::InnerIterator lhs_it ( lhsIdx, lhsOuter ) ;
				typename RhsIndex::InnerIterator rhs_it ( rhsIdx, rhsOuter ) ;

				while( lhs_it && rhs_it )
				{
					if( lhs_it.inner() > rhs_it.inner() ) ++rhs_it ;
					else if( lhs_it.inner() < rhs_it.inner() ) ++lhs_it ;
					else {
						currentBlock.first. push_back( std::make_pair( lhs_it.inner() > lhsOuter, lhs_it.ptr() ) ) ;
						currentBlock.second.push_back( std::make_pair( rhs_it.inner() > rhsOuter, rhs_it.ptr() ) ) ;
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

template < typename Index, typename BlockPtr, bool is_symmetric, bool is_col_major >
struct SparseBlockProductIndex< true, Index, BlockPtr, is_symmetric, is_col_major >
{
	typedef std::vector< std::pair< bool, BlockPtr > > BlockComputationFactor ;
	typedef std::pair< BlockComputationFactor, BlockComputationFactor > BlockComputation ;

	typedef std::map< Index, BlockComputation > InnerType;
	typedef typename InnerType::const_iterator InnerIterator;

	static const BlockComputation* get( const InnerIterator& iter ) { return &(iter->second) ; }

	SparseBlockIndex< true, Index, BlockPtr > compressed ;
	std::vector< InnerType > to_compute ;

	template< typename LhsIndex, typename RhsIndex >
	void compute(
			const Index outerSize,
			const Index innerSize,
			const LhsIndex &lhsIdx,
			const RhsIndex &rhsIdx )
	{
		assert( lhsIdx.outerSize() == rhsIdx.outerSize() ) ;
		( void ) innerSize ;

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
			for( int ii = 0 ; ii < (int) productSize ; ++ii )
			{
				const Index i = (Index) ii ;

				if( is_col_major )
				{
					for( typename RhsIndex::InnerIterator rhs_it ( rhsIdx, i ) ; rhs_it ; ++rhs_it )
					{
						for( typename LhsIndex::InnerIterator lhs_it ( lhsIdx, i ) ;
							 lhs_it && ( !is_symmetric || lhs_it.inner() <= rhs_it.inner() ) ;
							 ++lhs_it )
						{
							BlockComputation &bc = loc_compute[ rhs_it.inner() ][ lhs_it.inner() ] ;
							bc.first.push_back( std::make_pair( lhs_it.inner() > i, lhs_it.ptr() ) ) ;
							bc.second.push_back( std::make_pair( rhs_it.inner() > i, rhs_it.ptr() ) ) ;
						}
					}
				} else {
					for( typename LhsIndex::InnerIterator lhs_it ( lhsIdx, i ) ; lhs_it ; ++lhs_it )
					{
						for( typename RhsIndex::InnerIterator rhs_it ( rhsIdx, i ) ;
							 rhs_it && ( !is_symmetric || rhs_it.inner() <= lhs_it.inner() ) ;
							 ++rhs_it )
						{
							BlockComputation &bc = loc_compute[ lhs_it.inner() ][ rhs_it.inner() ] ;
							bc.first.push_back( std::make_pair( lhs_it.inner() > i, lhs_it.ptr() ) ) ;
							bc.second.push_back( std::make_pair( rhs_it.inner() > i, rhs_it.ptr() ) ) ;
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
template < bool ColWise, typename LhsIndex, typename RhsIndex, typename LhsBlock, typename RhsBlock,
		   typename LhsGetter, typename RhsGetter >
void SparseBlockMatrixBase<Derived>::setFromProduct(const LhsIndex &lhsIdx,
		const RhsIndex &rhsIdx,
		const LhsBlock *lhsData,
		const RhsBlock *rhsData,
		const LhsGetter &lhsGetter, const RhsGetter &rhsGetter,
		Scalar scaling)
{
	typedef SparseBlockProductIndex< ColWise, Index, BlockPtr, Traits::is_symmetric, Traits::is_col_major> ProductIndex ;
	typedef typename ProductIndex::BlockComputation BlockComputation ;

	clear() ;

	assert( lhsIdx.valid ) ;
	assert( rhsIdx.valid ) ;
	rowMajorIndex().resizeOuter( colMajorIndex().innerSize() ) ;
	colMajorIndex().resizeOuter( rowMajorIndex().innerSize() ) ;

	const unsigned outerSize = majorIndex().outerSize() ;

	ProductIndex productIndex  ;

	productIndex.compute( majorIndex().outerSize(), minorIndex().outerSize(), lhsIdx, rhsIdx ) ;


	prealloc( productIndex.compressed.outer[ outerSize ] ) ;

	std::vector< const BlockComputation* > flat_compute ( nBlocks() ) ;
	std::vector< Index > outerIndices ( nBlocks() ) ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( int i = 0 ; i < (int) outerSize ; ++i )
	{
		typename ProductIndex::InnerIterator j = productIndex.to_compute[i].begin() ;
		for( typename SparseBlockIndex< true, Index, BlockPtr >::InnerIterator c_it( productIndex.compressed, i ) ;
			 c_it ; ++c_it, ++j )
		{
			outerIndices[ c_it.ptr() ] = i  ;
			flat_compute[ c_it.ptr() ] = ProductIndex::get( j ) ;
		}
	}

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( long i = 0 ; i < (long) nBlocks() ; ++ i )
	{
		BlockType& b = block( i ) ;
		const BlockComputation &bc = *flat_compute[i] ;
		b = lhsGetter.get( lhsData[ bc.first[0].second ], bc.first[0].first)
				* rhsGetter.get( rhsData[ bc.second[0].second ], bc.second[0].first) ;
		for( unsigned j = 1 ; j != bc.first.size() ; ++ j)
		{
			b += lhsGetter.get( lhsData[ bc.first[j].second ], bc.first[j].first)
					* rhsGetter.get( rhsData[ bc.second[j].second ], bc.second[j].first) ;
		}

		if( scaling != 1 ) b *= scaling ;
	}

	productIndex.compressed.valid  = true ;
	m_majorIndex.move( productIndex.compressed );
	assert( m_majorIndex.valid ) ;

	Finalizer::finalize( *this ) ;

}

} //namespace bogus

#endif // SPARSEPRODUCT_IMPL_HPP
