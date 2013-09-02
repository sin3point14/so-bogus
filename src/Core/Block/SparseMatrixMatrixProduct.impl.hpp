/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_SPARSE_MATRIXMATRIX_PRODUCT_IMPL_HPP
#define BOGUS_SPARSE_MATRIXMATRIX_PRODUCT_IMPL_HPP

#include "Expressions.hpp"
#include "Access.hpp"

#include "SparseBlockMatrixBase.hpp"
#include "SparseBlockIndexComputer.hpp"

#include <map>

namespace bogus
{

template < typename Derived >
template < typename LhsT, typename RhsT >
Derived& SparseBlockMatrixBase<Derived>::operator=( const Product< LhsT, RhsT > &prod )
{
	setFromProduct< true >( prod ) ;
	return derived() ;
}


namespace mm_impl
{

template < bool ColWise, typename Index, typename BlockPtr, bool is_symmetric, bool is_col_major >
struct SparseBlockProductIndex
{
	typedef std::vector< std::pair< bool, BlockPtr > > BlockComputationFactor ;
	typedef std::pair< BlockComputationFactor, BlockComputationFactor > BlockComputation ;

	typedef std::vector< BlockComputation > InnerType;
	typedef typename InnerType::const_iterator InnerIterator;

	static const BlockComputation* get( const InnerIterator& iter ) { return &(*iter) ; }

	typedef SparseBlockIndex< true, Index, BlockPtr > CompressedIndexType ;
	CompressedIndexType compressed ;
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

	typedef SparseBlockIndex< true, Index, BlockPtr > CompressedIndexType ;
	CompressedIndexType compressed ;
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

template < bool LhsRuntimeTest, bool RhsRunTimeTest,
           bool LhsCompileTimeTranspose, bool RhsCompileTimeTranspose >
struct BinaryTransposeOption {

	template< typename LhsT, typename RhsT, typename ResT >
	static inline void mm_assign (
	        const LhsT &lhs, const RhsT &rhs, ResT &res,
	        bool transposeLhs, bool transposeRhs )
	{
		if( transposeLhs )
			if( transposeRhs )
				res = transpose_block( lhs ) * transpose_block( rhs ) ;
			else
				res = transpose_block( lhs ) * rhs ;
		else
			if( transposeRhs )
				res = lhs * transpose_block( rhs ) ;
			else
				res = lhs * rhs ;
	}

	template< typename LhsT, typename RhsT, typename ResT >
	static inline void mm_add (
	        const LhsT &lhs, const RhsT &rhs, ResT &res,
	        bool transposeLhs, bool transposeRhs )
	{
		if( transposeLhs )
			if( transposeRhs )
				res += transpose_block( lhs ) * transpose_block( rhs ) ;
			else
				res += transpose_block( lhs ) * rhs ;
		else
			if( transposeRhs )
				res += lhs * transpose_block( rhs ) ;
			else
				res += lhs * rhs ;
	}
} ;

// Lhs transpose known at compile time
template < bool LhsTranspose, bool RhsTranspose >
struct BinaryTransposeOption< false, true, LhsTranspose, RhsTranspose >
{
	template< typename LhsT, typename RhsT, typename ResT >
	static inline void mm_assign (
	        const LhsT &lhs, const RhsT &rhs, ResT &res,
	        bool, bool transposeRhs )
	{
		if( transposeRhs )
			res = TransposeIf< LhsTranspose >::get( lhs ) * transpose_block( rhs ) ;
		else
			res = TransposeIf< LhsTranspose >::get( lhs ) * rhs ;
	}
	template< typename LhsT, typename RhsT, typename ResT >
	static inline void mm_add (
	        const LhsT &lhs, const RhsT &rhs, ResT &res,
	        bool, bool transposeRhs )
	{
		if( transposeRhs )
			res += TransposeIf< LhsTranspose >::get( lhs ) * transpose_block( rhs ) ;
		else
			res += TransposeIf< LhsTranspose >::get( lhs ) * rhs ;
	}

} ;

// Rhs transpose known at compile time
template < bool LhsTranspose, bool RhsTranspose >
struct BinaryTransposeOption< true, false, LhsTranspose, RhsTranspose >
{
	template< typename LhsT, typename RhsT, typename ResT >
	static inline void mm_assign (
	        const LhsT &lhs, const RhsT &rhs, ResT &res,
	        bool transposeLhs, bool )
	{
		if( transposeLhs )
			res = transpose_block( lhs ) * TransposeIf< RhsTranspose >::get( rhs )  ;
		else
			res =                    lhs * TransposeIf< RhsTranspose >::get( rhs )  ;
	}
	template< typename LhsT, typename RhsT, typename ResT >
	static inline void mm_add (
	        const LhsT &lhs, const RhsT &rhs, ResT &res,
	        bool transposeLhs, bool )
	{
		if( transposeLhs )
			res += transpose_block( lhs ) * TransposeIf< RhsTranspose >::get( rhs )  ;
		else
			res +=                    lhs * TransposeIf< RhsTranspose >::get( rhs )  ;
	}

} ;

// Both transpose known at compile time
template < bool LhsTranspose, bool RhsTranspose >
struct BinaryTransposeOption< false, false, LhsTranspose, RhsTranspose >
{
	template< typename LhsT, typename RhsT, typename ResT >
	static inline void mm_assign (
	        const LhsT &lhs, const RhsT &rhs, ResT &res,
	        bool , bool )
	{
		res = TransposeIf< LhsTranspose >::get( lhs )
		        * TransposeIf< RhsTranspose >::get( rhs )  ;
	}

	template< typename LhsT, typename RhsT, typename ResT >
	static inline void mm_add (
	        const LhsT &lhs, const RhsT &rhs, ResT &res,
	        bool , bool )
	{
		res += TransposeIf< LhsTranspose >::get( lhs )
		        * TransposeIf< RhsTranspose >::get( rhs )  ;
	}


} ;

template< typename TransposeOption, typename ProductIndex,
          typename LhsBlock, typename RhsBlock, typename ResBlock >
static void compute_blocks(
        const ProductIndex &productIndex, const std::size_t nBlocks,
        const LhsBlock *lhsData, const RhsBlock *rhsData,
        ResBlock* resData, typename BlockTraits< ResBlock >::Scalar scaling)
{

	typedef typename ProductIndex::BlockComputation BlockComputation ;

	std::vector< const BlockComputation* > flat_compute ( nBlocks ) ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( int i = 0 ; i < (int) productIndex.to_compute.size() ; ++i )
	{
		typename ProductIndex::InnerIterator j = productIndex.to_compute[i].begin() ;
		for( typename ProductIndex::CompressedIndexType::InnerIterator c_it( productIndex.compressed, i ) ;
		     c_it ; ++c_it, ++j )
		{
			flat_compute[ c_it.ptr() ] = ProductIndex::get( j ) ;
		}
	}

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( long i = 0 ; i < (long) nBlocks ; ++ i )
	{
		ResBlock& b = resData[ i ] ;
		const BlockComputation &bc = *flat_compute[i] ;
		TransposeOption::mm_assign(
		            lhsData[ bc.first[0].second ], rhsData[ bc.second[0].second ],
		        b, bc.first[0].first, bc.second[0].first ) ;

		for( unsigned j = 1 ; j != bc.first.size() ; ++ j)
		{
			TransposeOption::mm_add(
			            lhsData[ bc.first[j].second ], rhsData[ bc.second[j].second ],
			        b, bc.first[j].first, bc.second[j].first ) ;
		}

		if( scaling != 1 ) b *= scaling ;
	}

}



} //namespace mm_impl

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

	rowMajorIndex().resizeOuter( colMajorIndex().innerSize() ) ;
	colMajorIndex().resizeOuter( rowMajorIndex().innerSize() ) ;

	{

		typedef mm_impl::SparseBlockProductIndex< ColWise, Index, BlockPtr,
		        Traits::is_symmetric, Traits::is_col_major> ProductIndex ;
		ProductIndex productIndex  ;

		{
			SparseBlockIndexComputer< typename Prod::PlainLhsMatrixType, LhsTraits::is_symmetric,
					ColWise, Prod::transposeLhs> lhsIndexComputer ( *lhs ) ;
			SparseBlockIndexComputer< typename Prod::PlainRhsMatrixType, RhsTraits::is_symmetric,
					!ColWise, Prod::transposeRhs> rhsIndexComputer ( *rhs ) ;

			productIndex.compute( majorIndex().outerSize(), minorIndex().outerSize(),
								  lhsIndexComputer.get(), rhsIndexComputer.get() ) ;
		}

		const unsigned outerSize = majorIndex().outerSize() ;
		prealloc( productIndex.compressed.outer[ outerSize ] ) ;

		typedef mm_impl::BinaryTransposeOption
				< LhsTraits::is_symmetric && !( BlockTraits< typename LhsTraits::BlockType >::is_self_transpose ),
				RhsTraits::is_symmetric && !( BlockTraits< typename RhsTraits::BlockType >::is_self_transpose ),
				Prod::transposeLhs, Prod::transposeRhs > TransposeOption ;

		mm_impl::template compute_blocks< TransposeOption >( productIndex, nBlocks(),
													lhs->data(), rhs->data(), this->data(),
													prod.lhs.scaling * prod.rhs.scaling ) ;

		productIndex.compressed.valid  = true ;
		m_majorIndex.move( productIndex.compressed );
	}

	assert( m_majorIndex.valid ) ;

	Finalizer::finalize( *this ) ;

}

} //namespace bogus

#endif // SPARSEPRODUCT_IMPL_HPP
