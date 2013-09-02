/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_SPARSE_BLOCK_MATRIX_BASE_IMPL_HPP
#define BOGUS_SPARSE_BLOCK_MATRIX_BASE_IMPL_HPP

#include "Access.hpp"

#include "SparseBlockMatrixBase.hpp"
#include "SparseBlockIndexComputer.hpp"

#include <algorithm>

namespace bogus {

// Finalizer

template < bool Symmetric >
struct SparseBlockMatrixFinalizer
{
	template <typename T>
	static void finalize( const T& ) { }
} ;
template < >
struct SparseBlockMatrixFinalizer<  true >
{
	template < typename Derived >
	static void finalize( SparseBlockMatrixBase< Derived >& matrix )
	{ matrix.computeMinorIndex() ; }
} ;

// Sparse Block Matrix
template < typename Derived >
const typename SparseBlockMatrixBase< Derived >::BlockPtr SparseBlockMatrixBase< Derived >::InvalidBlockPtr( -1 );

template < typename Derived >
typename SparseBlockMatrixBase< Derived >::RowIndexType &SparseBlockMatrixBase< Derived >::rowMajorIndex( )
{
	return SparseBlockIndexGetter< Derived, !Traits::is_col_major >::get( *this ) ;
}
template < typename Derived >
typename SparseBlockMatrixBase< Derived >::ColIndexType &SparseBlockMatrixBase< Derived >::colMajorIndex( )
{
	return SparseBlockIndexGetter< Derived, Traits::is_col_major >::get( *this ) ;
}
template < typename Derived >
const typename SparseBlockMatrixBase< Derived >::RowIndexType &SparseBlockMatrixBase< Derived >::rowMajorIndex( ) const
{
	return SparseBlockIndexGetter< Derived, !Traits::is_col_major >::get( *this ) ;
}
template < typename Derived >
const typename SparseBlockMatrixBase< Derived >::ColIndexType &SparseBlockMatrixBase< Derived >::colMajorIndex( ) const
{
	return SparseBlockIndexGetter< Derived, Traits::is_col_major >::get( *this ) ;
}

template < typename Derived >
SparseBlockMatrixBase< Derived >::SparseBlockMatrixBase()
	: Base()
{
	//Resize to zero
	setRows( 0, (const unsigned*) 0 ) ;
	setCols( 0, (const unsigned*) 0 ) ;
	m_transposeIndex.resizeOuter(0) ;
	m_transposeIndex.valid = false ;
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::setRows(
		const Index nBlocks,
		const unsigned* rowsPerBlock )
{
	setInnerOffets( colMajorIndex(), nBlocks, rowsPerBlock );
	m_rows = colMajorIndex().innerOffsets.back() ;
	rowMajorIndex().resizeOuter( nBlocks ) ;

	if( Traits::is_symmetric ) setCols( nBlocks, rowsPerBlock ) ;
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::setCols(
		const Index nBlocks,
		const unsigned* colsPerBlock )
{
	setInnerOffets( rowMajorIndex(), nBlocks, colsPerBlock );
	m_cols = rowMajorIndex().innerOffsets.back() ;

	colMajorIndex().resizeOuter( nBlocks ) ;
}

template < typename Derived >
template < typename IndexT >
void SparseBlockMatrixBase< Derived >::setInnerOffets(
		IndexT &index, const Index nBlocks, const unsigned *blockSizes) const
{
	index.innerOffsets.resize( nBlocks + 1 ) ;
	index.innerOffsets[0] = 0 ;
	for ( Index i = 0 ; i < nBlocks ; ++ i )
	{
		index.innerOffsets[ i+1 ] = index.innerOffsets[ i ] + blockSizes[i] ;
	}
}

template < typename Derived >
typename SparseBlockMatrixBase< Derived >::BlockType& SparseBlockMatrixBase< Derived >::insertBackOuterInner( Index outer, Index inner )
{
	BlockPtr ptr ;
	BlockType &b = allocateBlock( ptr ) ;

	m_majorIndex.insertBack( outer, inner, ptr ) ;
	m_minorIndex.valid = false ;

	return b ;
}

template < typename Derived >
typename SparseBlockMatrixBase< Derived >::BlockType&
SparseBlockMatrixBase< Derived >::allocateBlock( BlockPtr &ptr )
{
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp critical (allocateBlock)
#endif
	{
		ptr = m_blocks.size() ;
		m_blocks.resize( ptr+1 ) ;
	}
	return m_blocks[ptr] ;
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::prealloc ( std::size_t nBlocks )
{
	clear() ;
	m_blocks.resize( nBlocks ) ;
	m_minorIndex.valid = false ;
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::finalize()
{
	assert( m_majorIndex.valid ) ;
	m_majorIndex.finalize() ;
	m_minorIndex.valid = false ;

	Finalizer::finalize( *this ) ;
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::clear()
{
	m_majorIndex.clear() ;
	m_minorIndex.clear() ;
	m_transposeIndex.clear() ;
	m_transposeIndex.valid = false ;
	m_transposeBlocks.clear()  ;
	m_blocks.clear() ;
}

template < typename Derived >
Derived& SparseBlockMatrixBase< Derived >::prune( const Scalar precision )
{
	MajorIndexType oldIndex = m_majorIndex ;

	typename Traits::BlocksArrayType old_blocks ;
	old_blocks.swap( m_blocks ) ;

	reserve( old_blocks.size() ) ;
	clear() ;

	for( Index outer = 0 ; outer < oldIndex.outerSize() ; ++outer )
	{
		for( typename MajorIndexType::InnerIterator it( oldIndex, outer ) ; it ; ++it )
		{
			if( ! is_zero( old_blocks[ it.ptr() ], precision ) )
			{
				m_majorIndex.insertBack( outer, it.inner(), m_blocks.size() ) ;
				m_blocks.push_back( old_blocks[ it.ptr() ] ) ;
			}
		}
	}

	m_minorIndex.valid = empty() ;
	finalize() ;

	return derived() ;
}

template < typename Derived >
bool SparseBlockMatrixBase< Derived >::computeMinorIndex()
{
	if ( m_minorIndex.valid ) return true ;

	computeMinorIndex( m_minorIndex ) ;

	return m_minorIndex.valid ;
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::computeMinorIndex( UncompressedIndexType &cmIndex) const
{
	cmIndex.clear() ;
	cmIndex.innerOffsets = m_minorIndex.innerOffsets ;

	cmIndex.template setToTranspose< Traits::is_symmetric >( m_majorIndex ) ;
}

template < typename Derived >
const typename SparseBlockMatrixBase< Derived >::UncompressedIndexType&
SparseBlockMatrixBase< Derived >::getOrComputeMinorIndex( UncompressedIndexType &cmIndex) const
{
	if( m_minorIndex.valid ) return m_minorIndex ;
	computeMinorIndex( cmIndex ) ;
	return cmIndex ;
}

template < typename Derived >
void SparseBlockMatrixBase<Derived>::cacheTranspose()
{
	BOGUS_STATIC_ASSERT( IsTransposable< typename Derived::BlockType >::Value,
						 TRANSPOSE_IS_NOT_DEFINED_FOR_THIS_BLOCK_TYPE
	) ;

	if ( m_transposeIndex.valid ) return ;

	computeMinorIndex() ;

	BlockPtr base = 0 ;

	std::vector< BlockPtr > ptrOffsets( m_minorIndex.outerSize() ) ;
	for( Index i = 0 ; i < m_minorIndex.outerSize() ; ++i )
	{
		ptrOffsets[i] = base ;
		base += m_minorIndex.size( i ) ;
	}

	m_transposeBlocks.resize( base ) ;

	assert( m_minorIndex.valid ) ;
	m_transposeIndex = m_minorIndex ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for ( Index i = 0 ; i < m_minorIndex.outerSize() ; ++ i )
	{

		typename MinorIndexType::InnerIterator uncompressed_it
			 ( m_minorIndex , i ) ;
		for( typename TransposeIndexType::InnerIterator it( m_transposeIndex, i ) ;
			 it ; ++ it, ++uncompressed_it )
		{
			const BlockPtr ptr = ptrOffsets[i]++ ;
			m_transposeBlocks[ ptr ] = transpose_block( m_blocks[ uncompressed_it.ptr() ] ) ;
		}
	}

	m_transposeIndex.valid = true ;

}

template < typename Derived >
typename SparseBlockMatrixBase< Derived >::BlockPtr SparseBlockMatrixBase< Derived >::diagonalBlockPtr( const Index row ) const
{
	if( Traits::is_symmetric ) {
		const typename MajorIndexType::InnerIterator it = m_majorIndex.last( row ) ;
		return ( it && it.inner() == row ) ? it.ptr() : InvalidBlockPtr ;
	}

	return blockPtr( row, row ) ;
}

template < typename Derived >
typename SparseBlockMatrixBase< Derived >::BlockPtr SparseBlockMatrixBase< Derived >::blockPtr( Index row, Index col ) const
{
	if( Traits::is_col_major ) std::swap( row, col ) ;

	const typename MajorIndexType::InnerIterator innerIt( majorIndex(), row ) ;

	const typename MajorIndexType::InnerIterator
			found( std::lower_bound( innerIt, innerIt.end(), col ) ) ;

	return found && found.inner() == col ? found.ptr() : InvalidBlockPtr ;
}

template< typename Derived >
template< typename OtherDerived >
void SparseBlockMatrixBase<Derived>::cloneDimensions( const BlockMatrixBase< OtherDerived > &source )
{
	std::vector< unsigned > dims( source.rowsOfBlocks() ) ;

	for( unsigned i = 0 ; i < dims.size() ; ++i )
	{
		dims[i] = (unsigned) source.blockRows( i ) ;
	}

	setRows( dims ) ;
	dims.resize( source.colsOfBlocks() ) ;

	for( unsigned i = 0 ; i < dims.size() ; ++i )
	{
		dims[i] = (unsigned) source.blockCols( i )  ;
	}

	setCols( dims ) ;
}

template < typename Derived >
template < typename OtherDerived >
void SparseBlockMatrixBase<Derived>::cloneStructure( const SparseBlockMatrixBase< OtherDerived > &source )
{
	BOGUS_STATIC_ASSERT( static_cast<unsigned>(BlockMatrixTraits< Derived >::flags)
						 == static_cast< unsigned >(BlockMatrixTraits< OtherDerived >::flags),
						 OPERANDS_HAVE_INCONSISTENT_FLAGS ) ;

	rowMajorIndex() = source.rowMajorIndex() ;
	colMajorIndex() = source.colMajorIndex() ;

	m_cols = source.cols() ;
	m_rows = source.rows() ;
	m_blocks.resize( source.nBlocks() ) ;

	m_transposeBlocks.clear() ;
	m_transposeIndex.valid = false ;
}


}

#endif
