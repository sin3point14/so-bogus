#ifndef BOGUS_SPARSE_BLOCK_MATRIX_IMPL_HPP
#define BOGUS_SPARSE_BLOCK_MATRIX_IMPL_HPP

#include "SparseBlockMatrix.hpp"
#include "Expressions.hpp"

namespace bogus {

// Finalizer

template < bool Symmetric >
struct SparseBlockMatrixFinalizer
{
	template < typename Derived >
	static void finalize( SparseBlockMatrixBase< Derived >& ) { }
} ;
template < >
struct SparseBlockMatrixFinalizer<  true >
{
	template < typename Derived >
	static void finalize( SparseBlockMatrixBase< Derived >& matrix )
	{ matrix.computeMinorIndex() ; }
} ;

// Index getter

template < typename Derived, bool Major >
struct SparseBlockIndexGetter
{
	typedef SparseBlockMatrixBase< Derived > MatrixType ;
	typedef typename MatrixType::UncompressedIndexType ReturnType ;

	static ReturnType& get( MatrixType& matrix )
	{
		return matrix.m_minorIndex ;
	}

	static const ReturnType& get( const MatrixType& matrix )
	{
		return matrix.minorIndex() ;
	}

	static const ReturnType&
	getOrCompute( const MatrixType& matrix,
				  typename MatrixType::UncompressedIndexType& tempIndex
				  )
	{
		return matrix.getOrComputeMinorIndex( tempIndex) ;
	}
} ;

template < typename Derived >
struct SparseBlockIndexGetter< Derived, true >
{
	typedef SparseBlockMatrixBase< Derived > MatrixType ;
	typedef typename MatrixType::SparseIndexType ReturnType ;

	static ReturnType& get( MatrixType& matrix )
	{
		return matrix.m_majorIndex ;
	}
	static const ReturnType& get( const MatrixType& matrix )
	{
		return matrix.majorIndex() ;
	}

	static const ReturnType&
	getOrCompute( const MatrixType& matrix,
				  typename MatrixType::UncompressedIndexType& )
	{
		return matrix.majorIndex() ;
	}
} ;

// Sparse Block Matrix

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
	: m_nBlocks(0)
{
	setRows( 0, 0 ) ;
	setCols( 0, 0 ) ;
	m_transposeIndex.resizeOuter(0) ;
	m_transposeIndex.valid = false ;
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::setRows(
		const std::vector< Index > &rowsPerBlocks )
{
	setInnerOffets( colMajorIndex(), rowsPerBlocks );
	this->m_rows = colMajorIndex().innerOffsets.back() ;
	rowMajorIndex().resizeOuter( rowsPerBlocks.size() ) ;

	if( Traits::is_symmetric ) setCols( rowsPerBlocks ) ;
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::setCols(
		const std::vector< Index > &colsPerBlocks )
{
	setInnerOffets( rowMajorIndex(), colsPerBlocks );
	this->m_cols = rowMajorIndex().innerOffsets.back() ;

	colMajorIndex().resizeOuter( colsPerBlocks.size() ) ;
}

template < typename Derived >
template < typename IndexT >
void SparseBlockMatrixBase< Derived >::setInnerOffets(
		IndexT &index, const std::vector<Index> &blockSizes) const
{
	index.innerOffsets.resize( blockSizes.size() + 1 ) ;
	index.innerOffsets[0] = 0 ;
	for ( unsigned i = 0 ; i < blockSizes.size() ; ++ i )
	{
		index.innerOffsets[ i+1 ] = index.innerOffsets[ i ] + blockSizes[i] ;
	}
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::finalize()
{
	assert( m_majorIndex.valid ) ;
	m_majorIndex.finalize() ;

	Finalizer::finalize( *this ) ;
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::clear()
{
	m_majorIndex.clear() ;
	m_minorIndex.clear() ;
	m_nBlocks = 0 ;
	this->m_blocks.clear() ;
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

	if( Traits::is_symmetric )
	{
		cmIndex.resizeOuter( m_majorIndex.innerSize() );
		for ( unsigned i = 0 ; i < m_majorIndex.outerSize() ; ++ i )
		{
			// For a symmetric matrix, do not store diagonal block in col-major index
			for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( m_majorIndex, i ) ;
				 it && it.inner() != i ; ++ it )
			{
				cmIndex.insertBack( it.inner(), i, it.ptr() );
			}
		}

	} else {
		cmIndex.setToTranspose( m_majorIndex ) ;
	}

	cmIndex.finalize() ;
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
void SparseBlockMatrixBase< Derived >::cacheTranspose()
{
	if ( m_transposeIndex.valid ) return ;

	computeMinorIndex() ;

	this->m_blocks.resize( 2*m_nBlocks ) ;

	BlockPtr base = m_nBlocks ;
	std::vector< BlockPtr > ptrOffsets( m_minorIndex.outerSize() ) ;
	for( unsigned i = 0 ; i < m_minorIndex.outerSize() ; ++i )
	{
		ptrOffsets[i] = base ;
		base += m_minorIndex.size( i ) ;
	}

	m_transposeIndex = (const UncompressedIndexType& ) m_minorIndex ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for ( unsigned i = 0 ; i < m_minorIndex.outerSize() ; ++ i )
	{

		typename SparseBlockIndex< >::InnerIterator uncompressed_it
			 ( m_minorIndex , i ) ;
		for( typename SparseIndexType::InnerIterator it( m_transposeIndex, i ) ;
			 it ; ++ it, ++uncompressed_it )
		{
			const BlockPtr ptr = ptrOffsets[i]++ ;
			block( ptr ) = transpose_block( block( uncompressed_it.ptr() ) ) ;
			m_transposeIndex.setPtr( it, ptr ) ;
		}
	}

	assert( m_minorIndex.valid ) ;
	assert( m_transposeIndex.valid ) ;
}

template < typename Derived >
const typename SparseBlockMatrixBase< Derived >::BlockType& SparseBlockMatrixBase< Derived >::diagonal( const Index row ) const
{
	if( Traits::is_symmetric ) return this->block( m_majorIndex.last( row ) );

	for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( majorIndex(), row ) ;
		 it ; ++ it )
	{
		if( it.inner() == row )
			return this->m_blocks[ it.ptr() ] ;
	}
	assert( 0 && "No diagonal block" ) ;
	return this->m_blocks[ (BlockPtr)-1 ] ;
}

template < typename Derived >
const typename SparseBlockMatrixBase< Derived >::BlockType& SparseBlockMatrixBase< Derived >::block( Index row, Index col ) const
{
	if( Traits::is_col_major ) std::swap( row, col ) ;

	for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( majorIndex(), row ) ;
		 it ; ++ it )
	{
		if( it.inner() == col )
			return this->m_blocks[ it.ptr() ] ;
	}
	assert( 0 && "Bock does not exist" ) ;
	return this->m_blocks[ (BlockPtr)-1 ] ;
}


template < typename Derived >
template < typename BlockT2 >
void SparseBlockMatrixBase<Derived>::cloneStructure( const SparseBlockMatrix< BlockT2, SparseBlockMatrixBase<Derived>::Traits::flags > &source )
{
	m_nBlocks = source.nBlocks() ;
	rowMajorIndex() = source.rowMajorIndex() ;
	colMajorIndex() = source.colMajorIndex() ;
	m_transposeIndex.valid = false ;

	this->m_cols = source.cols() ;
	this->m_rows = source.rows() ;
	this->m_blocks.resize( m_nBlocks ) ;
}

template < typename Derived >
template < typename OtherDerived >
Derived& SparseBlockMatrixBase<Derived>::operator=( const SparseBlockMatrixBase< OtherDerived > &source )
{
	if( static_cast< const void* >( this ) == static_cast< const void* >( &source ) ) return this->derived() ;

	typedef typename SparseBlockMatrixBase< OtherDerived >::Traits OtherTraits ;
	const bool sameMajorness = ( (bool) Traits::is_col_major ) ^ !( OtherTraits::is_col_major ) ;

	m_nBlocks = source.nBlocks() ;

	if( Traits::is_symmetric || sameMajorness )
	{
		m_majorIndex = source.majorIndex() ;
		m_minorIndex = source.minorIndex() ;
		m_transposeIndex = source.transposeIndex() ;
	} else {
		rowMajorIndex() = source.rowMajorIndex() ;
		colMajorIndex() = source.colMajorIndex() ;
		m_transposeIndex.clear() ;
		m_transposeIndex.valid = false ;
	}

	this->m_cols = source.cols() ;
	this->m_rows = source.rows() ;

	if( m_majorIndex.valid )
	{
		this->m_blocks.resize( source.blocks().size() ) ;
		std::copy( source.blocks().begin(), source.blocks().end(), this->m_blocks.begin() ) ;
	} else {
		// If we're here, this means that :
		//  - either one matrix is column major, the other row major
		//  -     or the major index of the destination iscompressed and cannot accomodate the source

		clear() ;
		this->m_blocks.reserve( source.blocks().size() ) ;

		assert( source.majorIndex().valid ) ;

		SparseBlockIndex< > uncompressed ;
		if( sameMajorness )
		{
			// Same col-major-ness
			uncompressed = source.majorIndex() ;
		} else {
			uncompressed.setToTranspose( source.majorIndex() ) ;
		}

		for( unsigned i = 0 ; i < uncompressed.outerSize() ; ++i )
		{
			for( typename SparseBlockIndex<>::InnerIterator src_it( uncompressed, i ) ;
				 src_it ; ++ src_it )
			{
				insertBackOuterInner( i, src_it.inner() ) = source.block( src_it.ptr() ) ;
			}
		}
		finalize() ;
	}

	Finalizer::finalize( *this ) ;

	return this->derived();
}


}

#endif
