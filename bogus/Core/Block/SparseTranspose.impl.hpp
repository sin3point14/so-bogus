/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_SPARSETRANSPOSE_IMPL_HPP
#define BOGUS_SPARSETRANSPOSE_IMPL_HPP

#include "SparseBlockMatrix.hpp"
#include "BlockTranspose.hpp"

namespace bogus
{

template < typename Derived >
template < typename OtherDerived >
Derived& SparseBlockMatrixBase<Derived>::operator=( const Transpose< SparseBlockMatrixBase< OtherDerived > > &tSource )
{
	const SparseBlockMatrixBase< OtherDerived >& source = tSource.matrix ;

	assert( static_cast< const void* >( this ) != static_cast< const void* >( &source ) ) ;

	typedef typename SparseBlockMatrixBase< OtherDerived >::Traits OtherTraits ;

	if( OtherTraits::is_symmetric )
	{
		return derived() = source ;
	}

	// sameMajorness will be true if the majorIndex of the untransposed source can be reused for the dest
	const bool sameMajorness =
			( (bool) Traits::is_col_major ) != ( (bool) OtherTraits::is_col_major ) ;

	m_nBlocks = source.nBlocks() ;

	bool needTranspose = false ;

	if( source.transposeIndex().valid && !sameMajorness )
	{
		m_majorIndex = source.transposeIndex() ;
		m_transposeIndex = source.majorIndex() ;
		m_minorIndex.innerOffsets = m_transposeIndex.innerOffsets ;
		m_minorIndex.clear() ;
		m_minorIndex.valid = false ;
	} else {
		needTranspose = true ;
		rowMajorIndex() = source.colMajorIndex() ;
		colMajorIndex() = source.rowMajorIndex() ;

		// Will be valid if sameMajorness && source.transposeIndex().valid
		m_transposeIndex = source.transposeIndex()  ;
	}

	if( Traits::is_symmetric )
	{
		// Dest is symmetric and source is not -- we may have too many blocks
		// Assume nothing is valid
		m_majorIndex.valid = false ;
	}

	m_cols = source.rows() ;
	m_rows = source.cols() ;

	if( m_majorIndex.valid )
	{
		m_blocks.resize( source.blocks().size() ) ;

		if( needTranspose )
		{
			for( unsigned i = 0 ; i < m_blocks.size() ; ++i )
			{
				block( i ) = transpose_block( source.block(i) ) ;
			}
		} else {
			std::copy( source.blocks().begin(), source.blocks().end(), m_blocks.begin() ) ;
		}

		Finalizer::finalize( *this ) ;
	} else {
		// If we're here, this means that :
		//  - either both matrices have the same ordering
		//  -     or the major index of the destination iscompressed and cannot accomodate the source

		clear() ;
		m_blocks.reserve( source.blocks().size() ) ;

		assert( source.majorIndex().valid ) ;

		UncompressedIndexType uncompressed ;
		if( sameMajorness )
		{
			// Same col-major-ness
			uncompressed = source.majorIndex() ;
		} else {
			uncompressed.setToTranspose( source.majorIndex() ) ;
		}

		for( Index i = 0 ; i < uncompressed.outerSize() ; ++i )
		{
			for( typename UncompressedIndexType::InnerIterator src_it( uncompressed, i ) ;
				 src_it ; ++ src_it )
			{
				if( Traits::is_symmetric && i < src_it.inner() ) break ;
				insertBackOuterInner( i, src_it.inner() ) = transpose_block( source.block( src_it.ptr() ) ) ;
			}
		}
		finalize() ;
	}

	return derived();
}
}

#endif // SPARSETRANSPOSE_IMPL_HPP
