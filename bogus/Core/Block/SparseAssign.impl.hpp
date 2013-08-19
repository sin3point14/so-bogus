/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_SPARSEASSIGN_IMPL_HPP
#define BOGUS_SPARSEASSIGN_IMPL_HPP

#include "Access.hpp"

#include "SparseBlockMatrix.hpp"
#include "SparseBlockIndexComputer.hpp"

namespace bogus
{

template< bool Transpose >
struct BlockCopier
{
	template < typename BlockT1, typename BlockT2, typename ScalarT >
	static void copy( BlockT1* dest, const BlockT2* source, int n, ScalarT scale )
	{
		BlockGetter< Transpose > getter ;
		if( scale == 1 )
		{
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
			for( int i = 0 ; i < n ; ++i )
				dest[i] = getter.get( source[i] ) ;
		} else {
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
			for( int i = 0 ; i < n ; ++i )
				dest[i] = scale * getter.get( source[i] ) ;
		}
	}
} ;

template < typename Derived >
Derived& SparseBlockMatrixBase<Derived>::operator= ( const SparseBlockMatrixBase &source )
{
	m_rows = source.rows() ;
	m_cols = source.cols() ;
	m_blocks = source.blocks() ;

	m_nBlocks = source.nBlocks() ;
	m_majorIndex = source.majorIndex() ;
	m_minorIndex = source.minorIndex() ;
	m_transposeIndex = source.transposeIndex() ;

	return derived() ;
}

template < typename Derived >
template < bool Transpose, typename OtherDerived >
Derived& SparseBlockMatrixBase<Derived>::assign( const SparseBlockMatrixBase< OtherDerived > &source, Scalar scale )
{
	BOGUS_STATIC_ASSERT( !Transpose || IsTransposable< typename OtherDerived::BlockType >::Value,
						 TRANSPOSE_IS_NOT_DEFINED_FOR_THIS_BLOCK_TYPE
	) ;

	if( static_cast< const void* >( this ) == static_cast< const void* >( &source ) ) return derived() ;

	typedef typename SparseBlockMatrixBase< OtherDerived >::Traits OtherTraits ;
	const bool sameMajorness = Transpose ^ ( ( (bool) Traits::is_col_major )  == ( (bool) OtherTraits::is_col_major ) ) ;
	const bool sameSymmetry =  ( (bool) Traits::is_symmetric )  == ( (bool) OtherTraits::is_symmetric ) ;
	bool useTransposeIndex = false ;

	if( ( Traits::is_symmetric && OtherTraits::is_symmetric ) || sameMajorness )
	{
		m_majorIndex = source.majorIndex() ;
		m_minorIndex = source.minorIndex() ;
		m_transposeIndex = source.transposeIndex() ;
	}
	else if ( Transpose && source.transposeIndex().valid )
	{
		m_majorIndex = source.transposeIndex() ;
		m_transposeIndex = source.majorIndex() ;
		m_minorIndex.innerOffsets = m_transposeIndex.innerOffsets ;
		m_minorIndex.clear() ;
		m_minorIndex.valid = false ;
		useTransposeIndex = true ;
	} else {
		if( Transpose )
		{
			rowMajorIndex() = source.colMajorIndex() ;
			colMajorIndex() = source.rowMajorIndex() ;
		} else {
			rowMajorIndex() = source.rowMajorIndex() ;
			colMajorIndex() = source.colMajorIndex() ;
		}
		m_transposeIndex.clear() ;
		m_transposeIndex.valid = false ;
	}

	m_majorIndex.valid &= sameSymmetry ;

	m_cols = source.cols() ;
	m_rows = source.rows() ;
	if( Transpose ) std::swap( m_cols, m_rows ) ;

	if( m_majorIndex.valid )
	{
		m_nBlocks = source.nBlocks() ;
		m_blocks.resize( source.blocks().size() ) ;

		if( useTransposeIndex )
		{
			BlockCopier< false >::copy( this->data(), source.data(), source.blocks().size(), scale ) ;
		} else {
			const bool needTranspose = Transpose ^ ( Traits::is_symmetric && OtherTraits::is_symmetric && !sameMajorness ) ;
			BlockCopier< needTranspose >::copy( this->data(), source.data(), source.blocks().size(), scale ) ;
		}

		Finalizer::finalize( *this ) ;
	} else {

		clear() ;
		reserve( source.blocks().size() ) ;

		SparseBlockIndexComputer< OtherDerived, OtherTraits::is_symmetric, Traits::is_col_major, Transpose >
				indexComputer( source ) ;
		typedef typename SparseBlockIndexComputer< OtherDerived, OtherTraits::is_symmetric, Traits::is_col_major, Transpose >::ReturnType
				SourceIndexType ;
		const SourceIndexType &sourceIndex = indexComputer.get() ;

		BlockTransposeOption< OtherTraits::is_symmetric, Transpose > blockGetter ;

		assert( sourceIndex.valid ) ;

		for( Index i = 0 ; i < sourceIndex.outerSize() ; ++i )
		{
			for( typename SourceIndexType::InnerIterator src_it( sourceIndex, i ) ;
				 src_it && !( Traits::is_symmetric && i < src_it.inner() ) ; ++ src_it )
			{
				const bool afterDiag = ( (bool) Traits::is_col_major )  == ( (bool) OtherTraits::is_col_major )
						 ? (src_it.inner() > i) : (src_it.inner() < i ) ;
				insertBackOuterInner( i, src_it.inner() ) = scale *
						blockGetter.get( source.block( src_it.ptr() ), afterDiag ) ;
			}
		}
		finalize() ;
	}


	return derived();
}

} //namespace bogus

#endif // SPARSETRANSPOSE_IMPL_HPP
