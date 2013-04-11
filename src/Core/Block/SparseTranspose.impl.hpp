#ifndef BOGUS_SPARSETRANSPOSE_IMPL_HPP
#define BOGUS_SPARSETRANSPOSE_IMPL_HPP

#include "SparseBlockMatrix.hpp"

namespace bogus
{

template < typename Derived >
template < typename OtherDerived >
Derived& SparseBlockMatrixBase<Derived>::operator=( const Transpose< SparseBlockMatrixBase< OtherDerived > > &tSource )
{
	const SparseBlockMatrixBase< OtherDerived >& source = tSource.matrix ;

	assert( static_cast< const void* >( this ) != static_cast< const void* >( &source ) ) ;

	m_nBlocks = source.nBlocks() ;

	if( Traits::is_symmetric )
	{
		m_majorIndex = source.majorIndex() ;
		m_minorIndex = source.minorIndex() ;
	} else {
		rowMajorIndex() = source.colMajorIndex() ;
		colMajorIndex() = source.rowMajorIndex() ;
	}

	this->m_cols = source.rows() ;
	this->m_rows = source.cols() ;

	if( m_majorIndex.valid )
	{
		m_transposeCached = source.transposeCached() ;
		this->m_blocks.resize( source.blocks().size() ) ;

		if( m_transposeCached )
		{
			this->m_blocks.resize( source.blocks().size() ) ;
			std::copy( source.blocks().begin(), source.blocks().end(), this->m_blocks.begin() ) ;
		} else {
			for( unsigned i = 0 ; i < this->m_blocks.size() ; ++i )
			{
				block( i ) = source.block(i).transpose() ;
			}
		}
	} else {
		// If we're here, this means that :
		//  - either both matrices have the same ordering
		//  -     or the major index of the destination iscompressed and cannot accomodate the source

		clear() ;
		this->m_blocks.reserve( source.blocks().size() ) ;

		assert( source.majorIndex().valid ) ;

		SparseBlockIndex< > uncompressed ;
		if( ( Traits::is_col_major && BlockMatrixTraits< OtherDerived >::is_col_major ) ||
				!( Traits::is_col_major || BlockMatrixTraits< OtherDerived >::is_col_major ) )
		{
			// Same col-major-ness
			uncompressed.setToTranspose( source.majorIndex() ) ;
		} else {
			uncompressed = source.majorIndex() ;
		}

		for( unsigned i = 0 ; i < uncompressed.outerSize() ; ++i )
		{
			for( typename SparseBlockIndex<>::InnerIterator src_it( uncompressed, i ) ;
				 src_it ; ++ src_it )
			{
				insertBackOuterInner( i, src_it.inner() ) = source.block( src_it.ptr() ).transpose() ;
			}
		}
		finalize() ;
	}

	return this->derived();
}
}

#endif // SPARSETRANSPOSE_IMPL_HPP
