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
		// If we're here, this means that the matrices are either both row-major or both column-Major

		clear() ;
		this->m_blocks.reserve( source.blocks().size() ) ;

		typedef typename BlockMatrixTraits< OtherDerived >::SparseIndexType SourceIndexType ;
		const SourceIndexType & sourceIndex = source.majorIndex() ;

		for( unsigned i = 0 ; i < sourceIndex.outerSize() ; ++i )
		{
			for( typename SourceIndexType::InnerIterator src_it( sourceIndex, i ) ;
				 src_it ; ++ src_it )
			{
				insertBackOuterInner( src_it.inner(), i ) = source.block( src_it.ptr() ).transpose() ;
			}
		}
		finalize() ;
	}

	return this->derived();
}
}

#endif // SPARSETRANSPOSE_IMPL_HPP
