#ifndef BOGUS_SPARSE_BLOCK_MATRIX_IMPL_HPP
#define BOGUS_SPARSE_BLOCK_MATRIX_IMPL_HPP

#include "SparseBlockMatrix.hpp"

#include <iostream>

namespace bogus {

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
void SparseBlockMatrixBase< Derived >::setInnerOffets(
		SparseIndexType &index, const std::vector<Index> &blockSizes)
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
	// Warning if Traits::is_compressed, the col major index will not be
	// valid unless the transpose Matrix is cached as well

	if ( m_minorIndex.valid ) return true ;

	SparseBlockIndex< > cmIndex ;
	computeMinorIndex( cmIndex ) ;

	m_minorIndex = cmIndex ;

	return m_minorIndex.valid ;
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::computeMinorIndex(SparseBlockIndex< > &cmIndex) const
{

	cmIndex.innerOffsets = m_minorIndex.innerOffsets ;
	cmIndex.resizeOuter( m_majorIndex.innerSize() );

	for ( unsigned i = 0 ; i < m_majorIndex.outerSize() ; ++ i )
	{
		// For a symmetric matrix, do not store diagonal block in col-major index
		for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( m_majorIndex, i ) ;
			 it && !( Traits::is_symmetric && it.inner() == i ) ; ++ it )
		{
			cmIndex.insertBack( it.inner(), i, it.ptr() );
		}
	}

	cmIndex.finalize() ;
}

template < typename Derived >
const SparseBlockIndex< >& SparseBlockMatrixBase< Derived >::getUncompressedMinorIndex(SparseBlockIndex< > &cmIndex) const
{
	if( !Traits::is_compressed && m_minorIndex.valid && !m_transposeCached ) return m_minorIndex.asUncompressed() ;
	computeMinorIndex( cmIndex ) ;
	return cmIndex ;
}

template < typename Derived >
void SparseBlockMatrixBase< Derived >::cacheTranspose()
{
	if ( m_transposeCached ) return ;

	SparseBlockIndex< > uncompressedIndex ;
	const SparseBlockIndex< >& minorIndex = getUncompressedMinorIndex( uncompressedIndex ) ;

	if ( !m_minorIndex.valid )
	{
		m_minorIndex = minorIndex ;
	}

	this->m_blocks.reserve( 2*m_nBlocks ) ;

	for ( unsigned i = 0 ; i < minorIndex.outerSize() ; ++ i )
	{
		typename SparseBlockIndex< >::InnerIterator uncompressed_it
			 ( minorIndex , i ) ;
		for( typename SparseIndexType::InnerIterator it( m_minorIndex, i ) ;
			 it ; ++ it, ++uncompressed_it )
		{
			allocateBlock() = block( uncompressed_it.ptr() ).transpose() ;
			m_minorIndex.setPtr( it, ( BlockPtr ) ( this->m_blocks.size() - 1 )  ) ;
		}
	}

	m_transposeCached = true ;
}

template < typename Derived >
const SparseBlockIndexBase& SparseBlockMatrixBase< Derived >::getIndex(const bool transpose, const bool colWise,
		SparseBlockIndex< >& aux ) const
{
	if ( bool(colWise ^ transpose) == bool( Traits::is_col_major ) )
	{
		// Native order
		return m_majorIndex ;
	} else {
		// Minor order
		if( m_minorIndex.valid && !m_transposeCached ) {
			return m_minorIndex ;
		} else {
			computeMinorIndex( aux ) ;
			return aux ;
		}
	}
}

template < typename Derived >
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::multiply( const RhsT& rhs, ResT& res, bool transposed ) const
{
	if( Traits::is_symmetric )
	{
		if( m_transposeCached )
		{
			for( Index i = 0 ; i < majorIndex().outerSize() ; ++i )
			{
				typename ResT::SegmentReturnType seg( rowSegment( res, i ) ) ;
				innerRowMultiply( rowMajorIndex(), i, rhs, seg ) ;
				innerRowMultiply( colMajorIndex(), i, rhs, seg ) ;
			}
		} else {
			for( Index i = 0 ; i < majorIndex().outerSize() ; ++i )
			{
				typename RhsT::ConstSegmentReturnType rhs_seg( majorIndex().innerSegment( rhs, i ) ) ;
				typename ResT::SegmentReturnType res_seg( majorIndex().innerSegment( res, i ) ) ;
				for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( m_majorIndex, i ) ;
					 it ; ++ it )
				{
					const BlockType &b = block( it.ptr() ) ;
					res_seg += b * m_minorIndex.innerSegment( rhs, it.inner() ) ;
					if( it.inner() != i )
						m_minorIndex.innerSegment( res, it.inner() ) += b.transpose() * rhs_seg  ;
				}
			}
		}

	} else if( transposed )
	{
		if ( Traits::is_col_major )  {
			for( Index i = 0 ; i < colsOfBlocks() ; ++i )
			{
				typename ResT::SegmentReturnType seg( colSegment( res, i ) ) ;
				innerRowTransposedMultiply( colMajorIndex(), i, rhs, seg ) ;
			}
		} else {
			if( m_transposeCached )
			{
				for( Index i = 0 ; i < colsOfBlocks() ; ++i )
				{
					typename ResT::SegmentReturnType seg( colSegment( res, i ) ) ;
					innerRowMultiply( colMajorIndex(), i, rhs, seg ) ;
				}
			} else {
				for( Index i = 0 ; i < rowsOfBlocks() ; ++i )
				{
					typename RhsT::ConstSegmentReturnType seg( rowSegment( rhs, i ) ) ;
					innerColTransposedMultiply( rowMajorIndex(), i, seg, res ) ;
				}
			}
		}
	} else if ( Traits::is_col_major )  {
		for( Index i = 0 ; i < colsOfBlocks() ; ++i )
		{
			typename RhsT::ConstSegmentReturnType seg( colSegment( rhs, i ) ) ;
			innerColMultiply( colMajorIndex(), i, seg, res ) ;
		}
	} else {
		for( Index i = 0 ; i < rowsOfBlocks() ; ++i )
		{
			typename ResT::SegmentReturnType seg( rowSegment( res, i ) ) ;
			innerRowMultiply( rowMajorIndex(), i, rhs, seg ) ;
		}
	}
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
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::splitRowMultiply( const Index row, const RhsT& rhs, ResT& res ) const
{
	assert( Traits::is_symmetric || !Traits::is_col_major ) ;

	for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( m_majorIndex, row ) ;
		 it ; ++ it )
	{
		if( it.inner() != row )
			res += block( it.ptr() ) * colSegment( rhs, it.inner() ) ;
	}
	if( Traits::is_symmetric )
	{
		if( m_transposeCached )
		{
			innerRowMultiply( m_minorIndex, row, rhs, res ) ;
		} else {
			assert( m_minorIndex.valid ) ;
			innerRowTransposedMultiply( m_minorIndex, row, rhs, res ) ;
		}
	}
}

template < typename Derived >
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::innerRowMultiply( const SparseIndexType &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const
{
	for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( index, outerIdx ) ;
		 it ; ++ it )
	{
		res += block( it.ptr() ) * index.innerSegment( rhs, it.inner() ) ;
	}
}

template < typename Derived >
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::innerRowTransposedMultiply( const SparseIndexType &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const
{
	for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( index, outerIdx ) ;
		 it ; ++ it )
	{
		res += block( it.ptr() ).transpose() * index.innerSegment( rhs, it.inner() ) ;
	}
}

template < typename Derived >
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::innerColMultiply( const SparseIndexType &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const
{
	for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( index, outerIdx ) ;
		 it ; ++ it )
	{
		index.innerSegment( res, it.inner() ) += block( it.ptr() ) * rhs  ;
	}
}

template < typename Derived >
template < typename RhsT, typename ResT >
void SparseBlockMatrixBase< Derived >::innerColTransposedMultiply( const SparseIndexType &index, const Index outerIdx, const RhsT& rhs, ResT& res ) const
{
	for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( index, outerIdx ) ;
		 it ; ++ it )
	{
		index.innerSegment( res, it.inner() ) += block( it.ptr() ).transpose() * rhs  ;
	}
}

template < typename Derived >
std::ostream& operator<<( std::ostream &out, const SparseBlockMatrixBase< Derived > &sbm )
{
	out << " Total rows: " << sbm.rows() << " / cols: " << sbm.cols() << std::endl ;
	for ( unsigned i = 0 ; i < sbm.majorIndex().outerSize() ; ++ i )
	{
		out << i << ": " ;
		for( typename SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( sbm.majorIndex(), i ) ;
			 it ; ++ it )
		{
			out << "(" << it.inner() << ";" << it.ptr() << ")" ;
		}
		out << std::endl ;
	}
	out << " Blocks (" << sbm.nBlocks() << ")" << std::endl ;
	for ( unsigned i = 0 ; i < sbm.nBlocks() ; ++ i )
	{
		out << sbm.block(i) << std::endl ;
		out << " --- " << std::endl ;
	}
	return out ;
}

template < typename Derived >
template < typename BlockT2 >
void SparseBlockMatrixBase<Derived>::cloneStructure( const SparseBlockMatrix< BlockT2, SparseBlockMatrixBase<Derived>::Traits::flags > &source )
{
	m_nBlocks = source.nBlocks() ;
	rowMajorIndex() = source.rowMajorIndex() ;
	colMajorIndex() = source.colMajorIndex() ;
	m_minorIndex.valid = ! source.transposeCached() ;
	m_transposeCached = false ;

	this->m_cols = source.cols() ;
	this->m_rows = source.rows() ;
	this->m_blocks.resize( m_nBlocks ) ;
}

template < typename Derived >
template < typename OtherDerived >
Derived& SparseBlockMatrixBase<Derived>::operator=( const SparseBlockMatrixBase< OtherDerived > &source )
{
	if( static_cast< const void* >( this ) == static_cast< const void* >( &source ) ) return this->derived() ;

	m_nBlocks = source.nBlocks() ;

	if( Traits::is_symmetric )
	{
		m_majorIndex = source.majorIndex() ;
		m_minorIndex = source.minorIndex() ;
	} else {
		rowMajorIndex() = source.rowMajorIndex() ;
		colMajorIndex() = source.colMajorIndex() ;
	}

	this->m_cols = source.cols() ;
	this->m_rows = source.rows() ;

	if( m_majorIndex.valid )
	{
		m_transposeCached = m_minorIndex.valid && source.transposeCached() ;

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
		if( Traits::is_col_major == BlockMatrixTraits< OtherDerived >::is_col_major )
		{
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

	return this->derived();
}

template < typename BlockT, int Flags >
class SparseBlockMatrix : public  SparseBlockMatrixBase< SparseBlockMatrix< BlockT, Flags > >
{
	typedef SparseBlockMatrixBase< SparseBlockMatrix< BlockT, Flags > > Base ;
public:
	SparseBlockMatrix() : Base() {}

	template < typename RhsT >
	SparseBlockMatrix( const RhsT& rhs ) : Base()
	{
		Base::operator= ( rhs ) ;
	}

	template < typename RhsT >
	SparseBlockMatrix< BlockT, Flags>& operator=( const RhsT& rhs )
	{
		return Base::operator= ( rhs ) ;
	}
} ;

}

#endif
