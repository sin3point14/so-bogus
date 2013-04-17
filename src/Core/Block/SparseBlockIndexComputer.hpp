#ifndef BOGUS_SPARSE_BLOCK_INDEX_COMPUTER_HPP
#define BOGUS_SPARSE_BLOCK_INDEX_COMPUTER_HPP

#include "CompoundSparseBlockIndex.hpp"

namespace bogus
{


template < typename MatrixType, bool Symmetric, bool ColWise, bool Transpose >
struct SparseBlockIndexComputer
{
	typedef BlockMatrixTraits< MatrixType > Traits ;
	enum { is_major = bool(ColWise ^ Transpose) == bool( Traits::is_col_major ) } ;
	typedef SparseBlockIndexGetter< MatrixType, is_major > Getter ;
	typedef typename Getter::ReturnType ReturnType ;

	SparseBlockIndexComputer( const MatrixType &matrix )
		: m_matrix( matrix )
	{}

	const ReturnType& get( )
	{
		return Getter::getOrCompute( m_matrix, m_aux ) ;
	}

private:
	const MatrixType m_matrix ;
	typename Traits::UncompressedIndexType m_aux ;
} ;

template < typename MatrixType, bool ColWise, bool Transpose >
struct SparseBlockIndexComputer< MatrixType, true, ColWise, Transpose >
{
	typedef BlockMatrixTraits< MatrixType > Traits ;
	typedef CompoundSparseBlockIndex< typename Traits::SparseIndexType, typename Traits::UncompressedIndexType >
	ReturnType ;

	SparseBlockIndexComputer( const MatrixType &matrix )
		: m_index( matrix.majorIndex(), matrix.minorIndex() )
	{}

	const ReturnType& get( )
	{
		return m_index ;
	}

private:
	ReturnType m_index ;
} ;


}

#endif
