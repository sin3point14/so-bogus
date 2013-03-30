#ifndef BOGUS_SPARSEPRODUCT_IMPL_HPP
#define BOGUS_SPARSEPRODUCT_IMPL_HPP

#include "Expressions.hpp"
#include "SparseBlockMatrix.hpp"

namespace bogus
{

template < typename Derived >
template < typename LhsDerived, typename RhsDerived >
void SparseBlockMatrixBase<Derived>::setFromProduct( const SparseBlockMatrixBase< LhsDerived >& lhs,
					 const SparseBlockMatrixBase< RhsDerived >& rhs,
					 bool lhsTransposed, bool rhsTransposed )
{
	typedef BlockMatrixTraits< LhsDerived > LhsTraits ;
	typedef BlockMatrixTraits< RhsDerived > RhsTraits ;

	TransposeMode transposeModeLhs ; TransposeMode transposeModeRhs ;
	SparseBlockIndex< > auxIndexLhs, auxIndexRhs ;

	assert( ! LhsTraits::is_symmetric ) ;
	assert( ! RhsTraits::is_symmetric ) ;

	const SparseBlockIndexBase &lhsIdx = getIndex( lhsTransposed, true, transposeModeLhs, auxIndexLhs ) ;
	const SparseBlockIndexBase &rhsIdx = getIndex( rhsTransposed, true, transposeModeRhs, auxIndexRhs ) ;

	if( lhsIdx.isCompressed() )
	{
		if( rhsIdx.isCompressed() )
		{
			setFromProduct( lhsIdx.asCompressed(), rhsIdx.asCompressed(),
							lhs.blocks(), rhs.blocks(), transposeModeLhs, transposeModeRhs ) ;
		} else {
			setFromProduct( lhsIdx.asCompressed(), rhsIdx.asUncompressed(),
							lhs.blocks(), rhs.blocks(), transposeModeLhs, transposeModeRhs ) ;
		}
	} else {
		if( rhsIdx.isCompressed() )
		{
			setFromProduct( lhsIdx.asUncompressed(), rhsIdx.asCompressed(),
							lhs.blocks(), rhs.blocks(), transposeModeLhs, transposeModeRhs ) ;
		} else {
			setFromProduct( lhsIdx.asUncompressed(), rhsIdx.asUncompressed(),
							lhs.blocks(), rhs.blocks(), transposeModeLhs, transposeModeRhs ) ;
		}
	}

}

template < typename Derived >
template < typename LhsIndex, typename RhsIndex, typename LhsBlock, typename RhsBlock  >
void SparseBlockMatrixBase<Derived>::setFromProduct(
		const LhsIndex &lhsIdx,
		const RhsIndex &rhsIdx,
		const std::vector< LhsBlock > &lhsData,
		const std::vector< RhsBlock > &rhsData,
		TransposeMode transposeLhs, TransposeMode transposeRhs
		)
{
	(void) lhsIdx ;
	(void) rhsIdx ;
	(void) lhsData ;
	(void) rhsData;
	(void) transposeLhs;
	(void) transposeRhs;
}


}

#endif // SPARSEPRODUCT_IMPL_HPP
