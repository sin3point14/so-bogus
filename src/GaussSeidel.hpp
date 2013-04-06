#ifndef BOGUS_GAUSS_SEIDEL_HPP
#define BOGUS_GAUSS_SEIDEL_HPP

#include "Block/BlockMatrix.hpp"
#include "GaussSeidel/LocalProblem.hpp"

namespace bogus
{

template < typename BlockMatrixType >
class GaussSeidel
{
public:

	typedef typename BlockMatrixTraits< BlockMatrixType >::BlockType LocalMatrixType ;
	typedef LocalProblem< LocalMatrixType > LocalProblemType ;
	typedef typename LocalProblemType::Traits ProblemTraits ;

	explicit GaussSeidel( const BlockMatrixBase< BlockMatrixType > & M ) ;

	template < typename Derived >
	void setRhs( const Eigen::MatrixBase< Derived >&b,

	template < typename NSLaw, typename Derived, typename OtherDerived >
	NSLaw::ErrorType solve( const NSLaw &law,
							Eigen::MatrixBase< OtherDerived > &x ) const ;

private:
	const BlockMatrixBase< Derived > & m_matrix ;
	std::vector< LocalProblem > m_localProblems ;

	unsigned m_maxIters ;
	unsigned m_evalEvery ;
	ProblemTraits::Scalar m_skipTol ;
	ProblemTraits::Scalar m_skipIters ;
} ;

}


#endif
