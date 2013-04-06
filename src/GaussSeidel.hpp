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
	typedef LocalProblemTraits< LocalMatrixType > ProblemTraits ;
	typedef typename ProblemTraits::Scalar Scalar ;

	explicit GaussSeidel( const BlockMatrixBase< BlockMatrixType > & M ) ;

	template < typename NSLaw, typename Derived, typename OtherDerived >
	Scalar solve( const NSLaw &law,
							const Eigen::MatrixBase< Derived >&b,
							Eigen::MatrixBase< OtherDerived > &x ) const ;

private:
	const BlockMatrixBase< BlockMatrixType > & m_matrix ;
	std::vector< LocalMatrixType > m_localMatrices ;
	typename ProblemTraits::DynVector m_scaling ;

	unsigned m_maxIters ;
	Scalar m_tol ;

	unsigned m_evalEvery ;
	Scalar m_skipTol ;
	Scalar m_skipIters ;
} ;

} //namespace bogus

#include "GaussSeidel/GaussSeidel.impl.hpp"


#endif
