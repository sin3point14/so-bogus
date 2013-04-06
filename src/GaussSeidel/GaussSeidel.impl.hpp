#ifndef BOGUS_GAUSS_SEIDEL_IMPL_HPP
#define BOGUS_GAUSS_SEIDEL_IMPL_HPP

#include "../GaussSeidel.hpp"

namespace bogus
{


template < typename BlockMatrixType >
GaussSeidel< BlockMatrixType >::GaussSeidel( const BlockMatrixBase< BlockMatrixType > & M )
	: m_matrix( m ), m_maxIters( 1000 )
{
	const unsigned n = M.rowsOfBlocks() ;
	m_localProblems.resize( n ) ;

	for( unsigned i = 0 ; i < n ; ++i )
	{
		m_localProblems[i].A = M.diagonal( i ) ;
		m_localProblems[i].scaling = m_localProblems[i].trace() ;
	}

}

template < typename BlockMatrixType >
template < typename NSLaw, typename Derived, typename OtherDerived >
NSLaw::ErrorType GaussSeidel< BlockMatrixType >::solve( const NSLaw &law,
							const Eigen::MatrixBase< Derived >&b,
							Eigen::MatrixBase< OtherDerived > &x ) const
{
	const unsigned d = ProblemTraits::dimension ;
	const unsigned n = m_localProblems.size() ;
	assert( n*d == b.rows() ) ;
	assert( n*d == x.rows() ) ;

	std::vector< unsigned > skip( n, 0 ) ;

	for( unsigned GSIter = 0 ; GSIter < m_maxIters ; ++GSIter )
	{
		LocalProblemType::Vector lb, lx, ldx ;
		for( unsigned i = 0 ; i < n ; ++ i )
		{
			if( skip[i] ) {
				--skip[i] ;
				continue ;
			}

			lb = b.segment< d >( d*i ) ;
			m_matrix.splitRowMultiply( i, x, lb ) ;
			lx = x.segment< d >( d*i ) ;
			ldx = lx ;

			law.solveLocal( m_localProblems[i].A, lb, lx ) ;

			ldx -= lx ;
			const Traits::Scalar scaledSkipTol = scaling * scaling * m_skipTol ;
			if( ldx.squaredNorm() < scaledSkipTol || lx.squaredNorm() < scaledSkipTol )
			{
				skip[i] = m_skipIters ;
			}
		}

		if( 0 == (GSIter % m_evalEvery ) )
		{

		}

	}

}




}


#endif
