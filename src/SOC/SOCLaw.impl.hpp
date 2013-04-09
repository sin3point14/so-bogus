#ifndef BOGUS_SOCLAW_IMPL_HPP
#define BOGUS_SOCLAW_IMPL_HPP

#include "SOCLaw.hpp"
#include "../Utils/NumTraits.hpp"

#include "FischerBurmeister.hpp"
#include "LocalSOCSolver.hpp"

namespace bogus {

template < typename LocalMatrixType, bool DeSaxceCOV, local_soc_solver::Strategy Strat >
SOCLaw< LocalMatrixType, DeSaxceCOV, Strat >::SOCLaw( const std::vector< double >& mu )
	: m_mu(mu), m_localTol( std::pow( NumTraits< Scalar >::epsilon(), .75 ) )
{
}


template < typename LocalMatrixType, bool DeSaxceCOV, local_soc_solver::Strategy Strat >
bool SOCLaw< LocalMatrixType, DeSaxceCOV, Strat >::solveLocal(
			const unsigned problemIndex,
			const typename ProblemTraits::Matrix &A,
			const typename ProblemTraits::Vector &b,
			typename ProblemTraits::Vector &x ) const
{
	typedef LocalSOCSolver< ProblemTraits::dimension, typename ProblemTraits::Scalar, DeSaxceCOV, Strat > LocalSolver ;
	return m_localTol > LocalSolver::solve(  A, b, x, m_mu[ problemIndex ], m_localTol ) ;
}

}


#endif
