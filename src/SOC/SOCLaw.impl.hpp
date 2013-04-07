#ifndef BOGUS_SOCLAW_IMPL_HPP
#define BOGUS_SOCLAW_IMPL_HPP

#include "SOCLaw.hpp"
#include "../Utils/NumTraits.hpp"

#include "FischerBurmeister.hpp"
#include "LocalSOCSolver.hpp"

namespace bogus {

template < typename LocalMatrixType, bool DeSaxceCOV >
SOCLaw< LocalMatrixType, DeSaxceCOV >::SOCLaw()
	: localTol( std::pow( NumTraits< Scalar >::epsilon(), .75 ) )
{
}


template < typename LocalMatrixType, bool DeSaxceCOV >
bool SOCLaw< LocalMatrixType, DeSaxceCOV >::solveLocal(
			const unsigned problemIndex,
			const typename ProblemTraits::Matrix &A,
			const typename ProblemTraits::Vector &b,
            typename ProblemTraits::Vector &x ) const
{
	typedef LocalSOCSolver< ProblemTraits::dimension, typename ProblemTraits::Scalar, DeSaxceCOV > LocalSolver ;
    return localTol > LocalSolver::solve(  A, b, x, mu[ problemIndex ], localTol ) ;
}

}


#endif
