/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver          
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>                       
 *
 * This Source Code Form is subject to the terms of the Mozilla Public 
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_SOCLAW_IMPL_HPP
#define BOGUS_SOCLAW_IMPL_HPP

#include "SOCLaw.hpp"
#include "../Utils/NumTraits.hpp"

#include "FischerBurmeister.hpp"
#include "LocalSOCSolver.hpp"
#include "LocalSOCSolver.impl.hpp"

namespace bogus {

template < typename LocalMatrixType, bool DeSaxceCOV, local_soc_solver::Strategy Strat >
SOCLaw< LocalMatrixType, DeSaxceCOV, Strat >::SOCLaw(const unsigned n, const double *mu )
	: m_mu(mu), m_n(n), m_localTol( std::pow( NumTraits< Scalar >::epsilon(), .75 ) )
{
}


template < typename LocalMatrixType, bool DeSaxceCOV, local_soc_solver::Strategy Strat >
bool SOCLaw< LocalMatrixType, DeSaxceCOV, Strat >::solveLocal(const unsigned problemIndex,
			const typename ProblemTraits::Matrix &A,
			const typename ProblemTraits::Vector &b,
			typename ProblemTraits::Vector &xm , const Scalar scaling ) const
{
	typedef LocalSOCSolver< ProblemTraits::dimension, typename ProblemTraits::Scalar, DeSaxceCOV, Strat > LocalSolver ;
	return m_localTol > LocalSolver::solve(  A, b, xm, m_mu[ problemIndex ], m_localTol, scaling ) ;
}

}


#endif
