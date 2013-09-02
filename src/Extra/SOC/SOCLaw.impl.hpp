/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_SOCLAW_IMPL_HPP
#define BOGUS_SOCLAW_IMPL_HPP

#include "SOCLaw.hpp"
#include "../../Core/Utils/NumTraits.hpp"

#include "FischerBurmeister.hpp"
#include "LocalSOCSolver.hpp"
#include "LocalSOCSolver.impl.hpp"

namespace bogus {

template < unsigned Dimension, typename Scalar, bool DeSaxceCOV, local_soc_solver::Strategy Strat >
SOCLaw< Dimension, Scalar, DeSaxceCOV, Strat >::SOCLaw(const unsigned n, const double *mu )
	: m_mu(mu), m_n(n), m_localTol( std::pow( NumTraits< Scalar >::epsilon(), .75 ) )
{
}


template < unsigned Dimension, typename Scalar, bool DeSaxceCOV, local_soc_solver::Strategy Strat >
bool SOCLaw< Dimension, Scalar, DeSaxceCOV, Strat >::solveLocal(const unsigned problemIndex,
			const typename Traits::Matrix &A,
			const typename Traits::Vector &b,
			typename Traits::Vector &xm , const Scalar scaling ) const
{
	typedef LocalSOCSolver< Traits::dimension, typename Traits::Scalar, DeSaxceCOV, Strat > LocalSolver ;
	return m_localTol > LocalSolver::solve(  A, b, xm, m_mu[ problemIndex ], m_localTol, scaling ) ;
}

}


#endif
