/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_SECOND_ORDER_FWD_HPP
#define BOGUS_SECOND_ORDER_FWD_HPP

namespace bogus
{

//! Configuration properties of local Second Order Cone solver
namespace local_soc_solver
{
//! Strategy to be used by the local SOC solver.
/*! Note that some strategies may be unavailable for some loval problem types,
	in which case the solver will revert to the PureNewton strategy */
enum Strategy
{
	PureNewton                //!< Newton algorithm on the SOC FischerBurmeister function. \sa NonSmoothNewton
#ifndef BOGUS_WITHOUT_EIGEN
	,PureEnumerative          //!< Enumerative algorithm, such as describer in Appendix B of \cite DBB11
	,Hybrid                   //!< Newton algorithm, then Enumerative as failsafe
	,RevHybrid                //!< Enumerative algorithm, then Newton to refine the solution
#endif
} ;
}

template< unsigned Dimension, typename Scalar >
struct LocalProblemTraits ;

template < unsigned Dimension, typename Scalar, bool DeSaxceCOV,
#ifndef BOGUS_WITHOUT_EIGEN
		   local_soc_solver::Strategy Strat = local_soc_solver::RevHybrid  >
#else
		   local_soc_solver::Strategy Strat = local_soc_solver::PureNewton  >
#endif
class SOCLaw ;

}

#endif
