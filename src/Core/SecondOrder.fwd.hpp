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

namespace local_soc_solver
{
enum Strategy
{
	PureNewton
#ifndef BOGUS_WITHOUT_EIGEN
	,PureEnumerative
	,Hybrid
	,RevHybrid
#endif
} ;
}

template< unsigned Dimension, typename Scalar >
struct LocalProblemTraits ;

template < typename LocalMatrixType, bool DeSaxceCOV,
		   local_soc_solver::Strategy Strat = local_soc_solver::Hybrid  >
class SOCLaw ;

}

#endif
