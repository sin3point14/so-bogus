/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver          
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>                       
 *
 * This Source Code Form is subject to the terms of the Mozilla Public 
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_LOCAL_SOC_SOLVER_HPP
#define BOGUS_LOCAL_SOC_SOLVER_HPP

#include "../SecondOrder.fwd.hpp"

namespace bogus {


template< unsigned Dimension, typename Scalar, bool DeSaxceCOV,
		   local_soc_solver::Strategy Strat = local_soc_solver::Hybrid  >
struct LocalSOCSolver
{

  typedef LocalProblemTraits< Dimension, Scalar > Traits ;
  typedef typename Traits::Vector Vector ;
  typedef typename Traits::Matrix Matrix ;

  static Scalar solve(
		  const typename Traits::Matrix &A,
		  const typename Traits::Vector &b,
		  typename Traits::Vector &x,
		  const Scalar mu, const Scalar tol,
		  const Scalar scaling = 1
		  ) ;

} ;

}

#endif
