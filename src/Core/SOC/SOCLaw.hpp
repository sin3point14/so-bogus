/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver          
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>                       
 *
 * This Source Code Form is subject to the terms of the Mozilla Public 
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_SOCLAW_HPP
#define BOGUS_SOCLAW_HPP

#include "../BlockSolvers.fwd.hpp"
#include "FischerBurmeister.hpp"
#include "LocalSOCSolver.hpp"

#include <vector>

namespace bogus
{

template < typename LocalMatrixType, bool DeSaxceCOV,
		   local_soc_solver::Strategy Strat >
class SOCLaw
{
public:
	typedef ProblemTraits< LocalMatrixType > GlobalProblemTraits ;
	typedef LocalProblemTraits< GlobalProblemTraits::dimension, typename GlobalProblemTraits::Scalar > ProblemTraits ;
	typedef typename ProblemTraits::Scalar Scalar ;

	SOCLaw( const unsigned n, const double * mu ) ;

	template< typename VectorT, typename OtherVectorT >
	Scalar eval( const VectorT &x, const OtherVectorT &y ) const
	{
		typedef FischerBurmeister< ProblemTraits::dimension, typename ProblemTraits::Scalar, DeSaxceCOV > FBFunction ;

		assert( (unsigned) x.rows() == m_n * ProblemTraits::dimension ) ;
		assert( (unsigned) y.rows() == m_n * ProblemTraits::dimension ) ;

		Scalar sum = 0. ;
		typename ProblemTraits::Vector lx, ly, fb ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for private( lx, ly, fb ) reduction ( + : sum )
#endif
		for( int i = 0 ; i < (int) m_n ; ++ i )
		{
			lx = GlobalProblemTraits::segment( i, x ) ;
			ly = GlobalProblemTraits::segment( i, y ) ;
			FBFunction::compute( m_mu[i], lx, ly, fb ) ;
			sum += fb.squaredNorm() ;
		}

		return sum / ( 1 + m_n );
	}

	bool solveLocal(
			const unsigned problemIndex,
			const typename ProblemTraits::Matrix &A,
			const typename ProblemTraits::Vector &b,
			typename ProblemTraits::Vector &xm,
			const Scalar scaling
			) const ;

private:

	const double * m_mu ;
	const unsigned m_n ;
	Scalar m_localTol ;

} ;

}

#endif // SOCLAW_HPP
