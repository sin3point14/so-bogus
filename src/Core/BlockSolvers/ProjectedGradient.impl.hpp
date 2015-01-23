/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_PROJECTED_GRADIENT_IMPL_HPP
#define BOGUS_PROJECTED_GRADIENT_IMPL_HPP


#include "ProjectedGradient.hpp"
#include "ConstrainedSolverBase.impl.hpp"

namespace bogus
{

template < typename BlockMatrixType >
template < typename NSLaw, typename RhsT, typename ResT >
typename ProjectedGradient< BlockMatrixType >::Scalar
ProjectedGradient< BlockMatrixType >::solve(
		const NSLaw &law, const RhsT &b, ResT &x ) const
{
	return solve< projected_gradient::Standard, NSLaw, RhsT, ResT >( law, b, x ) ;
}

template < typename BlockMatrixType >
template < projected_gradient::Variant variant, typename NSLaw, typename RhsT, typename ResT >
typename ProjectedGradient< BlockMatrixType >::Scalar
ProjectedGradient< BlockMatrixType >::solve(
		const NSLaw &law, const RhsT &b, ResT &x ) const
{

	typename GlobalProblemTraits::DynVector
			Mx ( b.rows() ),
			y  ( b.rows() ), // = Mx +b
			xs ( x.rows() )  // tentative new value for x
			;

	//Conjugation
	typename GlobalProblemTraits::DynVector
			dir,
			prev_proj,
			prev_dir ;
	Scalar prev_n2 = 0. ;


	//APGD
	Scalar theta_prev = 1. ;


	projectOnConstraints( law, x ) ;
	prev_proj = x ;

	// Unconstrained objective function
	Mx = (*m_matrix)*x ;
	Scalar J = x.dot( .5 * Mx + b ) ;

	Scalar res = -1,
		   alpha = 1. ;

	// Eval optimal alpha
	xs.setOnes() ;
	const Scalar nMx = ( (*m_matrix)*xs ).squaredNorm() ;
	if( nMx > 1.e-4 )
		alpha = b.rows() / nMx ;

	for( unsigned pgIter = 0 ; pgIter < m_maxIters ; ++pgIter )
	{
		// y = grad J = Mx+b
		y = Mx + b ;
		res = Base::eval( law, y, x ) ;

		this->m_callback.trigger( pgIter, res );
		if( res < m_tol ) break ;

		if(variant == projected_gradient::Conjugated)
		{

			xs = x - y ;
			projectOnConstraints( law, xs ) ;
			xs = x - xs ;

			//Conjugation
			const Scalar ng2 = xs.squaredNorm() ;
			if( prev_n2 == 0. ) {
				dir = -y ;
			} else {
				const Scalar num = 1 ? (ng2 - xs.dot(prev_proj)) : ng2 ;
				const Scalar den = 0 ? (prev_dir.dot( xs - prev_proj )) : prev_n2 ;
				const Scalar beta = std::max( 0., num/den );
				// Polak-Ribiere
				//const Scalar beta = std::max( 0., (ng2 - xs.dot(prev_proj)) ) / prev_n2 ;

				dir = beta * prev_dir - y;
			}

			prev_proj = xs ;
			prev_n2 = ng2 ;

		}


		Scalar Js = J ;

		alpha *= m_lsOptimisticFactor ;

		// Line-search
		for( unsigned lsIter = 0 ;
			 lsIter < m_lsIters ;
			 ++ lsIter, alpha *= m_lsPessimisticFactor )
		{
			if(variant == projected_gradient::Conjugated)
				xs = x + alpha * dir ;
			else
				xs = x - alpha * y ;

			projectOnConstraints( law, xs ) ;

			Mx = (*m_matrix) * xs ;
			Js = xs.dot( .5 * Mx + b ) ;

			const Scalar decr = (variant == projected_gradient::Conjugated)
					? - dir.dot( xs - x )
					:     y.dot(xs - x) ;

			if( Js < J + decr + .5 / alpha * ( xs - x ).squaredNorm() )
				break ;

		}

		if( variant == projected_gradient::APGD && (xs - x).dot( y ) <= 0. )
		{
			const Scalar theta_p2 = theta_prev * theta_prev ;
			const Scalar delta = std::sqrt( theta_p2 + 4 ) ;
			Scalar theta = .5 * ( theta_prev * delta - theta_p2 ) ;

			const Scalar beta = ( theta_prev - theta_p2 ) / (theta_p2 + theta ) ;
			x = xs + beta * (xs - prev_proj) ;

			Mx = (*m_matrix) * x ;
			J = x.dot( .5 * Mx + b ) ;

			prev_proj = xs ;
			theta_prev = theta ;
		} else {

			if(variant  == projected_gradient::Conjugated) {
				prev_dir = (xs - x) / alpha ;
			}

			x = xs ;
			J = Js ;

			theta_prev = 1. ; //APGD
		}

	}


	return res ;
}

template < typename BlockMatrixType >
template < typename NSLaw, typename VectorT >
void ProjectedGradient< BlockMatrixType >::projectOnConstraints(
		const NSLaw &law, VectorT &x ) const
{
	Segmenter< NSLaw::dimension, VectorT, typename BlockMatrixType::Index >
			xSegmenter( x, m_matrix->rowOffsets() ) ;

	const Index n = m_matrix->rowsOfBlocks() ;
	typename NSLaw::Traits::Vector lx ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for private( lx )
#endif
	for( Index i = 0 ; i < n ; ++ i )
	{
		lx = xSegmenter[ i ] ;
		law.projectOnConstraint( i, lx ) ;
		xSegmenter[ i ] = lx ;
	}

}

} //namespace bogus

#endif
