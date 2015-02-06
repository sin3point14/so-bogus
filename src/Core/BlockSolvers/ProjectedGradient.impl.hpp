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
	switch ( m_defaultVariant )
	{
	case projected_gradient::Standard:
		return solve< projected_gradient::Standard, NSLaw, RhsT, ResT >( law, b, x ) ;
	case projected_gradient::Conjugated:
		return solve< projected_gradient::Conjugated, NSLaw, RhsT, ResT >( law, b, x ) ;
	case projected_gradient::APGD:
		return solve< projected_gradient::APGD, NSLaw, RhsT, ResT >( law, b, x ) ;
	}

	return -1 ;
}

template < typename BlockMatrixType >
template < projected_gradient::Variant variant, typename NSLaw, typename RhsT, typename ResT >
typename ProjectedGradient< BlockMatrixType >::Scalar
ProjectedGradient< BlockMatrixType >::solve(
		const NSLaw &law, const RhsT &b, ResT &x ) const
{

	typename GlobalProblemTraits::DynVector
			Mx ( b.rows() ),
			y  ( b.rows() ), // = Mx +b  (gradient)
			xs ( x.rows() )  // tentative new value for x
			;

	//Conjugation
	typename GlobalProblemTraits::DynVector	dir, prev_dir, prev_proj ;
	Scalar prev_n2 = 0. ;


	//APGD
	Scalar theta_prev = 1. ;


	projectOnConstraints( law, x ) ;
	prev_proj = x ;

	// Unconstrained objective function
	Mx = (*m_matrix)*x ;
	Scalar J = x.dot( .5 * Mx + b ) ;

	// Eval optimal alpha (step size )
	Scalar alpha = 1. ;
	xs.setOnes() ;
	const Scalar nMx = ( (*m_matrix)*xs ).squaredNorm() ;
	if( nMx > 1.e-4 ) //Don't try too big alphas
		alpha = b.rows() / nMx ;

	// Best iterate storage
	typename GlobalProblemTraits::DynVector x_best ;
	Scalar res, min_res = -1 ;

	for( unsigned pgIter = 0 ; pgIter < m_maxIters ; ++pgIter )
	{
		// y = grad J = Mx+b
		y = Mx + b ;
		res = Base::eval( law, y, x ) ;

		this->m_callback.trigger( pgIter, res );
		if( 0 == pgIter || res < min_res ) {
			x_best = x ;
			min_res = res ;
		}
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
				const Scalar den = prev_dir.dot( xs - prev_proj ) ;
				const Scalar beta = ( den < 1.e-12 )
						// Polak-Ribiere
						? std::max( 0., (ng2 - xs.dot(prev_proj))) / prev_n2
						// Hestness-Stiefel
						: std::max( 0., (ng2 - xs.dot(prev_proj))) / den  ;


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

		if( variant == projected_gradient::APGD && (xs - prev_proj).dot( y ) < 0. )
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

		alpha = std::max(1.e-4, alpha) ;

	}

	x = x_best ;
	return min_res ;
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
