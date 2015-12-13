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


namespace pg_impl 
{

template<projected_gradient::Variant variant>
struct PgMethod {

	// variant = Descent, APGD

	template < typename BlockMatrixType, typename NSLaw, typename RhsT, typename ResT >
	static typename ProjectedGradient< BlockMatrixType >::Scalar
	solve( const ProjectedGradient< BlockMatrixType >& pg, 
			const NSLaw &law, const RhsT &b, ResT &x ) 
	{
		typedef ProjectedGradient< BlockMatrixType > PgType ;
		typedef typename PgType::Scalar Scalar ;
		typename PgType::GlobalProblemTraits::DynVector
			Mx ( b.rows() ),
			y  ( b.rows() ), // = Mx +b  (gradient)
			xs ( x.rows() ), // tentative new value for x
			prev_proj,
			x_best
			;
		
		pg.projectOnConstraints( law, x ) ;

		// Unconstrained objective function
		Mx = pg.matrix()*x ;
		Scalar J = x.dot( .5 * Mx + b ) ;

		//APGD
		Scalar theta_prev = 1. ;

		prev_proj = x ;

		// Eval optimal alpha (step size )
		Scalar alpha = 1. ;
		xs.setOnes() ;
		const Scalar nMx = ( pg.matrix()*(x-xs) ).squaredNorm() ;
		if( nMx > 1.e-16 ) //Don't try too big alphas
			alpha = std::min( 1., (x-xs).squaredNorm() / nMx ) ;
		

		// Best iterate storage
		Scalar res, min_res = -1 ;

		for( unsigned pgIter = 0 ; pgIter < pg.maxIters() ; ++pgIter )
			{
			// y = grad J = Mx+b
			// The residual should be evaluated in prev_proj instead of x
			// However this would been one more matrix product per iterations
			y = Mx + b ;
			res = pg.eval( law, y, x ) ;

			pg.callback().trigger( pgIter, res );
			if( 0 == pgIter || res < min_res ) {
				x_best = x ;
				min_res = res ;
			}
			if( res < pg.tol() ) break ;


			Scalar Js = J ;
			alpha *= pg.lineSearchOptimisticFactor() ;

			// Line-search
			for( unsigned lsIter = 0 ;
				 lsIter < pg.lineSearchIterations() ;
				 ++ lsIter, alpha *= pg.lineSearchPessimisticFactor() )
			{
				xs = x - alpha * y ;

				pg.projectOnConstraints( law, xs ) ;

				Mx = pg.matrix() * xs ;
				Js = xs.dot( .5 * Mx + b ) ;

				const Scalar decr = y.dot(xs - x) ;

				if( Js < J + decr + .5 / alpha * ( xs - x ).squaredNorm() )
					break ;

			}

			const Scalar decr = (xs - prev_proj).dot( y ) ;

			if( variant == projected_gradient::APGD && decr < 0 )
			{
				const Scalar theta_p2 = theta_prev * theta_prev ;
				const Scalar delta = std::sqrt( theta_p2 + 4 ) ;
				Scalar theta = .5 * ( theta_prev * delta - theta_p2 ) ;

				const Scalar beta = ( theta_prev - theta_p2 ) / (theta_p2 + theta ) ;
				x = xs + beta * (xs - prev_proj) ;

				theta_prev = theta ;

				Mx = pg.matrix() * x ;
				J = x.dot( .5 * Mx + b ) ;

			} else {

				x = xs ;
				J = Js ;

				theta_prev = 1. ; //APGD
			}

			prev_proj = xs ;

		}

		x = x_best ;
		return min_res ;

	}

} ;

template< >
struct PgMethod< projected_gradient::Standard > {

	template < typename BlockMatrixType, typename NSLaw, typename RhsT, typename ResT >
	static typename ProjectedGradient< BlockMatrixType >::Scalar
	solve( const ProjectedGradient< BlockMatrixType >& pg, 
			const NSLaw &law, const RhsT &b, ResT &x ) 
	{
		typedef ProjectedGradient< BlockMatrixType > PgType ;
		typedef typename PgType::Scalar Scalar ;
		typename PgType::GlobalProblemTraits::DynVector
			Mx ( b.rows() ),
			y  ( b.rows() ), // = Mx +b  (gradient)
			xs ( x.rows() ), // tentative new value for x
			proj_grad( b.rows() ),
			x_best
			;
		
		pg.projectOnConstraints( law, x ) ;

		// Unconstrained objective function
		Mx = pg.matrix()*x ;
		Scalar J = x.dot( .5 * Mx + b ) ;

		Scalar alpha = 1 ;
		
		// Best iterate storage
		Scalar res, min_res = -1 ;

		for( unsigned pgIter = 0 ; pgIter < pg.maxIters() ; ++pgIter )
		{
			y = Mx + b ;
			res = pg.eval( law, y, x ) ;

			pg.callback().trigger( pgIter, res );
			if( 0 == pgIter || res < min_res ) {
				x_best = x ;
				min_res = res ;
			}
			if( res < pg.tol() ) break ;
			
			xs = x - y ; 
			pg.projectOnConstraints( law, xs ) ;

			proj_grad = (x - xs) ;
			const Scalar beta = - y.dot(proj_grad)
					/ ( proj_grad.dot( pg.matrix()*proj_grad ) );

			// Line-search
			Scalar Js = J ;
			alpha = std::min(1., alpha*pg.lineSearchOptimisticFactor() ) ;
			for( unsigned lsIter = 0 ;
				 lsIter < pg.lineSearchIterations() ;
				 ++ lsIter, alpha *= pg.lineSearchPessimisticFactor() )
			{
				xs = x + alpha * beta * proj_grad ;
				pg.projectOnConstraints( law, xs ) ;

				Mx = pg.matrix() * xs ;
				Js = xs.dot( .5 * Mx + b ) ;

				const Scalar decr = y.dot(xs - x) ;

				if( Js < J + pg.lineSearchArmijoCoefficient() *decr)
					break ;
			}

			x = xs ;
			J = Js ;
		}

		x = x_best ;
		return min_res ;
	}

} ;

template<>
struct PgMethod< projected_gradient::Conjugated > {

	template < typename BlockMatrixType, typename NSLaw, typename RhsT, typename ResT >
	static typename ProjectedGradient< BlockMatrixType >::Scalar
	solve( const ProjectedGradient< BlockMatrixType >& pg, 
			const NSLaw &law, const RhsT &b, ResT &x ) 
	{
		typedef ProjectedGradient< BlockMatrixType > PgType ;
		typedef typename PgType::Scalar Scalar ;
		typename PgType::GlobalProblemTraits::DynVector
			Mx ( b.rows() ),
			y  ( b.rows() ), // = Mx +b  (gradient)
			xs ( x.rows() ), // tentative new value for x
			dir, prev_dir, prev_proj,
			x_best
			;
		
		pg.projectOnConstraints( law, x ) ;

		// Unconstrained objective function
		Mx = pg.matrix()*x ;
		Scalar J = x.dot( .5 * Mx + b ) ;

		prev_proj = x ;

		Scalar prev_n2 = 0 ;
		Scalar alpha = 1. ;

		// Best iterate storage
		Scalar res, min_res = -1 ;

		for( unsigned pgIter = 0 ; pgIter < pg.maxIters() ; ++pgIter )
		{
			// y = grad J = Mx+b
			y = Mx + b ;
			res = pg.eval( law, y, x ) ;

			pg.callback().trigger( pgIter, res );
			if( 0 == pgIter || res < min_res ) {
				x_best = x ;
				min_res = res ;
			}
			if( res < pg.tol() ) break ;

			{

				xs = x - y ;
				pg.projectOnConstraints( law, xs ) ;
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

					if( dir.dot(y) >= 0. )
						dir = -y ;
				}

				prev_proj = xs ;
				prev_n2 = ng2 ;

			}

			Scalar Js = J ;
			alpha *= pg.lineSearchOptimisticFactor() ;

			// Line-search
			for( unsigned lsIter = 0 ;
				 lsIter < pg.lineSearchIterations() ;
				 ++ lsIter, alpha *= pg.lineSearchPessimisticFactor() )
			{
				xs = x + alpha * dir ;

				pg.projectOnConstraints( law, xs ) ;

				Mx = pg.matrix() * xs ;
				Js = xs.dot( .5 * Mx + b ) ;

				const Scalar decr = y.dot( xs - x ) ;

				if( Js < J + pg.lineSearchArmijoCoefficient()*decr ) 
					break ;
			}


			prev_dir = (xs - x) / alpha ;

			x = xs ;
			J = Js ;

		}

		x = x_best ;
		return min_res ;

	}

} ;

template< >
struct PgMethod< projected_gradient::SPG > {

	template < typename BlockMatrixType, typename NSLaw, typename RhsT, typename ResT >
	static typename ProjectedGradient< BlockMatrixType >::Scalar
	solve( const ProjectedGradient< BlockMatrixType >& pg, 
			const NSLaw &law, const RhsT &b, ResT &x ) 
	{
		typedef ProjectedGradient< BlockMatrixType > PgType ;
		typedef typename PgType::Scalar Scalar ;
		typename PgType::GlobalProblemTraits::DynVector
			Mx ( b.rows() ),
			y  ( b.rows() ), // = Mx +b  (gradient)
			xs ( x.rows() ), // tentative new value for x
			z ( x.rows() ), 
			s ( x.rows() ), 
			prev_proj,
			x_best
			;
		
		pg.projectOnConstraints( law, x ) ;
		prev_proj = x ;

		// Unconstrained objective function
		Mx = pg.matrix()*x ;
		Scalar J = x.dot( .5 * Mx + b ) ;

		Scalar theta_prev = 1, q = 0 ;
		Scalar alpha = 1, lambda = 1 ;
		Scalar a_min = 1.e-6, a_max = 1.e6 ;
		
		// Eval optimal alpha
		xs.setOnes() ;
		const Scalar nMx = ( pg.matrix()*(x-xs) ).squaredNorm() ;
		if( nMx > 1.e-16 ) //Don't try too big alphas
			alpha = std::min( 1., (x-xs).squaredNorm() / nMx ) ;
		
		// Best iterate storage
		Scalar res, min_res = -1 ;

		for( unsigned pgIter = 0 ; pgIter < pg.maxIters() ; ++pgIter )
		{
			y = Mx + b ;
			res = pg.eval( law, y, x ) ;

			pg.callback().trigger( pgIter, res );
			if( 0 == pgIter || res < min_res ) {
				x_best = x ;
				min_res = res ;
			}
			if( res < pg.tol() ) break ;
			
			xs = x - alpha * y ; 
			pg.projectOnConstraints( law, xs ) ;

			if( (xs - prev_proj).dot( y ) < 0. )
			{
				const Scalar theta_p2 = theta_prev * theta_prev ;
				const Scalar bb = theta_p2 - q ;
				const Scalar delta = std::sqrt( bb*bb + 4*theta_p2 ) ;
				Scalar theta = .5 * ( delta - bb )  ;

				const Scalar beta = ( theta_prev - theta_p2 ) / (theta_p2 + theta ) ;

				s = ( xs + beta * (xs - prev_proj) ) - x ;
				x += s ;

				theta_prev = theta ;

			} else {
				s = xs - x ;
				x = xs ;

				theta_prev = 1. ;
			}

			Mx = pg.matrix() * x ;
			prev_proj = xs ;

			z = Mx + b - y ;
			alpha = s.dot(z) / z.squaredNorm() ;
			alpha = std::min( a_max, std::max(a_min, alpha ) );

		}

		x = x_best ;
		return min_res ;
	}

} ;


} //ns pg_impl




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
	case projected_gradient::Descent:
		return solve< projected_gradient::Descent, NSLaw, RhsT, ResT >( law, b, x ) ;
	case projected_gradient::Conjugated:
		return solve< projected_gradient::Conjugated, NSLaw, RhsT, ResT >( law, b, x ) ;
	case projected_gradient::APGD:
		return solve< projected_gradient::APGD, NSLaw, RhsT, ResT >( law, b, x ) ;
	case projected_gradient::SPG:
		return solve< projected_gradient::SPG, NSLaw, RhsT, ResT >( law, b, x ) ;
	}

	return -1 ;
}

template < typename BlockMatrixType >
template < projected_gradient::Variant variant, typename NSLaw, typename RhsT, typename ResT >
typename ProjectedGradient< BlockMatrixType >::Scalar
ProjectedGradient< BlockMatrixType >::solve(
		const NSLaw &law, const RhsT &b, ResT &x ) const
{
	return pg_impl::PgMethod< variant >::solve( *this, law, b, x ) ;
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
