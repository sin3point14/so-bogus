/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_ADMM_IMPL_HPP
#define BOGUS_ADMM_IMPL_HPP


#include "ADMM.hpp"

#include "ConstrainedSolverBase.impl.hpp"
#include "ProjectedGradient.impl.hpp"

#include "../Utils/NumTraits.hpp"

namespace bogus
{

//! Evaluation of prox_{1/l} (J) with J = 1/2 xM'x + f'x
/*! Requires way to solve linear system with lhs ( M + l I)
 */
template< typename ObjectType >
struct QuadraticProxOp
{
	// Usual defs
	typedef typename BlockMatrixTraits< typename ObjectType::PlainObjectType >::BlockType LocalMatrixType ;
	typedef ProblemTraits< LocalMatrixType > GlobalProblemTraits ;
	typedef typename GlobalProblemTraits::Scalar Scalar ;

	typedef BlockObjectBase< ObjectType > LinearOp ;
	typedef typename GlobalProblemTraits::DynVector AffineVec ;

	//! Mutliplying by linearExpr should be equivalent to solving with lhs ( M + lambda I )
	QuadraticProxOp( const LinearOp& linearExpr, const Scalar lambda,
					 const AffineVec& affinePart )
		: m_coefficient( lambda ), m_linearOp( linearExpr ), m_affinePart( affinePart )
	{
		//assert( ( m_coefficient > 0 ) ^ ( ForAMA ) ) ;
	}

	//! Evaluates prox_{J, 1/coeff}( rhs/coeff )
	/*! \warning may modify rhs */
	//  =  min J(x) + coeff/{2} || x  - tmp/coeff ||^2
	//  =  min J(x) + coeff/{2} || x ||^2 - <tmp, x>
	//  =  min 1/2 ( x'(M + coeff I)x  ) + <f - tmp, x>
	//  = (M + coeff I)^{-1} ( tmp - f )
	template < typename RhsT, typename ResT >
	void eval( RhsT & rhs, ResT & res ) const {
		rhs -= m_affinePart ;
		m_linearOp.template multiply< false >( rhs, res, 1, 0 ) ;
	}

	Scalar coefficient() const {
		return m_coefficient ;
	}

private:
	const Scalar     m_coefficient ;
	const LinearOp&  m_linearOp ;
	const AffineVec& m_affinePart ;
};

template < typename BlockMatrixType >
template < typename NSLaw, typename ProxOp, typename RhsT, typename ResT >
typename ADMM< BlockMatrixType >::Scalar
ADMM< BlockMatrixType >::solve(
		const NSLaw &law, const ProxOp& op,
		const RhsT &w, ResT &v, ResT &r ) const
{
	switch ( m_defaultVariant )
	{
	case admm::Standard:
		return solve< admm::Standard, NSLaw, ProxOp, RhsT, ResT >( law, op, w, v, r ) ;
	case admm::Accelerated:
		return solve< admm::Accelerated, NSLaw, ProxOp, RhsT, ResT >( law, op, w, v, r ) ;
	}

	return -1 ;
}

template < typename BlockMatrixType >
template < admm::Variant variant, typename NSLaw, typename ProxOp, typename RhsT, typename ResT >
typename ADMM< BlockMatrixType >::Scalar
ADMM< BlockMatrixType >::solve(
		const NSLaw &law, const ProxOp& op,
		const RhsT &w, ResT &x, ResT &r ) const
{

	const Scalar lambda = op.coefficient() ;
	const Scalar gamma  = stepSize() ; // gamma/(lambda*lambda)
	const Scalar inv_gamma = 1./ gamma ;

	Scalar res = -1 ;

	typename GlobalProblemTraits::DynVector ut, prox_arg( x.rows() ), z ;

	ut = w ;
	m_matrix->template multiply<false>( x, ut, 1, 1 ) ;
	z  = ut - inv_gamma * r ;
	this->projectOnConstraints( law, z ) ;

	// Nesterov acceleration
	typename GlobalProblemTraits::DynVector y, y_prev = - inv_gamma * r ;
	typename GlobalProblemTraits::DynVector zacc = z, z_prev = z ;
	Scalar theta_prev = 1. ; // Previous Nesterov acceleration
	Scalar res_prev = -1 ;

	for( unsigned adIter = 0 ; adIter < m_maxIters ; ++ adIter ) {

		if( lambda < NumTraits< Scalar >::epsilon() ) {
			// AMA
			m_matrix->template multiply<true>( r, prox_arg, 1, 0 ) ;
		} else {
			// ADMM
			prox_arg = x ;
			ut = r + gamma*(zacc - ut) ; //Re-use ut as we no-longer need its current value
			m_matrix->template multiply<true>( ut, prox_arg, 1, 1./lambda ) ;
		}
		op.eval( prox_arg, x ) ;

		ut = w ;
		m_matrix->template multiply<false>( x, ut, 1, 1 ) ;

		res = this->eval( law, r, ut )
				+ (	this->usesInfinityNorm()
					? (ut-z).template lpNorm< Eigen::Infinity >()
					: (ut-z).squaredNorm() / (1 + ut.rows() ) ) ;
				;
		this->callback().trigger( adIter, res );

		if( res < this->tol() )
			break ;


		z  = ut - inv_gamma * r ;
		this->projectOnConstraints( law, z ) ;

		if( variant == admm::Accelerated && res < res_prev ) {

			y = ut - z - inv_gamma * r ;
			const Scalar beta = bogus::pg_impl::nesterov_inertia( theta_prev, 0. ) ;

			r = - gamma * ( y + beta * ( y - y_prev ) ) ; // Over-relaxation
			y_prev = y ;

			zacc   = z + beta * (z - z_prev) ;
			z_prev = z ;

		} else {
			r += gamma * ( z - ut) ;

			zacc = z ; z_prev = z ;
			y_prev = - inv_gamma * r ;
			theta_prev = 1 ;
		}

//		const Scalar g = (ut - z).squaredNorm() ;
//		std::cout << adIter << " \t gap: " << g << " \t res " << res << std::endl ;

		res_prev = res ;

	}

	return res ;

}

template < typename BlockMatrixType >
template < typename NSLaw, typename MatrixT, typename RhsT, typename ResT >
typename DualAMA< BlockMatrixType >::Scalar
DualAMA< BlockMatrixType>::solve(
		const NSLaw &law, const BlockObjectBase< MatrixT >& A,
		const RhsT &f, const RhsT &w, ResT &v, ResT &r ) const
{
	switch( m_defaultVariant ) {
		case admm::Accelerated:
			return solve<admm::Accelerated>( law, A, f, w, v, r) ;
		case admm::Standard:
			return solve<admm::Standard>( law, A, f, w, v, r) ;
	}

	return -1 ;
}

template < typename BlockMatrixType >
template < admm::Variant variant, typename NSLaw, typename MatrixT, typename RhsT, typename ResT >
typename DualAMA< BlockMatrixType >::Scalar
DualAMA< BlockMatrixType>::solve(
		const NSLaw &law, const BlockObjectBase< MatrixT >& A,
		const RhsT &f, const RhsT &w, ResT &v, ResT &r ) const
{
	typedef typename GlobalProblemTraits::DynVector DynVec ;

	Scalar lambda = projStepSize() ;
	const Scalar gamma  = fpStepSize() ;

	Scalar res = -1 ;

	typename GlobalProblemTraits::DynVector x, ut, gap, Hr ( v.rows() ), g, s( r.rows() ) ;

	const Segmenter< NSLaw::dimension, const DynVec, typename BlockMatrixType::Index >
			 utSegmenter( ut, m_matrix->rowOffsets() ) ;
	Segmenter< NSLaw::dimension, DynVec, typename BlockMatrixType::Index >
			 sSegmenter( s, m_matrix->rowOffsets() ) ;

	// Nesterov acceleration
	DynVec y, y_prev = v ;
	Scalar theta_prev = 1. ; // Previous Nesterov acceleration
	Scalar res_prev = -1 ;

	// Line search
	DynVec rs ;

	m_matrix->template multiply< true >( r, Hr, 1, 0 ) ;

	for( unsigned adIter = 0 ; adIter < m_maxIters ; ++ adIter )
	{

		x  = f ;
		A.template multiply< false >( v, x, 1, 1 ) ;
		gap = Hr - x ;

		ut = w ;
		m_matrix->template multiply<false>( v, ut, 1, 1 ) ;

		// Eval current reisual,  exit if small enough
		res = this->eval( law, ut, r ) +
				( this->usesInfinityNorm()
				  ? gap.template lpNorm< Eigen::Infinity >()
				  : gap.squaredNorm() / (1 + gap.rows() ) ) ;

		this->callback().trigger( adIter, res );

		if( res < this->tol() )
			break ;


		// Acceleration
		Scalar beta = 0 ;
		if ( variant == admm::Accelerated && res < res_prev ) {
			beta =  bogus::pg_impl::nesterov_inertia( theta_prev, 0. ) ;
		} else {
			theta_prev = 1 ;
		}

		this->dualityCOV( law, ut, s ) ;
		g = ut + s ;
		m_matrix->template multiply< false >( gap, g, gamma, 1 ) ;

		// (Optional) Line search
		if( lineSearchIterations() == 0 ) {
			r -= lambda * g ;
			this->projectOnConstraints( law, r ) ;

			m_matrix->template multiply< true >( r, Hr, 1, 0 ) ;

		} else {

			// TODO avoid performing extra A multiplication by keeping LS results

			const Scalar h0 = gap.squaredNorm() ;


			lambda *= lineSearchOptimisticFactor() ;
			for( unsigned lsIter = 0 ; lsIter < lineSearchIterations() ; ++lsIter ) {

				rs = r - lambda * g ;
				this->projectOnConstraints( law, rs ) ;
				m_matrix->template multiply< true >( rs, Hr, 1, 0 ) ;

				y = v + gamma * ( Hr - x ) ;
				y += beta * ( y - y_prev ) ; //Over Relaxation


				gap = Hr - f ;
				A.template multiply< false >( y, gap, -1, 1 ) ;

				Scalar h = gap.squaredNorm() ;

				if( h < h0 ) {
//					std::cout << lsIter << " [ "  << lambda<< " ] " << " \t " << h << " vs " << h0 << std::endl ;
					break ;
				}

				lambda *= lineSearchPessimisticFactor() ;
			}

			r = rs ;

		}

		// (end line-search)


		y = v + gamma * ( Hr - x ) ;
		v = y +  beta * ( y  - y_prev ) ;  //Over Relaxation

		y_prev   = y ;
		res_prev = res ;


	}

	return res ;

}


} //namespace bogus

#endif
