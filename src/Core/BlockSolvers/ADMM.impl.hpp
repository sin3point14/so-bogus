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

//template < typename BlockMatrixType >
//template < typename NSLaw, typename ProxOp, typename RhsT, typename ResT >
//typename ADMM< BlockMatrixType >::Scalar
//ADMM< BlockMatrixType >::solve(
//		const NSLaw &law, const ProxOp& op,
//		const RhsT &b, const RhsT &f,
//		ResT &x, ResT &r ) const
//{
//	switch ( m_defaultVariant )
//	{
//	case admm::Standard:
//		return solve< admm::Standard, NSLaw, RhsT, ResT >( law, b, x ) ;
//	case admm::Fast_ADMM:
//		return solve< admm::Fast_ADMM, NSLaw, RhsT, ResT >( law, b, x ) ;
//	case admm::AMA:
//		return solve< admm::AMA, NSLaw, RhsT, ResT >( law, b, x ) ;
//	case admm::Fast_AMA:
//		return solve< admm::Fast_AMA, NSLaw, RhsT, ResT >( law, b, x ) ;
//	}

//	return -1 ;
//}

template < typename BlockMatrixType >
template < admm::Variant variant, typename NSLaw, typename ProxOp, typename RhsT, typename ResT >
typename ADMM< BlockMatrixType >::Scalar
ADMM< BlockMatrixType >::solve(
		const NSLaw &law, const ProxOp& op,
		const RhsT &b, ResT &x, ResT &r ) const
{

	const Scalar lambda = op.coefficient() ;
	const Scalar gamma  = stepSize() ; // gamma/(lambda*lambda)
	const Scalar inv_gamma = 1./ gamma ;

	Scalar res = -1 ;

	typename GlobalProblemTraits::DynVector ut, prox_arg( x.rows() ), z ;

	ut = b ;
	m_matrix->template multiply<false>( x, ut, 1, 1 ) ;
	z  = ut - inv_gamma * r ;
	this->projectOnConstraints( law, z ) ;

	// Nesterov acceleration
	typename GlobalProblemTraits::DynVector y, y_prev = - inv_gamma * r ;
	Scalar theta_prev = 1. ; // Previous Nesterov acceleration
	Scalar res_prev = -1 ;

	for( unsigned adIter = 0 ; adIter < m_maxIters ; ++ adIter ) {

		if( lambda < NumTraits< Scalar >::epsilon() ) {
			// AMA
			m_matrix->template multiply<true>( r, prox_arg, 1, 0 ) ;
		} else {
			// ADMM
			prox_arg = x ;
			ut = r + gamma*(z - ut) ; //Re-use ut as we no-longer need its current value
			m_matrix->template multiply<true>( ut, prox_arg, 1, 1./lambda ) ;
		}
		op.eval( prox_arg, x ) ;

		ut = b ;
		m_matrix->template multiply<false>( x, ut, 1, 1 ) ;


		res = this->eval( law, r, ut ) ;
		this->callback().trigger( adIter, res );

		if( res < this->tol() )
			break ;


		z  = ut - inv_gamma * r ;
		this->projectOnConstraints( law, z ) ;

		if( variant == admm::Accelerated ) {

			if( res_prev < res ) theta_prev = 1 ; //Restart

			y = ut - z - inv_gamma * r ;
			const Scalar beta = bogus::pg_impl::nesterov_inertia( theta_prev, 0. ) ;

			r = - gamma * ( y + beta * ( y - y_prev ) ) ; // Over-relaxation
			y_prev = y ;

		} else {
			r += gamma * ( z - ut) ;
		}

//		const Scalar g = (ut - z).squaredNorm() ;
//		std::cout << adIter << " \t gap: " << g << " \t res " << res << std::endl ;

		res_prev = res ;

	}

	return res ;

}

} //namespace bogus

#endif
