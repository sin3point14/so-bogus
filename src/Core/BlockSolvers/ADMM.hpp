/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_ADMM_HPP
#define BOGUS_ADMM_HPP

#include "ConstrainedSolverBase.hpp"

#include <vector>

namespace bogus
{

//! Options for ADMM solvers
namespace admm {
	//! Variants of ADMM algorithm
	enum Variant {
		//! Standard ADMM
		Standard,
		//! ADDM with \cite Nesterov1983 acceleration (see Goldstein 2014)
		Accelerated
	} ;
}

//! ADMM iterative solver.
/*!
	Minimizes J(x) with ( M x + b ) in C

	Requires ability to evaluate prox_J( x )
*/
template < typename BlockMatrixType >
class ADMM : public ConstrainedSolverBase< ADMM, BlockMatrixType >
{
public:
	typedef ConstrainedSolverBase< bogus::ADMM, BlockMatrixType > Base ;

	typedef typename Base::GlobalProblemTraits GlobalProblemTraits ;
	typedef typename GlobalProblemTraits::Scalar Scalar ;

	//! Default constructor -- you will have to call setMatrix() before using the solve() function
	ADMM( ) : Base() { init() ; }
	//! Constructor with the system matrix
	explicit ADMM( const BlockObjectBase< BlockMatrixType > & matrix ) : Base()
	{ init() ; Base::setMatrix( matrix ) ; }

	/*!
	 *  Solve J(x),  ( Mx + b \in C )
	 *
	 */
	template < admm::Variant variant, typename NSLaw, typename ProxOp, typename RhsT, typename ResT >
	Scalar solve(
			const NSLaw &law, const ProxOp& op,
			const RhsT &b, ResT &x, ResT &r ) const ;

	ADMM& setMatrix( const BlockObjectBase< BlockMatrixType > & matrix )
	{
		m_matrix = &matrix ;
		Base::updateScalings() ;
		return *this ;
	}

	//! Sets the step size for updating the dual variable (forces).
	void setStepSize( const Scalar size )
	{ m_stepSize = size ; }

	//! Sets the variant that will be used when calling solve() without template arguments
	void setDefaultVariant( admm::Variant variant )
	{ m_defaultVariant = variant ; }

	Scalar stepSize() const { return m_stepSize ; }

protected:

	typedef typename Base::Index Index ;

	//! Sets up the default values for all parameters
	void init()
	{
		m_tol = 1.e-6 ;
		m_maxIters = 300 ;

		m_stepSize = 1.e-2 ;

		m_defaultVariant = admm::Accelerated ;
	}

	using Base::m_matrix ;
	using Base::m_maxIters ;
	using Base::m_tol ;

	Scalar m_stepSize ;

	admm::Variant m_defaultVariant ;

} ;

} //namespace bogus


#endif

