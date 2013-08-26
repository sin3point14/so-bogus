/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef BOGUS_KRYLOV_HPP
#define BOGUS_KRYLOV_HPP

#include "BlockSolverBase.hpp"
#include "Preconditioners.hpp"

#include <vector>

namespace bogus
{


//! Preconditionned Krylov Solvers
/*!
  \tparam BlockMatrixT The type of system matrix, which should be a subclass of BlockMatrixBase
  \tparam PreconditionerType The preconditioner type. It should accept BlockMatrixT as a template parameter.
	The default value, TrivialPreconditioner, means that no preconditioning will be done.
	\sa TrivialPreconditioner, DiagonalPreconditioner, DiagonalLUPreconditioner, DiagonalLDLTPreconditioner
	\sa krylov
  */

template < typename BlockMatrixType,
		   template< typename BlockMatrixT > class PreconditionerType >
class Krylov : public BlockSolverBase< BlockMatrixType >
{
public:
	typedef BlockSolverBase< BlockMatrixType > Base ;

	typedef typename Base::LocalMatrixType LocalMatrixType ;
	typedef typename Base::GlobalProblemTraits GlobalProblemTraits ;
	typedef typename GlobalProblemTraits::Scalar Scalar ;

	//! Default constructor -- you will have to call setMatrix() before using any of the solve() functions
    Krylov( ) ;
	//! Constructor with the system matrix -- initializes preconditioner
    explicit Krylov( const BlockMatrixBase< BlockMatrixType > & matrix ) ;

	//! Sets the system matrix and initializes the preconditioner
	void setMatrix( const BlockMatrixBase< BlockMatrixType > & matrix ) ;

	//! Solves ( m_matrix * \p x = \p b ) using the Conjugate Gradient algorithm
    /*! For symmetric matrices only. Converges for positive definite linear systems.

        <b>Matrix-vector mults/iter: </b> 1
        <b>Preconditionner calls/iter: </b> 1
        <b>Storage requirements: </b> 4n
    */
	template < typename RhsT, typename ResT >
	Scalar solve_CG( const RhsT &b, ResT &x ) const ;

	//! Solves ( m_matrix * \p x = \p b ) using the BiConjugate Gradient algorithm
	/*! Works for non-symmetric linear systems. Convergence not guaranteed.
        Requires ability to perform transpose multiplication and preconditioning

        <b>Matrix-vector mults/iter: </b> 2 ( inc. 1 transpose )
        <b>Preconditionner calls/iter: </b> 2 ( inc. 1 transpose )
        <b>Storage requirements: </b> 8n
    */
	template < typename RhsT, typename ResT >
	Scalar solve_BiCG( const RhsT &b, ResT &x ) const ;

	//! Solves ( m_matrix * \p x = \p b ) using the BiConjugate Gradient stabilized algorithm
    /*! Works for non-symmetric linear systems. Convergence not guaranteed.
        Supposedly less erratic convergence than BiCG method.

        <b>Matrix-vector mults/iter: </b> 2
        <b>Preconditionner calls/iter: </b> 2
        <b>Storage requirements: </b> 8n
    */
	template < typename RhsT, typename ResT >
	Scalar solve_BiCGSTAB( const RhsT &b, ResT &x ) const ;

	//! Solves ( m_matrix * \p x = \p b ) using the Conjugate Gradient Squared algorithm
    /*! Works for non-symmetric linear systems. Convergence not guaranteed.
        Supposedly faster convergence than BiCG \a when \a converging.

        <b>Matrix-vector mults/iter: </b> 2
        <b>Preconditionner calls/iter: </b> 2
        <b>Storage requirements: </b> 7n
    */
	template < typename RhsT, typename ResT >
	Scalar solve_CGS( const RhsT &b, ResT &x ) const ;

    //! Solves ( m_matrix * \p x = \p b ) using the (restarted) Generalized Minimum Residual
    /*!
        \param restart If non-zero, use the GMRES(m) restarted algorithm. Lower the storage cost,
        but slows-down or even forbid the convergence of the algorithm.

        Works for non-symmetric linear systems.
        Probably the more robust method for non symmetric systems, but with the highest storage
        cost.

        <b>Matrix-vector mults/iter: </b> 1
        <b>Preconditionner calls/iter: </b> 1
        <b>Other ops/iter: </b> 1 k*k triangular solve, 2 k*n m-v mult, 1 k*k m-v mult
        <b>Storage requirements: </b> 2*m*( n + m )

    */
    template < typename RhsT, typename ResT >
	Scalar solve_GMRES( const RhsT &b, ResT &x, unsigned restart = 0 ) const ;

	//! Solves ( m_matrix * \p x = \p b ) using the transpose-free Quasi Minimal Reisual method
    /*! Works for non-symmetric linear systems. Convergence not guaranteed.
        Supposedly less erratic convergence than BiCG method, but faster than BiCGSTAB.

        <b>Matrix-vector mults/iter: </b> 2
        <b>Preconditionner calls/iter: </b> 2
        <b>Storage requirements: </b> 7n

        \warning This function returns an approximation of the residual instead of the real one
    */
	template < typename RhsT, typename ResT >
	Scalar solve_TFQMR( const RhsT &b, ResT &x ) const ;

	//! Solve function that takes the method to use as an argument
	template < typename RhsT, typename ResT >
	Scalar solve(  const RhsT &b, ResT &x,
				   krylov::Method method = krylov::CG ) const ;

protected:
	using Base::m_matrix ;
	using Base::m_maxIters ;
	using Base::m_tol ;

	//! Check init guess, reset it to zero if that would give a lower residual
	template < typename RhsT, typename ResT >
	Scalar init( const RhsT &b, ResT &x, typename GlobalProblemTraits::DynVector &r0 ) const ;

	PreconditionerType< BlockMatrixBase< BlockMatrixType > > m_preconditioner ;
	Scalar m_scale ;
} ;

} //namesoace bogus

#endif

