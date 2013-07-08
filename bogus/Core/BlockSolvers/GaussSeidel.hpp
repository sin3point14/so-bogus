/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_BLOCK_GAUSS_SEIDEL_HPP
#define BOGUS_BLOCK_GAUSS_SEIDEL_HPP

#include "BlockSolverBase.hpp"

#include <vector>

namespace bogus
{

//! Projected Gauss-Seidel iterative solver.
/*!
   Works by taking into account only one block-row of the system at a time, and iterating
   several times over the whole set of rows several times until convergence has been achieved.

   Each inner iteration of the algorithm will try to solve the local problem
	  \f[
		\left\{
		  \begin{array}{rcl}
			y_i^{k+1} &=& M_{i,i} x_i^{k+1}  + b_i^{k} \\
			&s.t.& law (x^{k+1},y^{k+1})
		  \end{array}
		\right.
	  \f]
	where \b k is the current global iteration, \b i the current row
	and \f[ b_i^{k} := b_i + \sum_{ j < i }{ M_{i,j}x_j^{k+1} } +  \sum_{ j > i }{ M_{i,j}x_j^{k} } \f]

   See also solve() and \cite JAJ98.
  */
template < typename BlockMatrixType >
class GaussSeidel : public BlockSolverBase< BlockMatrixType >
{
public:
	typedef BlockSolverBase< BlockMatrixType > Base ;

	typedef typename Base::GlobalProblemTraits GlobalProblemTraits ;
	typedef typename GlobalProblemTraits::Scalar Scalar ;

	//! Default constructor -- you will have to call setMatrix() before using the solve() function
	GaussSeidel( ) ;
	//! Constructor with the system matrix
	explicit GaussSeidel( const BlockMatrixBase< BlockMatrixType > & matrix ) ;

	//! Sets the system matrix and initializes internal structures
	void setMatrix( const BlockMatrixBase< BlockMatrixType > & matrix ) ;

	//! Solves a constrained linear system
	/*!
	  Implements Algorithm 1. from \cite DBB11 to solve
	   \f[
		\left\{
		  \begin{array}{rcl}
			y &=& M x + b \\
			&s.t.& law (x,y)
		  \end{array}
		\right.
	  \f]
	  \param law The (non-smooth) law that should define:
		- An error function for the global problem
		- A local solver for each row of the system ( e.g. 1 contact solver )
		\sa SOCLaw
	  \param b the const part of the right hand side
	  \param x the unknown. Can be warm-started
	  */
	template < typename NSLaw, typename RhsT, typename ResT >
	Scalar solve( const NSLaw &law, const RhsT &b, ResT &x ) const ;

	//! Sets whether the solver is allowed to trade off determiniticity for speed
	void setDeterministic( bool deterministic ) { m_deterministic = deterministic ; }

	//! Sets the auto-regularization coefficient
	/*!
	  The regularization works by slightly altering the local problems, so at each iteration
	  we try to solve
	  \f[
		\left\{
		  \begin{array}{rcl}
			y^{k+1} &=& ( M + \alpha I ) x^{k+1} - \alpha x^k + b^{k} \\
			&s.t.& law (x^{k+1},y^{k+1})
		  \end{array}
		\right.
	  \f]
	  where \f$\alpha\f$ is the regularization coefficient.

	  Note that as \f$ | x^{k+1} - x^{k} | \rightarrow 0 \f$ when the algorithm converges, we are still
	  trying to find a solution of the same global problem.

	  For under-determined problems, regularization might helps preventing \b x reaching problematically high values.
	  Setting \f$\alpha\f$ to a too big value will however degrade the convergence of the global algorithm.

	  \param maxRegul If greater than zero, then positive terms will be added to the diagonal
	  of the local matrices so that their minimal eigenvalue become greater than \p maxRegul.
	  */
	void setAutoRegularization( Scalar maxRegul ) { m_autoRegularization = maxRegul ; }

	// Debug

	//! Sets the number of iterations that should be performed between successive evaluations of the global error function
	/*! ( Those evaluations require a full matrix/vector product, and are therfore quite costly ) */
	void setEvalEvery( unsigned evalEvery ) { m_evalEvery = evalEvery ; }
	//! Sets the minimum iteration step size under which local problems are temporarily frozen
	void setSkipTol  ( unsigned skipTol   ) { m_skipTol   = skipTol   ; }
	//! Sets the number of iterations for temporarily freezing local problems
	void setSkipIters( unsigned skipIters ) { m_skipIters = skipIters ; }

protected:

	//! Eval the current global reisual
	/*! \p y should be such that \p y = m_matrix * \p x */
	template < typename NSLaw, typename RhsT, typename ResT >
	Scalar eval ( const NSLaw &law, const RhsT &x, const ResT &y ) const ;

	using Base::m_matrix ;
	using Base::m_maxIters ;
	using Base::m_tol ;

	typedef typename LocalProblemTraits< GlobalProblemTraits::dimension, Scalar >::Matrix DiagonalMatrixType ;
	typename BlockContainerTraits< DiagonalMatrixType >::Type m_localMatrices ;
	typename GlobalProblemTraits::DynVector m_scaling ;
	typename GlobalProblemTraits::DynVector m_regularization ;

	bool m_deterministic ;

	unsigned m_evalEvery ;
	Scalar m_skipTol ;
	unsigned m_skipIters ;
	Scalar m_autoRegularization ;
} ;

} //namespace bogus


#endif
