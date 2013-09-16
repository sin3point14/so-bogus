/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_BLOCK_GAUSS_SEIDEL_HPP
#define BOGUS_BLOCK_GAUSS_SEIDEL_HPP

#include "ConstrainedSolverBase.hpp"
#include "Coloring.hpp"

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
class GaussSeidel : public ConstrainedSolverBase< GaussSeidel, BlockMatrixType >
{
public:
	typedef ConstrainedSolverBase< bogus::GaussSeidel, BlockMatrixType > Base ;

	typedef typename Base::GlobalProblemTraits GlobalProblemTraits ;
	typedef typename GlobalProblemTraits::Scalar Scalar ;

	//! Default constructor -- you will have to call setMatrix() before using the solve() function
	GaussSeidel( ) : Base() { init() ; }
	//! Constructor with the system matrix
	explicit GaussSeidel( const BlockMatrixBase< BlockMatrixType > & matrix ) : Base()
	{ init() ; setMatrix( matrix ) ; }

	//! Sets the system matrix and initializes internal structures
	void setMatrix( const BlockMatrixBase< BlockMatrixType > & matrix ) ;

	//! Finds an approximate solution for a constrained linear problem
	/*!
	  Stops when the residual computed in eval() is below \ref m_tol, of the number
	  of iterations exceed \ref m_maxIters

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
		- An error function for the local problem
		- A local solver for each row of the system ( e.g. 1 contact solver )
		\sa SOCLaw
	  \param b the const part of the right hand side
	  \param x the unknown. Can be warm-started
	  \param tryZeroAsWell If true, the algorithm will reset r to zero if that would result in a lower residual
	  */
	template < typename NSLaw, typename RhsT, typename ResT >
	Scalar solve( const NSLaw &law, const RhsT &b, ResT &x, bool tryZeroAsWell = true ) const ;


	//! Access to the current Coloring. Will be reset whenever the matrix is changed.
	/*! Determiniticy is achieved through the mean of contact coloring ;
	  contacts that do not interact directly together can chare the same color,
	  and all contacts within a given color can be solver in parallel */
	Coloring& coloring( ) { return m_coloring ; }

	//! Sets the maximum number of threads that the solver can use.
	/*! If \p maxThreads is zero, then it will use the current OpenMP setting.

	\warning If multi-threading is enabled without coloring,
	  the result will not be deterministic, as it will depends on the
	  order in which threads solve contacts.

	  On the other hand, the algorithm will run much faster.
	*/
	void setMaxThreads( unsigned maxThreads = 0 ) {
		m_maxThreads = maxThreads ;
	}


	//! Sets the auto-regularization (a.k.a. proximal point) coefficient
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

	//! Sets up the default values for all parameters
	void init()
	{
		m_tol = 1.e-6 ;
		m_maxIters = 250 ;
		m_maxThreads =  0 ;
		m_evalEvery = 25  ;
		m_skipTol = m_tol * m_tol ;
		m_skipIters = 10 ;
		m_autoRegularization = 0. ;
	}

	using Base::m_matrix ;
	using Base::m_maxIters ;
	using Base::m_tol ;
	using Base::m_scaling ;

	typedef typename Base::Index Index ;
	typedef typename LocalProblemTraits< GlobalProblemTraits::dimension, Scalar >::Matrix DiagonalMatrixType ;
	typename ResizableSequenceContainer< DiagonalMatrixType >::Type m_localMatrices ;
	typename GlobalProblemTraits::DynVector m_regularization ;

	//! See setMaxThreads(). Defaults to 0 .
	unsigned m_maxThreads ;

	//! See setEvalEvery(). Defaults to 25
	unsigned m_evalEvery ;
	//! See setSkipTol(). Defaults to (\ref m_tol)²
	Scalar m_skipTol ;
	//! See setSkipIters() Defaults to 10
	unsigned m_skipIters ;

	//! \sa setAutoRegularization(). Defaults to 0.
	Scalar m_autoRegularization ;

	//! \sa coloring()
	Coloring m_coloring ;
} ;

} //namespace bogus


#endif
