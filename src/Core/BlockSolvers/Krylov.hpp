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
#include "KrylovMethods.hpp"

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
	Krylov& setMatrix( const BlockMatrixBase< BlockMatrixType > & matrix ) ;

	// For each value of the krylov::Method enum, create a MethodName typedef
	// and the asMethodName() -> MethodName and solve_MethodName -> Scalar methods
#define BOGUS_PROCESS_KRYLOV_METHOD( MethodName )\
	typedef krylov::MethodName<  \
		BlockMatrixType, PreconditionerType< BlockMatrixBase< BlockMatrixType > >, \
		GlobalProblemTraits > MethodName ; \
	\
	MethodName as##MethodName() const { \
		return MethodName(  \
			Base::m_matrix->derived(), m_preconditioner, \
			Base::m_callback, \
			Base::m_tol, Base::m_maxIters ) ; } \
	\
	template < typename RhsT, typename ResT > \
	Scalar solve_##MethodName( const RhsT &b, ResT &x ) const  \
	{ return as##MethodName().solve( b, x ) ; }

BOGUS_KRYLOV_METHODS
#undef BOGUS_PROCESS_KRYLOV_METHOD

	//! Solve function that takes the method to use as an argument
	template < typename RhsT, typename ResT >
	Scalar solve(  const RhsT &b, ResT &x,
				   krylov::Method method = krylov::kCG ) const ;

protected:

	PreconditionerType< BlockMatrixBase< BlockMatrixType > > m_preconditioner ;
} ;

} //namesoace bogus

#endif

