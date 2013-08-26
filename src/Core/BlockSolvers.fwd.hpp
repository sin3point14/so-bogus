/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_BLOCK_SOLVERS_FWD_HPP
#define BOGUS_BLOCK_SOLVERS_FWD_HPP

namespace bogus {

//! Configuration properties for Krylov solvers
namespace krylov
{
enum Method
{
	CG,				//!< Conjugate Gradient. \sa Krylov::solve_CG()
	BiCG,			//!< BiConjugate Gradient \sa Krylov::solve_BiCG()
	BiCGSTAB, 		//!< BiConjugate Gradient Stabilized \sa Krylov::solve_BiCGSTAB()
	CGS, 			//!< Conjugate Gradient Squared \sa Krylov::solve_CGS()
	GMRES,			//!< Generalized Minimal Residual \sa Krylov::solve_GMRES()
	TFQMR			//!< Tranpose-free Quasi Minimal Residual \sa Krylov::solve_TFQMR()
} ;

} // namespace iterative_linear_solvers

template < typename MatrixType >
struct ProblemTraits ;

template < typename BlockMatrixType >
class GaussSeidel ;

template < typename MatrixType >
class TrivialPreconditioner ;

template < typename BlockMatrixType,
		   template< typename BlockMatrixT > class PreconditionerType = TrivialPreconditioner >
class Krylov ;

}

#endif
