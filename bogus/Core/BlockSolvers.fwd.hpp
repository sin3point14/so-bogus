/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_BLOCK_SOLVERS_FWD_HPP
#define BOGUS_BLOCK_SOLVERS_FWD_HPP

namespace bogus {

//! Configuration properties Iterative Linear Solvers
namespace iterative_linear_solvers
{
enum Method
{
	CG,				//!< Conjugate Gradient
	BiCG,			//!< BiConjugate Gradient
	BiCG_STAB, 		//!< BiConjugate Gradient Stabilized
	GMRES,			//!< Generalized Minimal Residual
	TFQMR			//!< Tranpose-free Quasi Minimum Residual
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
class IterativeLinearSolver ;

}

#endif
