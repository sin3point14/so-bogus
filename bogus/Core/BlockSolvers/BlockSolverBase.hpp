/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef BOGUS_BLOCK_SOLVER_BASE_HPP
#define BOGUS_BLOCK_SOLVER_BASE_HPP

#include "../Block.fwd.hpp"
#include "../BlockSolvers.fwd.hpp"

namespace bogus
{

template < typename BlockMatrixType >
class BlockSolverBase
{
public:

	typedef typename BlockMatrixTraits< BlockMatrixType >::BlockType LocalMatrixType ;
	typedef ProblemTraits< LocalMatrixType > GlobalProblemTraits ;
	typedef typename GlobalProblemTraits::Scalar Scalar ;

	virtual ~BlockSolverBase() { }

	void setMaxIters( unsigned maxIters ) { m_maxIters = maxIters ; }
	void setTol( double tol ) { m_tol = tol ; }

protected:

	BlockSolverBase( const BlockMatrixBase< BlockMatrixType > & matrix,
	                 unsigned maxIters, Scalar tol ) ;

	const BlockMatrixBase< BlockMatrixType > & m_matrix ;

	unsigned m_maxIters ;
	Scalar m_tol ;
} ;

} //namespace bogus

#endif
