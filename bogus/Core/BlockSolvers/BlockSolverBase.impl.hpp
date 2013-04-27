/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver          
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>                       
 *
 * This Source Code Form is subject to the terms of the Mozilla Public 
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef BOGUS_BLOCK_SOLVER_BASE_IMPL_HPP
#define BOGUS_BLOCK_SOLVER_BASE_IMPL_HPP

#include "BlockSolverBase.hpp"

namespace bogus {

template < typename BlockMatrixType >
BlockSolverBase< BlockMatrixType >::BlockSolverBase(
		const BlockMatrixBase< BlockMatrixType > & M,
		unsigned maxIters, Scalar tol )
	: m_matrix( M ), m_maxIters( maxIters ), m_tol( tol )
{
}

} // namespace bogus

#endif
