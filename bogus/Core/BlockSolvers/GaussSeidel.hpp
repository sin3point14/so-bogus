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

//! Non-smooth Gauss-Seidel iterative solver
template < typename BlockMatrixType >
class GaussSeidel : public BlockSolverBase< BlockMatrixType >
{
public:
	typedef BlockSolverBase< BlockMatrixType > Base ;

	typedef typename Base::LocalMatrixType LocalMatrixType ;
	typedef typename Base::GlobalProblemTraits GlobalProblemTraits ;
	typedef typename GlobalProblemTraits::Scalar Scalar ;

	GaussSeidel( ) ;
	explicit GaussSeidel( const BlockMatrixBase< BlockMatrixType > & matrix ) ;

	void setMatrix( const BlockMatrixBase< BlockMatrixType > & matrix ) ;

	template < typename NSLaw, typename RhsT, typename ResT >
	Scalar solve( const NSLaw &law, const RhsT &b, ResT &x ) const ;

	void setDeterministic( bool deterministic ) { m_deterministic = deterministic ; }


	// Debug
	void setEvalEvery( unsigned evalEvery ) { m_evalEvery = evalEvery ; }
	void setSkipTol  ( unsigned skipTol   ) { m_skipTol   = skipTol   ; }
	void setSkipIters( unsigned skipIters ) { m_skipIters = skipIters ; }

protected:
	using Base::m_matrix ;
	using Base::m_maxIters ;
	using Base::m_tol ;

	typename BlockContainerTraits< LocalMatrixType >::Type m_localMatrices ;
	typename GlobalProblemTraits::DynVector m_scaling ;

	bool m_deterministic ;

	unsigned m_evalEvery ;
	Scalar m_skipTol ;
	Scalar m_skipIters ;
} ;

} //namespace bogus


#endif
