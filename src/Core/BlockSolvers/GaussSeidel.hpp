/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver          
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>                       
 *
 * This Source Code Form is subject to the terms of the Mozilla Public 
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_BLOCK_GAUSS_SEIDEL_HPP
#define BOGUS_BLOCK_GAUSS_SEIDEL_HPP

#include "../Block.fwd.hpp"
#include "../BlockSolvers.fwd.hpp"

namespace bogus
{

template < typename BlockMatrixType >
class GaussSeidel
{
public:

	typedef typename BlockMatrixTraits< BlockMatrixType >::BlockType LocalMatrixType ;
	typedef ProblemTraits< LocalMatrixType > ProblemTraits ;
	typedef typename ProblemTraits::Scalar Scalar ;

	explicit GaussSeidel( const BlockMatrixBase< BlockMatrixType > & M ) ;

	template < typename NSLaw, typename RhsT, typename ResT >
	Scalar solve( const NSLaw &law, const RhsT &b, ResT &x ) const ;

	void setMaxIters( unsigned maxIters ) { m_maxIters = maxIters ; }
	void setTol( double tol ) { m_tol = tol ; }
	void setDeterministic( bool deterministic ) { m_deterministic = deterministic ; }

private:
	const BlockMatrixBase< BlockMatrixType > & m_matrix ;
	std::vector< LocalMatrixType > m_localMatrices ;
	typename ProblemTraits::DynVector m_scaling ;

	unsigned m_maxIters ;
	Scalar m_tol ;
	bool m_deterministic ;

	unsigned m_evalEvery ;
	Scalar m_skipTol ;
	Scalar m_skipIters ;
} ;

} //namespace bogus


#endif
