/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef BOGUS_KRYLOV_IMPL_HPP
#define BOGUS_KRYLOV_IMPL_HPP

#include "Krylov.hpp"
#include "Preconditioners.impl.hpp"
#include "BlockSolverBase.impl.hpp"

#include "../Utils/NumTraits.hpp"

namespace bogus {



template < typename BlockMatrixType, template< typename BlockMatrixT > class PreconditionerType >
Krylov< BlockMatrixType, PreconditionerType >::Krylov(
		const BlockMatrixBase< BlockMatrixType > & matrix )
	: Base( &matrix, matrix.rows(), NumTraits< Scalar >::epsilon() )
{
	setMatrix( matrix ) ;
}

template < typename BlockMatrixType, template< typename BlockMatrixT > class PreconditionerType >
Krylov< BlockMatrixType, PreconditionerType >::Krylov()
	: Base( NULL, 100, NumTraits< Scalar >::epsilon() ), m_scale(0)
{
}

template < typename BlockMatrixType, template< typename BlockMatrixT > class PreconditionerType >
void Krylov< BlockMatrixType, PreconditionerType >::setMatrix(
		const BlockMatrixBase< BlockMatrixType > & matrix )
{
	m_matrix = &matrix ;
	m_preconditioner.setMatrix( matrix ) ;
	m_scale = 1. / ( 1 + m_matrix->size() ) ;

	return *this ;
}

template < typename BlockMatrixType, template< typename BlockMatrixT > class PreconditionerType >
template < typename RhsT, typename ResT >
typename Krylov< BlockMatrixType, PreconditionerType >::Scalar
Krylov< BlockMatrixType, PreconditionerType >::init( const RhsT &b, ResT &x,
																	typename GlobalProblemTraits::DynVector &r0  ) const
{
	r0 = b - (*m_matrix)*x ;

	Scalar res = r0.squaredNorm() ;
	const Scalar resAt0 = b.squaredNorm()  ;

	if( res > resAt0 ) {
		r0 = b;
		x.setZero() ;
		res = resAt0 ;
	}
	return res * m_scale ;
}



template < typename BlockMatrixType, template< typename BlockMatrixT > class PreconditionerType >
template < typename RhsT, typename ResT >
typename Krylov< BlockMatrixType, PreconditionerType >::Scalar
Krylov< BlockMatrixType, PreconditionerType >::solve( const RhsT &b, ResT &x,
																	 krylov::Method method ) const
{
	switch(method)
	{
		case krylov::CG:
			return solve_CG( b, x ) ;
		case krylov::BiCG:
			return solve_BiCG( b, x ) ;
		case krylov::BiCGSTAB:
			return solve_BiCGSTAB( b, x ) ;
		case krylov::CGS:
			return solve_CGS( b, x ) ;
		case krylov::GMRES:
			return solve_GMRES( b, x ) ;
		case krylov::TFQMR:
			return solve_TFQMR( b, x ) ;
	}
	return -1 ;
}

} // namespace bogus

#endif
