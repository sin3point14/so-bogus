/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_LINEAR_SOLVER_HPP
#define BOGUS_LINEAR_SOLVER_HPP

#include <cassert>

#include "../Block.fwd.hpp"

namespace bogus {

template < typename LSDerived >
struct LinearSolverTraits {} ;

//! Base class for linear solvers on base ( i.e. non-block ) matrices
template < typename Derived >
struct LinearSolverBase
{

	//! Returns the solution \b x of the linear system \b M \c * \b x \c = \c rhs
	template < typename RhsT >
	typename LinearSolverTraits< Derived >::template Result< RhsT >::Type
	solve( const RhsT& rhs ) const
	{
	   return derived().solve( rhs ) ;
	}

	//! Finds the solution \b x of the linear system \b M \c * \b x \c = \c b
	template < typename ResT, typename RhsT >
	void solve( const RhsT& rhs, ResT& x ) const
	{
	   return derived().solve( rhs, x ) ;
	}

	const Derived& derived() const
	{
		return static_cast< const Derived& >( *this ) ;
	}

	Derived& derived() {
		return static_cast< Derived& >( *this ) ;
	}


	typedef BlockTraits< typename LinearSolverTraits< Derived >::MatrixType > UnderlyingBlockTraits ;
	typedef typename UnderlyingBlockTraits::Scalar Scalar ;
	enum {
		RowsAtCompileTime = UnderlyingBlockTraits::RowsAtCompileTime,
		ColsAtCompileTime = UnderlyingBlockTraits::ColsAtCompileTime,
		IsRowMajor = UnderlyingBlockTraits::is_row_major
	} ;

	//! Tag: mv_add() and mv_set() are redefined for LinearSolverBase< Derived >
	typedef char MVOverload ;
	//! Tag: mm_add() and mm_set() are redefined for LinearSolverBase< Derived >
	typedef char MMOverload ;
} ;

//! Base class for LU factorizations
template < typename MatrixType >
struct LU : public LinearSolverBase< LU< MatrixType > >
{
} ;

//! Base class for LDLT factorizations
template < typename MatrixType >
struct LDLT : public LinearSolverBase< LDLT< MatrixType > >
{
} ;


//! Block product type deductions

template< typename LhsBlockT, typename RhsBlockT, bool TransposeLhs, bool TransposeRhs >
struct BlockBlockProductTraits < LU < LhsBlockT >, RhsBlockT, TransposeLhs, TransposeRhs >
{
	typedef typename LinearSolverTraits< LU < LhsBlockT > >::MatrixType LhsMatrixT ;
	typedef typename BlockBlockProductTraits < LhsMatrixT, RhsBlockT, TransposeLhs, TransposeRhs >::ReturnType
	ReturnType ;
} ;

template< typename LhsBlockT, typename RhsBlockT, bool TransposeLhs, bool TransposeRhs >
struct BlockBlockProductTraits < LDLT < LhsBlockT >, RhsBlockT, TransposeLhs, TransposeRhs >
{
	typedef typename LinearSolverTraits< LDLT < LhsBlockT > >::MatrixType LhsMatrixT ;
	typedef typename BlockBlockProductTraits < LhsMatrixT, RhsBlockT, TransposeLhs, TransposeRhs >::ReturnType
	ReturnType ;
} ;

}// namespace bogus

#endif
