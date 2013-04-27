/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver          
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>                       
 *
 * This Source Code Form is subject to the terms of the Mozilla Public 
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef BOGUS_PRECONDITIONERS_HPP
#define BOGUS_PRECONDITIONERS_HPP

#include "../Block/BlockMatrix.hpp"

namespace bogus {

template < typename MatrixType >
class TrivialPreconditioner
{
public:
	explicit TrivialPreconditioner( const MatrixType & )
	{}

	template < bool transpose, typename ResT, typename RhsT >
	void apply( const RhsT& rhs, ResT &res ) const
	{
		res = rhs ;
	}
} ;

template < typename MatrixType >
class DiagonalPreconditioner
{
} ;

template < typename MatrixType >
class DiagonalLUPreconditioner
{
} ;

template < typename MatrixType >
class DiagonalLDLTPreconditioner
{
} ;

}

#endif

