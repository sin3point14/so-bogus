/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef BOGUS_BLOCK_TRAITS_HPP
#define BOGUS_BLOCK_TRAITS_HPP

#include <vector>

namespace bogus {

template< typename Derived >
struct BlockMatrixTraits { } ;

template< typename BlockType >
struct BlockContainerTraits {
	typedef std::vector< BlockType > Type ;
} ;

template< typename BlockType >
struct BlockTraits
{
   typedef typename BlockType::Scalar Scalar ;
   enum { RowsAtCompileTime = BlockType::RowsAtCompileTime,
		  ColsAtCompileTime = BlockType::ColsAtCompileTime }  ;
} ;

// Transpose and matrix/vector product return types
// Specialization of these structures should define a ReturnType if the operation is allowed

template< typename BlockT >
struct BlockTransposeTraits {} ;

template< typename BlockT >
struct SelfTransposeTraits {} ;

template< typename BlockT >
struct BlockVectorProductTraits {} ;

template< typename LhsBlockT, typename RhsBlockT, bool TransposeLhs, bool TransposeRhs >
struct BlockBlockProductTraits {
	typedef LhsBlockT ReturnType ;
} ;

} // namespace bogus

#endif
