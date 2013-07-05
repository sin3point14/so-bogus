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
   enum {
	   //! Number of rows spanned by a block at compile time ; useful for efficient segmentation
	   RowsAtCompileTime = BlockType::RowsAtCompileTime,
	   //! Number of cols spanned by a block at compile time ; useful for efficient segmentation
	   ColsAtCompileTime = BlockType::ColsAtCompileTime,

	   //! Can be set to true if data_pointer( const BlockType& ) exist.
	   uses_plain_array_storage = 0,
	   //! Ordering inside the block ; only useful is_plain_array is true
	   is_row_major = BlockType::IsRowMajor
	}  ;
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
