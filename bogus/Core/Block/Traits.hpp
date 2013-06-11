/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef BOGUS_BLOCK_TRAITS_HPP
#define BOGUS_BLOCK_TRAITS_HPP

#include <vector>

#if !( defined( _OPENMP ) || defined( BOGUS_DONT_PARALLELIZE ) )
#define BOGUS_DONT_PARALLELIZE
#endif

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
} ;

} // namespace bogus

#endif
