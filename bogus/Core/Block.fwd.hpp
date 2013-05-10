/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_BLOCK_FWD_HPP
#define BOGUS_BLOCK_FWD_HPP

#include "Block/Traits.hpp"

namespace bogus {

//! Flags for compile-time tuning of the behavior of objects such as SparseBlockMatrix
/*! Any combination if those is theoretically possible, using the binary or '|' operator */
namespace flags
{
	enum {
		NONE = 0,
		//! Use a compressed index instead of a sparse one>
		//! This adds some restrictions on the order in which elements can be inserted,
		//! but can be more efficient and allow interoperability with other formats
		//! such as MKL's BSR
		COMPRESSED = 0x1,
		//! Store and index blocks in a column major way
		COL_MAJOR = 0x2,
		//! Store only half the matrix, or rather the triangular part which verifies ( \c inner <= \c outer ),
		//! \c outer being the row and \c inner the column for row-major matrices
		SYMMETRIC = 0x4
	} ;
}

template < typename Derived >
struct BlockObjectBase ;

template < typename Derived >
class BlockMatrixBase ;

template < typename Derived >
class SparseBlockMatrixBase ;

template < typename BlockT, int Flags = flags::NONE >
class SparseBlockMatrix  ;


}

#endif
