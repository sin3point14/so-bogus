/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_MAPPED_SPARSEBLOCKMATRIX_HPP
#define BOGUS_MAPPED_SPARSEBLOCKMATRIX_HPP

#include "CompressedSparseBlockIndex.hpp"
#include "SparseBlockMatrixBase.hpp"

namespace bogus
{

//! Specialization of BlockMatrixTraits for SparseBlockMatrix
template < typename BlockT, int Flags >
struct BlockMatrixTraits< MappedSparseBlockMatrix< BlockT, Flags > >
        : public BlockMatrixTraits< BlockObjectBase< MappedSparseBlockMatrix< BlockT, Flags > > >
{
    typedef BlockMatrixTraits< BlockObjectBase< MappedSparseBlockMatrix< BlockT, Flags > > > BaseTraits ;
    typedef typename BaseTraits::Index      Index;
    typedef typename BaseTraits::BlockPtr   BlockPtr;

    typedef BlockT BlockType ;
    typedef typename BlockTraits< BlockT >::Scalar Scalar ;
    typedef BlockType* BlocksArrayType ;

    enum {
        is_transposed  = 0,
        is_temporary   = 0,

        is_compressed  = Flags & flags::COMPRESSED,
        is_symmetric   = Flags & flags::SYMMETRIC,
        is_col_major   = Flags & flags::COL_MAJOR,
        flags          = Flags
    } ;

    typedef SparseBlockIndex< is_compressed, Index, BlockPtr > MajorIndexType ;

} ;

//! Mapped Sparse Block Matrix
/*!
  \tparam BlockT the type of the blocks of the matrix. Can be scalar, Eigen dense of sparse matrices,
  or basically anything provided a few functions are specialized
  \tparam Flags a combination of the values defined in \ref bogus::flags
  */
template < typename BlockT, int Flags >
class MappedSparseBlockMatrix : public  SparseBlockMatrixBase< MappedSparseBlockMatrix< BlockT, Flags > >
{
public:
    typedef SparseBlockMatrixBase< MappedSparseBlockMatrix< BlockT, Flags > > Base ;

    MappedSparseBlockMatrix() : Base() {}

    template < typename RhsT >
    MappedSparseBlockMatrix( const BlockObjectBase< RhsT >& rhs ) : Base()
    {
        Base::operator= ( rhs.derived() ) ;
    }

    template < typename RhsT >
    MappedSparseBlockMatrix& operator=( const BlockObjectBase< RhsT >& rhs )
    {
        return ( Base::operator= ( rhs.derived() ) ).derived() ;
    }

} ;

// Specialization for block matrix of MappedSparseBlockMatrix
template < typename BlockT, int Flags >
struct BlockTraits< MappedSparseBlockMatrix< BlockT, Flags > >
{
    typedef MappedSparseBlockMatrix< BlockT, Flags > BlockType ;
    typedef typename BlockType::Scalar Scalar ;

    enum {
        RowsAtCompileTime = internal::DYNAMIC,
        ColsAtCompileTime = internal::DYNAMIC,
        uses_plain_array_storage = 0,
        is_row_major = !BlockMatrixTraits< BlockType >::is_col_major
    }  ;
} ;

}

#endif // SPARSEBLOCKMATRIX_HH
