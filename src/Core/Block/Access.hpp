/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_BLOCK_ACCESS_HPP
#define BOGUS_BLOCK_ACCESS_HPP

#include "Constants.hpp"
#include "Traits.hpp"

#include "../Utils/CppTools.hpp"

namespace bogus {

//! Specialization of transpose_block() for self-adjoint types
template < typename SelfTransposeT >
inline typename EnableIf< BlockTraits< SelfTransposeT >::is_self_transpose, const SelfTransposeT& >::ReturnType
transpose_block( const SelfTransposeT &block  ) { return  block ; }

//! Specialization of transpose_block() for types that define a ConstTransposeReturnType
template< typename BlockType >
typename BlockType::ConstTransposeReturnType
inline transpose_block ( const BlockType& block )
{
	return block.transpose() ;
}

//! Defines the transpose type of a \p BlockType using self-introspection
/*! Process af follow:
        - If the BlockType is self transpose, the transpose type is the BlockType itself
        - If the BlockType defines a ConstTransposeReturnType, use it
        - If BlockTransposeTraits< BlockType > defines a ReturnType, use it
        - If the BlockType defines a Base, retry with this Base
  */
template< typename BlockType,
          bool IsSelfTranspose = BlockTraits< BlockType >::is_self_transpose,
          bool DefinesConstTranspose = HasConstTransposeReturnType< BlockType >::Value,
          bool DefinesTransposeTraits = HasReturnType< BlockTransposeTraits< BlockType > >::Value,
          bool DefinesBase = HasBase< BlockType >::Value >
struct BlockTranspose {
    enum { is_defined= 0 } ;
} ;

// Self-transpose
template< typename BlockType, bool DCT, bool DTT, bool DB >
struct BlockTranspose< BlockType, true, DCT, DTT, DB > {
    typedef const BlockType& ReturnType ;
    enum { is_defined = 1 } ;
} ;
// ConstTransposeReturnType
template< typename BlockType, bool DTT, bool DB >
struct BlockTranspose< BlockType, false, true, DTT, DB > {
    typedef typename BlockType::ConstTransposeReturnType ReturnType ;
    enum { is_defined = 1 } ;
} ;
// BlockTransposeTraits
template< typename BlockType, bool DB >
struct BlockTranspose< BlockType, false, false, true, DB > {
    typedef typename BlockTransposeTraits< BlockType >::ReturnType ReturnType ;
    enum { is_defined = 1 } ;
} ;
// Base
template< typename BlockType >
struct BlockTranspose< BlockType, false, false, false, true >
        : public BlockTranspose< typename BlockType::Base >
{} ;

template < typename BlockType >
struct IsTransposable
{
    enum {
        Value = BlockTranspose< BlockType >::is_defined
    } ;
} ;


//! Utility struct for expressing a compile-time conditional transpose of a block
// In all of the following get functions, the dummy "bool = false" argument is there so
// that the specialization of BlockTransposeOption that does not perform a runtime check
// can just inherit from BlockGetter, and does not have to know the type returned by the get() function
template < bool DoTranspose >
struct BlockGetter {
	template < typename BlockT >
	inline static const BlockT& get( const BlockT& src, bool = false )
	{ return src ; }
} ;
template < >
struct BlockGetter< true > {
	template < typename BlockT >
    inline static typename BlockTranspose< BlockT >::ReturnType get( const BlockT& src, bool = false )
    { return transpose_block( src ) ; }
} ;

template < bool RuntimeCheck, bool DoTranspose >
struct BlockTransposeOption : public BlockGetter< DoTranspose >
{} ;

template < bool IgnoredDoTranspose >
struct BlockTransposeOption< true, IgnoredDoTranspose > {
	template < typename BlockT >
	inline static BlockT get( const BlockT& src, bool doTranspose = false )
	{ return doTranspose ? BlockT( transpose_block( src ) ) : src ; }
} ;

//! Access to the dimensions of a block
template< typename BlockT, bool Transpose_ = false >
struct BlockDims
{
	typedef BlockTraits< BlockT > Traits ;
	typedef SwapIf< Transpose_, Traits::RowsAtCompileTime, Traits::ColsAtCompileTime > Dims ;

		enum { Rows = Dims::First,
			   Cols = Dims::Second } ;
} ;

//! Access to segment of a vector corresponding to a given block-row
template < int DimensionAtCompileTime, typename VectorType, typename Index >
struct Segmenter
{
	enum { dimension = DimensionAtCompileTime } ;

	typedef typename VectorType::template NRowsBlockXpr< dimension >::Type
	ReturnType ;
	typedef typename VectorType::template ConstNRowsBlockXpr< dimension >::Type
	ConstReturnType ;

	Segmenter( VectorType &vec, const Index* ) : m_vec( vec ) {}

	inline ReturnType operator[]( const Index inner )
	{
		return m_vec.template middleRows< dimension >( dimension*inner ) ;
	}

	inline ConstReturnType operator[]( const Index inner ) const
	{
		return m_vec.template middleRows< dimension >( dimension*inner ) ;
	}

private:
	VectorType &		 m_vec ;
} ;

template < typename VectorType, typename Index >
struct Segmenter< internal::DYNAMIC, VectorType, Index >
{
	typedef typename VectorType::RowsBlockXpr			ReturnType ;
	typedef typename VectorType::ConstRowsBlockXpr		ConstReturnType ;

	Segmenter( VectorType &vec, const Index* offsets ) : m_vec( vec ), m_offsets( offsets ) { }

	inline ReturnType operator[]( const Index inner )
	{
		return m_vec.middleRows( m_offsets[ inner ], m_offsets[ inner + 1 ] - m_offsets[ inner ] ) ;
	}

	inline ConstReturnType operator[]( const Index inner ) const
	{
		return m_vec.middleRows( m_offsets[ inner ], m_offsets[ inner + 1 ] - m_offsets[ inner ] ) ;
	}


private:
	VectorType &		 m_vec ;
	const Index* m_offsets ;
} ;

}

#endif
