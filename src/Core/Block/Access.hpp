/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_BLOCK_ACCESS_HPP
#define BOGUS_BLOCK_ACCESS_HPP

#include "Constants.hpp"

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
	//SFINAE, as we can't use decltype()

	template < typename BlockT >
	inline static typename BlockT::ConstTransposeReturnType get( const BlockT& src, bool = false )
	{ return src.transpose() ; }

	template < typename BlockT >
	inline static typename BlockTransposeTraits< BlockT >::ReturnType get( const BlockT& src, bool = false )
	{ return transpose_block( src ) ; }

	// Retry with BlockT::Base, but only if the previous one was not successful ( to avoid ambiguities )
	template < typename BlockT >
	inline static typename DisableIf< HasReturnType< BlockTransposeTraits< BlockT > >::Value,
	typename BlockTransposeTraits< typename BlockT::Base >::ReturnType >::ReturnType get( const BlockT& src, bool = false )
	{ return transpose_block( src ) ; }

	template < typename BlockT >
	inline static typename EnableIf< BlockTraits< BlockT >::is_self_transpose, const BlockT& >::ReturnType
	get( const BlockT& src, bool = false )
	{ return src ; }

	typedef char NonTransposableReturnType ;

	template < typename T >
	static NonTransposableReturnType get( const T &, ... ) ;

} ;

//! Check if there is a non-trivial specialization of BlockGetter<true>::get.
/*! Only useful to display better compilation error messages.
	False positive if sizeof( BlockType ) = sizeof( char ) and BlockType is not char.
*/
template < typename BlockType >
struct IsTransposable
{
public:
	enum {
		Value = ( sizeof(BlockGetter< true >::NonTransposableReturnType ) == sizeof( BlockType ) )
			|| 	( sizeof(BlockGetter< true >::NonTransposableReturnType ) !=
				  sizeof(BlockGetter< true >::template get< BlockType >( * static_cast< BlockType* >( 0 ) ) ) )
	} ;

} ;

template < bool RuntimeCheck, bool DoTranspose >
struct BlockTransposeOption : public BlockGetter< DoTranspose >{
} ;

template < bool IgnoredDoTranspose >
struct BlockTransposeOption< true, IgnoredDoTranspose > {
	template < typename BlockT >
	inline static BlockT get( const BlockT& src, bool doTranspose = false )
	{ return doTranspose ? BlockT( transpose_block( src ) ) : src ; }
} ;

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
