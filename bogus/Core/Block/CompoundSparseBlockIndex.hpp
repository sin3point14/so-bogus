/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_COMPOUND_SPARSE_BLOCK_INDEX_HPP
#define BOGUS_COMPOUND_SPARSE_BLOCK_INDEX_HPP

#include <vector>

#include "SparseBlockIndex.hpp"

namespace bogus
{

template < typename FirstIndexType, typename SecondIndexType >
struct CompoundSparseBlockIndex : public SparseBlockIndexBase< CompoundSparseBlockIndex< FirstIndexType, SecondIndexType > >
{
	typedef typename FirstIndexType::Index Index ;
	typedef typename FirstIndexType::BlockPtr BlockPtr ;

	typedef SparseBlockIndexBase< CompoundSparseBlockIndex< FirstIndexType, SecondIndexType > > Base ;
	typedef typename Base::InnerOffsetsType InnerOffsetsType ;
	using Base::valid ;

	CompoundSparseBlockIndex( const SparseBlockIndexBase<FirstIndexType >& index1,
							  const SparseBlockIndexBase<SecondIndexType>& index2 )
		: Base( index1.valid && index2.valid ),
		  first( index1.derived() ), second( index2.derived() ),
		  innerOffsets( index1.innerOffsetsArray() )
	{
		assert( first.outerSize() == first.outerSize() ) ;
		assert( second.innerSize() == second.innerSize() ) ;
	}

	Index outerSize() const { return first.outerSize() ; }
	Index nonZeros() const { return first.nonZeros() + second.nonZeros() ; }

	const InnerOffsetsType& innerOffsetsArray() const { return innerOffsets ; }

	struct InnerIterator
	{
		InnerIterator( const CompoundSparseBlockIndex& index, Index outer )
			: m_it1( index.first, outer ), m_it2( index.second, outer )
		{
		}

		operator bool() const
		{
			return m_it1 || m_it2 ;
		}

		InnerIterator& operator++()
		{
			if( m_it1 ) ++ m_it1 ;
			else        ++ m_it2 ;
			return *this ;
		}

		Index inner() const {
			return m_it1 ? m_it1.inner() : m_it2.inner() ;
		}
		BlockPtr ptr() const {
			return m_it1 ? m_it1.ptr()   : m_it2.ptr() ;
		}

	private:
		typename FirstIndexType::InnerIterator m_it1 ;
		typename SecondIndexType::InnerIterator m_it2 ;
	} ;

	const FirstIndexType& first ;
	const SecondIndexType& second ;
	const InnerOffsetsType &innerOffsets ;
} ;

template < typename FirstIndexType, typename SecondIndexType >
struct SparseBlockIndexTraits< CompoundSparseBlockIndex< FirstIndexType, SecondIndexType > >
		: public SparseBlockIndexTraits< FirstIndexType >
{
} ;


}

#endif

