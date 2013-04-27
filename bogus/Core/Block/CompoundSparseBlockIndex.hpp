/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver          
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>                       
 *
 * This Source Code Form is subject to the terms of the Mozilla Public 
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_COMPOUND_SPARSE_BLOCK_INDEX_HPP
#define BOGUS_COMPOUND_SPARSE_BLOCK_INDEX_HPP

namespace bogus
{

template < typename FirstIndexType, typename SecondIndexType >
struct CompoundSparseBlockIndex
{
	typedef typename FirstIndexType::Index Index ;
	typedef typename FirstIndexType::BlockPtr BlockPtr ;

	CompoundSparseBlockIndex( const FirstIndexType& index1, const SecondIndexType &index2 )
		: first( index1 ), second( index2 ), valid( first.valid && second.valid )
	{
		assert( index1.outerSize() == index2.outerSize() ) ;
		assert( index1.innerSize() == index2.innerSize() ) ;
	}

	Index outerSize() const { return first.outerSize() ; }
	Index innerSize() const { return first.innerSize() ; }

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
	bool valid ;
} ;

}

#endif

