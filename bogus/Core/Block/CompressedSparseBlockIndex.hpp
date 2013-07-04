/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef COMPRESSED_SPARSE_BLOCK_INDEX_HPP
#define COMPRESSED_SPARSE_BLOCK_INDEX_HPP

#include "SparseBlockIndex.hpp"

namespace bogus
{

//! Compressed index, compatible with usual BSR/BSC formats
template< typename _Index, typename _BlockPtr >
struct SparseBlockIndex< true, _Index, _BlockPtr > : public SparseBlockIndexBase< SparseBlockIndex< true, _Index, _BlockPtr > >
{
	typedef _Index Index ;
	typedef _BlockPtr BlockPtr ;

	typedef SparseBlockIndexBase< SparseBlockIndex< true, _Index, _BlockPtr > > Base ;
	typedef typename Base::InnerOffsetsType InnerOffsetsType ;
	using Base::valid ;

	typedef std::vector< Index > Inner ;
	typedef std::vector< Index > Outer ;

	InnerOffsetsType innerOffsets ;

	//! Vector of inner indices
	Inner inner ;
	//! Vector encoding the start and end of inner vectors
	Outer outer ;
	//! Constant offset to add to block pointers (i.e. pointer to first block )
	BlockPtr base ;

	SparseBlockIndex( )
		: Base(), base(0)
	{}

	void resizeOuter( Index size )
	{
		outer.assign( size+1, 0 ) ;
	}
	void reserve( Index nnz)
	{
		inner.reserve( nnz ) ;
	}

	Index outerSize( ) const { return outer.size() - 1 ; }
	Index nonZeros() const { return inner.size() ; }


	const InnerOffsetsType& innerOffsetsArray() const { return innerOffsets ; }

	//! \warning Only works for back insertion, and a call to \ref finalize()
	//! is always required once insertion is finished
	void insertBack( Index outIdx, Index inIdx, BlockPtr ptr )
	{
		valid &= ( ptr == (BlockPtr) ( base + inner.size() ) )
				&& ( 0 == outer[ outIdx+1 ] || inIdx > inner.back() ) ;
		++outer[ outIdx+1 ] ;
		inner.push_back( inIdx ) ;
	}

	//! Finalizes the outer indices vector
	/*! Before calling this function, \c outer[i] contains the number of blocks
		 in the \c i th inner vector

		 After calling this functions, \c outer[i] contains the index of the start
		 of the \c i th inner vector in \ref inner, and \c outer[i+1] its end
	*/
	void finalize()
	{
		for( unsigned i = 1 ; i < outer.size() ; ++i )
		{
			outer[i] += outer[i-1] ;
		}
	}

	void clear()
	{
		outer.assign( outer.size(), 0 ) ;
		inner.clear() ;

		valid = true ;
		base = 0 ;
	}

	SparseBlockIndex &operator=( const SparseBlockIndex &compressed )
	{
		if( &compressed != this )
		{
			outer = compressed.outer ;
			inner = compressed.inner ;
			if( !compressed.innerOffsets.empty() )
				innerOffsets = compressed.innerOffsets ;
			base  = compressed.base ;
			valid = compressed.valid ;
		}
		return *this ;
	}

	template < typename SourceDerived >
	SparseBlockIndex &operator=( const SparseBlockIndexBase< SourceDerived > &source )
	{
		resizeOuter( source.outerSize() ) ;
		reserve( source.nonZeros() ) ;

		inner.clear() ;
		if( source.hasInnerOffsets() ) {
			innerOffsets.resize( source.innerOffsetsArray().size() ) ;
			std::copy( source.innerOffsetsArray().begin(), source.innerOffsetsArray().end(), innerOffsets.begin() ) ;
		}
		valid = source.valid ;

		for( typename SourceDerived::Index i = 0 ; i < source.outerSize() ; ++i )
		{
			for( typename SourceDerived::InnerIterator it( source.derived(), i ) ;
				 it ; ++ it )
			{
				if( inner.empty() ) base = it.ptr() ;
				insertBack( i, it.inner(), it.ptr() ) ;
			}
		}

		finalize() ;

		return *this ;
	}

	SparseBlockIndex & move( SparseBlockIndex &compressed )
	{
		if( &compressed != this )
		{
			// Want to have fun with gcc 4.6.3 ? Just swap the following statements !
			// Note: Swapping works perfectly with (at least) gcc 4.6.3 in non-optimized mode, gcc 4.8, gcc 4.2, clang 4.2
			inner.swap( compressed.inner );
			outer.swap( compressed.outer );

			if( !compressed.innerOffsets.empty() )
				innerOffsets.swap( compressed.innerOffsets ) ;
			base  = compressed.base ;
			valid = compressed.valid ;
			compressed.valid = false ;

		}
		return *this ;
	}

	template < typename SourceDerived >
	SparseBlockIndex &move( const SparseBlockIndexBase< SourceDerived > &source )
	{
		return ( *this = source ) ;
	}

	BlockPtr last( const Index outerIdx ) const
	{
		return base + outer[ outerIdx + 1 ] - 1  ;
	}

	Index size( const Index outerIdx ) const
	{
		return  outer[ outerIdx + 1 ] - outer[ outerIdx ] ;
	}

	//! Forward iterator
	struct InnerIterator
	{
		// Warning: This class does not implement the full RandomAccessIterator concept ;
		// only the operations that are required by std::lower_bound are implemented
		typedef std::random_access_iterator_tag iterator_category;
		typedef Index                           value_type;
		typedef ptrdiff_t                       difference_type;
		typedef const Index*                    pointer;
		typedef const Index&                    reference;

		InnerIterator( ) : m_inner( NULL ) { }

		InnerIterator( const SparseBlockIndex& index, Index outer )
			: m_it( index.outer[ outer ] ), m_end( index.outer[ outer + 1] ),
			  m_base( index.base ), m_inner( &index.inner[0] )
		{
		}

		operator bool() const
		{
			return m_it != m_end ;
		}

		InnerIterator& operator++()
		{
			++ m_it ;
			return *this ;
		}

		InnerIterator& operator+= ( const std::size_t n )
		{
			m_it = std::min( m_it + (Index) n, m_end ) ;
			return *this ;
		}

		difference_type operator- ( const InnerIterator& other ) const
		{
			return ( (difference_type) m_it ) - ( difference_type ) other.m_it ;
		}

		Index operator* () const
		{
			return inner() ;
		}

		InnerIterator end() const
		{
			return InnerIterator( m_end, m_end, m_base, m_inner ) ;
		}

		Index inner() const { return m_inner[ m_it ] ; }
		BlockPtr ptr() const { return m_it + m_base ; }

		BlockPtr rawIndex() const { return m_it ; }
	private:

		InnerIterator( Index it, Index end,
					   BlockPtr base, const Index* inner )
			: m_it( it ), m_end( end ), m_base( base ), m_inner( inner )
		{}

		Index m_it ;
		Index m_end ;
		BlockPtr m_base ;
		const Index* m_inner ;
	} ;

	void setPtr( const InnerIterator& it, BlockPtr ptr )
	{
		base = ptr - it.rawIndex() ;
		valid = true ;
	}

} ;

} //namespace bogus


#endif
