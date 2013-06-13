/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_SPARSEBLOCKINDEX_HPP
#define BOGUS_SPARSEBLOCKINDEX_HPP


#include <vector>
#include <cassert>
#include <algorithm>

namespace bogus
{

template< typename Derived >
struct SparseBlockIndexTraits
{
} ;

template< typename Derived >
struct SparseBlockIndexBase
{
	typedef typename SparseBlockIndexTraits< Derived >::Index Index ;

	//! Type of the array encoding the size of each block of the inner dimension.
	/*! which can be retrieved as \code innerOffsets[inner+1] - innerOffsets[inner] \endcode */
	typedef std::vector< Index > InnerOffsetsType ;

	//! Whether this index is currently valid
	bool valid ;

	SparseBlockIndexBase( bool _valid = true ) : valid( _valid )
	{}

	Derived& derived() ;
	const Derived& derived() const ;

	Index innerSize( ) const ;
	Index outerSize( ) const ;

	const InnerOffsetsType& innerOffsetsArray() const ;

	bool hasInnerOffsets() const;
} ;

//! Uncompressed sparse block index
template < bool Compressed, typename _Index, typename _BlockPtr = _Index  >
struct SparseBlockIndex : public SparseBlockIndexBase< SparseBlockIndex< Compressed, _Index, _BlockPtr > >
{

	typedef _Index Index ;
	typedef _BlockPtr BlockPtr ;

	typedef SparseBlockIndexBase<  SparseBlockIndex< Compressed, _Index, _BlockPtr > > Base ;
	typedef typename Base::InnerOffsetsType InnerOffsetsType ;
	using Base::valid ;

	//! Vector of ( inner index ; block pointer ) tuples encoding an inner vector
	typedef std::vector < std::pair< Index, BlockPtr > > Inner ;
	//! Vector of inner vectors
	typedef std::vector < Inner > Outer ;

	InnerOffsetsType innerOffsets ;
	Outer outer ;

	SparseBlockIndex() : Base( )
	{}

	void resizeOuter( Index size )
	{
		outer.resize( size ) ;
	}
	Index outerSize( ) const { return outer.size() ; }
	const InnerOffsetsType& innerOffsetsArray() const { return innerOffsets ; }

	void insertBack( Index outIdx, Index inIdx, BlockPtr ptr )
	{
		outer[ outIdx ].push_back( std::make_pair( inIdx, ptr ) ) ;
	}

	void finalize()
	{
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
		for( int i = 0 ; i < (int) outer.size() ; ++i )
		{
			std::sort( outer[i].begin(), outer[i].end() ) ;
		}
	}

	void clear()
	{
		std::vector< Inner >( outer.size() ).swap( outer ) ;
		valid = true ;
	}

	SparseBlockIndex &operator=( const SparseBlockIndex &o )
	{
		if( &o != this )
		{
			outer = o.outer ;
			valid = o.valid ;
			if( !o.innerOffsets.empty() )
				innerOffsets = o.innerOffsets ;
		}
		return *this ;
	}

	SparseBlockIndex &operator=( SparseBlockIndex &uncompressed )
	{
		outer.swap( uncompressed.outer ) ;
		if( !uncompressed.innerOffsets.empty() )
			innerOffsets.swap( uncompressed.innerOffsets ) ;
		valid = uncompressed.valid ;
		uncompressed.valid = false ;
		return *this ;
	}

	template < typename SourceDerived >
	SparseBlockIndex &operator=( const SparseBlockIndexBase< SourceDerived > &source ) ;

	template < typename SourceDerived >
	SparseBlockIndex& setToTranspose( const SparseBlockIndexBase< SourceDerived > &source )
	{
		clear() ;
		resizeOuter( source.innerSize() ) ;
		valid = source.valid ;

		for(  typename SourceDerived::Index i = 0 ; i < source.outerSize() ; ++i )
		{
			for( typename SourceDerived::InnerIterator it( source.derived(), i ) ;
				 it ; ++ it )
			{
				insertBack( it.inner(), i, it.ptr() ) ;
			}
		}
		return *this ;
	}

	BlockPtr last( const Index outerIdx ) const
	{
		return outer[ outerIdx ].back().second ;
	}

	Index size( const Index outerIdx ) const
	{
		return outer[ outerIdx ].size() ;
	}

	//! Forward iterator
	struct InnerIterator
	{
		// Warning: This class does not implement the full RandomAccessIterator concept ;
		// only the operations that are required by std::lower_bound are implemented
		typedef std::random_access_iterator_tag iterator_category;
		typedef Index                           value_type;
		typedef std::ptrdiff_t                  difference_type;
		typedef const Index*                    pointer;
		typedef const Index&                    reference;

		InnerIterator() {}

		InnerIterator( const SparseBlockIndex& index, Index outer )
			: m_it( index.outer[ outer ].begin() ), m_end( index.outer[ outer ].end() )
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
			m_it += n ;
			return *this ;
		}

		difference_type operator- ( const InnerIterator& other ) const
		{
			return m_it - other.m_it ;
		}

		Index operator* () const
		{
			return inner() ;
		}

		InnerIterator end() const
		{
			return InnerIterator( m_end, m_end ) ;
		}

		Index inner() const { return m_it->first ; }
		BlockPtr ptr() const { return m_it->second ; }

		typename Inner::const_iterator asStdIterator() const { return m_it ; }
	private:

		InnerIterator( const typename Inner::const_iterator &it,
				const typename Inner::const_iterator &end )
			: m_it( it ), m_end( end )
		{}

		typename Inner::const_iterator m_it ;
		typename Inner::const_iterator m_end ;
	} ;

	void setPtr( const InnerIterator& it, BlockPtr ptr )
	{
		const_cast< BlockPtr& >( it.asStdIterator()->second ) = ptr ;
	}

	template < typename VecT >
	typename VecT::SegmentReturnType innerSegment( VecT& v, Index idx ) const
	{
		return v.segment( innerOffsets[ idx ], innerOffsets[ idx + 1 ] - innerOffsets[ idx ] ) ;
	}
	template < typename VecT >
	typename VecT::ConstSegmentReturnType innerSegment( const VecT& v, Index idx ) const
	{
		return v.segment( innerOffsets[ idx ], innerOffsets[ idx + 1 ] - innerOffsets[ idx ] ) ;
	}

} ;

template < bool Compressed, typename _Index, typename _BlockPtr >
struct SparseBlockIndexTraits<  SparseBlockIndex< Compressed, _Index, _BlockPtr > >
{
	typedef _Index Index;
} ;



}


#endif // SPARSEBLOCKINDEX_HPP
