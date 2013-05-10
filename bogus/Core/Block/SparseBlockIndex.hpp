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

namespace bogus
{

template< typename _Index >
struct SparseBlockIndexBase
{
	//! Vector encoding the size of each block of the inner dimension.
	/*! which can be retrieved as \code innerOffsets[inner+1] - innerOffsets[inner] \endcode */
	std::vector< _Index > innerOffsets ;
	//! Whether this index is currently valid
	bool valid ;

	SparseBlockIndexBase() : valid( true )
	{}
} ;

//! Uncompressed sparse block index
template < bool Compressed, typename _Index, typename _BlockPtr = _Index  >
struct SparseBlockIndex : public SparseBlockIndexBase< _Index >
{
	typedef SparseBlockIndexBase< _Index > Base ;
	using Base::innerOffsets ;
	using Base::valid ;

	typedef _Index Index ;
	typedef _BlockPtr BlockPtr ;

	//! Vector of ( inner index ; block pointer ) tuples encoding an inner vector
	typedef std::vector < std::pair< Index, BlockPtr > > Inner ;
	//! Vector of inner vectors
	typedef std::vector < Inner > Outer ;

	Outer outer ;

	SparseBlockIndex() : Base( )
	{}

	void resizeOuter( Index size )
	{
		outer.resize( size ) ;
	}
	Index outerSize( ) const { return outer.size() ; }
	Index innerSize( ) const { return innerOffsets.size() - 1 ; }

	void insertBack( Index outIdx, Index inIdx, BlockPtr ptr )
	{
		outer[ outIdx ].push_back( std::make_pair( inIdx, ptr ) ) ;
	}

	void finalize()
	{
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

	template < bool OCompressed, typename OIndex, typename OBlockPtr  >
	SparseBlockIndex &operator=( const SparseBlockIndex< OCompressed, OIndex, OBlockPtr > &compressed ) ;

	template < bool OCompressed, typename OIndex, typename OBlockPtr  >
	SparseBlockIndex& setToTranspose( const SparseBlockIndex< OCompressed, OIndex, OBlockPtr > &source )
	{
		clear() ;
		resizeOuter( source.innerSize() ) ;
		valid = source.valid ;

		for( OIndex i = 0 ; i < source.outerSize() ; ++i )
		{
			for( typename SparseBlockIndex< OCompressed, OIndex, OBlockPtr >::InnerIterator it( source, i ) ;
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

//! Compressed index, compatible with usual BSR/BSC formats
template< typename _Index, typename _BlockPtr >
struct SparseBlockIndex< true, _Index, _BlockPtr > : public SparseBlockIndexBase< _Index >
{
	typedef SparseBlockIndexBase< _Index > Base ;
	typedef _Index Index ;
	typedef _BlockPtr BlockPtr ;

	using Base::innerOffsets ;
	using Base::valid ;

	typedef std::vector< Index > Inner ;
	typedef std::vector< Index > Outer ;

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
	Index outerSize( ) const { return outer.size() - 1 ; }
	Index innerSize( ) const { return innerOffsets.size() - 1 ; }

	//! \warning Only works for back insertion, and a call to \ref finalize()
	//! is always required once insertion is finished
	void insertBack( Index outIdx, Index inIdx, BlockPtr ptr )
	{
		valid &= ( ptr == (BlockPtr) ( base + inner.size() ) ) ;
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

	SparseBlockIndex &operator=( SparseBlockIndex &compressed )
	{
		if( &compressed != this )
		{
			// We would like to swap this vector as well, but there seems to be a bug in 4.6.3
			// that prevent memormy ownership to be properly transfered
			// Note: Swapping works perfectly with gcc 4.6.3 in non-optimized mode, gcc 4.8, gcc 4.2, clang 4.2
			outer = compressed.outer ;
//			outer.swap( compressed.outer );
			inner.swap( compressed.inner );
			if( !compressed.innerOffsets.empty() )
				innerOffsets.swap( compressed.innerOffsets ) ;
			base  = compressed.base ;
			valid = compressed.valid ;
			compressed.valid = false ;
		}
		return *this ;
	}

	template < bool OCompressed, typename OIndex, typename OBlockPtr  >
	SparseBlockIndex &operator=( const SparseBlockIndex< OCompressed, OIndex, OBlockPtr > &source )
	{
		resizeOuter( source.outerSize() ) ;
		inner.clear() ;
		if( !source.innerOffsets.empty() ) {
			innerOffsets.resize( source.innerOffsets.size() ) ;
			std::copy( source.innerOffsets.begin(), source.innerOffsets.end(), innerOffsets.begin() ) ;
		}
		valid = source.valid ;

		for( Index i = 0 ; i < source.outerSize() ; ++i )
		{
			for( typename SparseBlockIndex< OCompressed, OIndex, OBlockPtr >::InnerIterator it( source, i ) ;
				 it ; ++ it )
			{
				if( inner.empty() ) base = it.ptr() ;
				insertBack( i, it.inner(), it.ptr() ) ;
			}
		}

		finalize() ;

		return *this ;
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

template < bool Compressed, typename Index, typename BlockPtr  >
template < bool OCompressed, typename OIndex, typename OBlockPtr  >
SparseBlockIndex< Compressed, Index, BlockPtr > & SparseBlockIndex< Compressed, Index, BlockPtr >::operator=(
		const SparseBlockIndex< OCompressed, OIndex, OBlockPtr> &source )
{
	clear() ;
	resizeOuter( source.outerSize() ) ;

	for( OIndex i = 0 ; i < source.outerSize() ; ++i )
	{
		for( typename SparseBlockIndex< OCompressed, OIndex, OBlockPtr >::InnerIterator it( source, i ) ;
			 it ; ++ it )
		{
			insertBack( i, it.inner(), it.ptr() ) ;
		}
	}

	finalize() ;
	valid = source.valid ;
	if( !source.innerOffsets.empty() ) {
		innerOffsets.resize( source.innerOffsets.size() ) ;
		std::copy( source.innerOffsets.begin(), source.innerOffsets.end(), innerOffsets.begin() ) ;
	}

	return *this ;
}

}


#endif // SPARSEBLOCKINDEX_HPP
