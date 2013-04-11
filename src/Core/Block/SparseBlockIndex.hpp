#ifndef BOGUS_SPARSEBLOCKINDEX_HPP
#define BOGUS_SPARSEBLOCKINDEX_HPP


#include <vector>
#include <cassert>

namespace bogus
{

template < bool Compressed = false >
struct SparseBlockIndex ;

struct SparseBlockIndexBase
{
	virtual bool isCompressed() const = 0 ;
	virtual const SparseBlockIndex<true>& asCompressed() const = 0 ;
	virtual const SparseBlockIndex<>& asUncompressed() const = 0 ;
} ;

template < bool Compressed  >
struct SparseBlockIndex : public SparseBlockIndexBase
{
	typedef unsigned Index ;
	typedef Index BlockPtr ;

	SparseBlockIndex() : valid( true )
	{}

	typedef std::vector < std::pair< Index, BlockPtr > > Inner ;
	typedef std::vector < Inner > Outer ;

	Outer outer ;
	std::vector< Index > innerOffsets ;
	bool valid ;

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

	SparseBlockIndex< Compressed > &operator=( const SparseBlockIndex< Compressed > &o )
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

	SparseBlockIndex< Compressed > &operator=( SparseBlockIndex< > &uncompressed )
	{
		outer.swap( uncompressed.outer ) ;
		if( !uncompressed.innerOffsets.empty() )
			innerOffsets.swap( uncompressed.innerOffsets ) ;
		valid = uncompressed.valid ;
		uncompressed.valid = false ;
		return *this ;
	}

	SparseBlockIndex< Compressed > &operator=( const SparseBlockIndex< true > &compressed ) ;

	template < bool OtherCompressed >
	SparseBlockIndex< Compressed >& setToTranspose( const SparseBlockIndex< OtherCompressed > &source )
	{
		clear() ;
		resizeOuter( source.innerSize() ) ;
		valid = source.valid ;

		for( unsigned i = 0 ; i < source.outerSize() ; ++i )
		{
			// For a symmetric matrix, do not store diagonal block in col-major index
			for( typename SparseBlockIndex< OtherCompressed >::InnerIterator it( source, i ) ;
				 it ; ++ it )
			{
				insertBack( it.inner(), i, it.ptr() ) ;
			}
		}
		return *this ;
	}

	bool isCompressed() const { return false ; }

	const SparseBlockIndex< > &asUncompressed () const
	{
		return *this ;
	}

	const SparseBlockIndex< true > &asCompressed () const
	{
		assert( 0 && "as Uncompressed should never be called on this object, segfaulting" ) ;
		return *static_cast< const SparseBlockIndex< true > * > ( 0 ) ;
	}

	BlockPtr last( const Index outerIdx ) const
	{
		return outer[ outerIdx ].back().second ;
	}

	Index size( const Index outerIdx ) const
	{
		return outer[ outerIdx ].size() ;
	}

	struct InnerIterator
	{
		InnerIterator( const SparseBlockIndex& index, Index outer )
			: m_it( index.outer[ outer ].begin() ), m_end( index.outer[ outer ].end() )
		{
		}

		operator bool()
		{
			return m_it != m_end ;
		}

		InnerIterator& operator++()
		{
			++ m_it ;
			return *this ;
		}

		Index inner() const { return m_it->first ; }
		BlockPtr ptr() const { return m_it->second ; }

		typename Inner::const_iterator asStdIterator() const { return m_it ; }
	private:
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

template<>
struct SparseBlockIndex< true > : public SparseBlockIndexBase
{
	typedef unsigned Index ;
	typedef Index BlockPtr ;

	SparseBlockIndex( )
		: base(0), valid( true )
	{}

	typedef std::vector< Index > Inner ;
	typedef std::vector< Index > Outer ;

	Inner inner ;
	Outer outer ;
	BlockPtr base ;
	std::vector< Index > innerOffsets ;
	bool valid ;

	void resizeOuter( Index size )
	{
		outer.assign( size+1, 0 ) ;
	}
	Index outerSize( ) const { return outer.size() - 1 ; }
	Index innerSize( ) const { return innerOffsets.size() - 1 ; }

	void insertBack( Index outIdx, Index inIdx, BlockPtr ptr )
	{
		valid &= ( ptr == base + inner.size() ) ;
		++outer[ outIdx+1 ] ;
		inner.push_back( inIdx ) ;
	}

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

	SparseBlockIndex< true > &operator=( const SparseBlockIndex< true > &compressed )
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

	SparseBlockIndex< true > &operator=( SparseBlockIndex< true > &compressed )
	{
		if( &compressed != this )
		{
			outer.swap( compressed.outer );
			inner.swap( compressed.inner );
			if( !compressed.innerOffsets.empty() )
				innerOffsets.swap( compressed.innerOffsets ) ;
			base  = compressed.base ;
			valid = compressed.valid ;
		}
		return *this ;
	}

	SparseBlockIndex< true > &operator=( const SparseBlockIndex< > &uncompressed )
	{
		resizeOuter( uncompressed.outerSize() ) ;
		inner.clear() ;
		if( !uncompressed.innerOffsets.empty() )
			innerOffsets = uncompressed.innerOffsets ;
		valid = uncompressed.valid ;

		for( unsigned i = 0 ; i < uncompressed.outerSize() ; ++i )
		{
			for( SparseBlockIndex< >::InnerIterator it( uncompressed, i ) ;
				 it ; ++ it )
			{
				if( inner.empty() ) base = it.ptr() ;
				insertBack( i, it.inner(), it.ptr() ) ;
			}
		}

		finalize() ;

		return *this ;
	}

	bool isCompressed() const { return true ; }

	const SparseBlockIndex< true > &asCompressed () const
	{
		return *this ;
	}

	const SparseBlockIndex< > &asUncompressed () const
	{
		assert( 0 && "as Uncompressed should never be called on this object, segfaulting" ) ;
		return *static_cast< const SparseBlockIndex< > * > ( 0 ) ;
	}

	BlockPtr last( const Index outerIdx ) const
	{
		return base + outer[ outerIdx + 1 ] - 1  ;
	}

	Index size( const Index outerIdx ) const
	{
		return  outer[ outerIdx + 1 ] - outer[ outerIdx ] ;
	}

	struct InnerIterator
	{
		InnerIterator( const SparseBlockIndex& index, Index outer )
			: m_it( index.outer[ outer ] ), m_end( index.outer[ outer + 1] ),
			  m_base( index.base ), m_inner( index.inner )
		{
		}

		operator bool()
		{
			return m_it != m_end ;
		}

		InnerIterator& operator++()
		{
			++ m_it ;
			return *this ;
		}

		Index inner() const { return m_inner[ m_it ] ; }
		BlockPtr ptr() const { return m_it + m_base ; }

		BlockPtr rawIndex() const { return m_it ; }
	private:
		Index m_it ;
		Index m_end ;
		BlockPtr m_base ;
		const Inner& m_inner ;
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

template < bool Compressed  >
SparseBlockIndex< Compressed > & SparseBlockIndex< Compressed >::operator=(
		const SparseBlockIndex< true > &compressed )
{
	clear() ;
	resizeOuter( compressed.outerSize() ) ;

	for( unsigned i = 0 ; i < compressed.outerSize() ; ++i )
	{
		for( typename SparseBlockIndex< true >::InnerIterator it( compressed, i ) ;
			 it ; ++ it )
		{
			insertBack( i, it.inner(), it.ptr() ) ;
		}
	}

	finalize() ;
	valid = compressed.valid ;
	if( !compressed.innerOffsets.empty() )
		innerOffsets = compressed.innerOffsets ;

	return *this ;
}

}

#endif // SPARSEBLOCKINDEX_HPP
