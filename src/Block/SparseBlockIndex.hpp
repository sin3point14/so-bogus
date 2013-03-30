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

	SparseBlockIndex()
	{}

	typedef std::vector < std::pair< Index, BlockPtr > > Inner ;
	typedef std::vector < Inner > Outer ;

	Outer outer ;
	std::vector< Index > innerOffsets ;

	void resizeOuter( Index size )
	{
		outer.resize( size ) ;
	}

	void insertBack( Index outIdx, Index inIdx, BlockPtr ptr )
	{
		outer[ outIdx ].push_back( std::make_pair( inIdx, ptr ) ) ;
	}

	void finalize()
	{
	}

	void assign( SparseBlockIndex< > &uncompressed )
	{
		outer.swap( uncompressed.outer ) ;
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
		: base(0)
	{}

	typedef std::vector< Index > Inner ;
	typedef std::vector< Index > Outer ;

	Inner inner ;
	Outer outer ;
	BlockPtr base ;
	std::vector< Index > innerOffsets ;

	void resizeOuter( Index size )
	{
		outer.resize( size+1 ) ;
	}

	void insertBack( Index outIdx, Index inIdx, BlockPtr ptr )
	{
		assert( ptr == base + inner.size() ) ;
		(void) ptr ;
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

	void assign( SparseBlockIndex< > &uncompressed )
	{
		resizeOuter( uncompressed.outer.size() ) ;

		for( unsigned i = 0 ; i < uncompressed.outer.size() ; ++i )
		{
			for( typename SparseBlockIndex< >::InnerIterator it( uncompressed, i ) ;
				 it ; ++ it )
			{
				insertBack( i, it.inner(), base+inner.size() ) ;
			}
		}

		finalize() ;
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
		std::cout << base << std::endl ;
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


}

#endif // SPARSEBLOCKINDEX_HPP
