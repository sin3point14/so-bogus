#ifndef SPARSEBLOCKINDEX_HPP
#define SPARSEBLOCKINDEX_HPP


#include <vector>
#include <cassert>

namespace bogus
{


template < typename BlockT, bool Compressed = false >
struct SparseBlockIndex
{
	typedef typename BlockMatrixTraits< SparseBlockIndex< BlockT, Compressed > >::Index Index ;
	typedef Index BlockPtr ;

	SparseBlockIndex( const std::vector< Index > &innerOffsets_ )
		: innerOffsets( innerOffsets_ )
	{}

	typedef const std::vector< BlockT >& Blocks ;
	typedef std::vector < std::pair< Index, BlockPtr > > Inner ;
	typedef std::vector < Inner > Outer ;

	Outer outer ;
	const std::vector< Index > &innerOffsets ;

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

	void assign( SparseBlockIndex< BlockT > &uncompressed, BlockPtr )
	{
		outer.swap( uncompressed.outer ) ;
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

	private:
		typename Inner::const_iterator m_it ;
		typename Inner::const_iterator m_end ;
	} ;

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

template < typename BlockT >
struct SparseBlockIndex< BlockT, true >
{
	typedef typename BlockMatrixTraits< SparseBlockIndex< BlockT, true > >::Index Index ;
	typedef Index BlockPtr ;

	SparseBlockIndex( const std::vector< Index > &innerOffsets_ )
		: base(0), innerOffsets( innerOffsets_ )
	{}

	typedef std::vector< Index > Inner ;
	typedef std::vector< Index > Outer ;

	Inner inner ;
	Outer outer ;
	BlockPtr base ;
	const std::vector< Index > &innerOffsets ;

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

	void assign( SparseBlockIndex< BlockT > &uncompressed, BlockPtr base_ptr )
	{
		base = base_ptr ;

		resizeOuter( uncompressed.outer.size() ) ;

		for( unsigned i = 0 ; i < uncompressed.outer.size() ; ++i )
		{
			for( typename SparseBlockIndex< BlockT >::InnerIterator it( uncompressed, i ) ;
				 it ; ++ it )
			{
				insertBack( i, it.inner(), base+inner.size() ) ;
			}
		}

		finalize() ;
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

	private:
		Index m_it ;
		Index m_end ;
		BlockPtr m_base ;
		const Inner& m_inner ;
	} ;

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
