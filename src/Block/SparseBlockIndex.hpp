#ifndef SPARSEBLOCKINDEX_HPP
#define SPARSEBLOCKINDEX_HPP


#include <vector>
#include <cassert>

namespace bogus
{


template < typename BlockT, bool Compressed >
struct SparseBlockIndex
{
	typedef typename BlockMatrixTraits< SparseBlockIndex< BlockT, Compressed > >::Index Index ;
	typedef Index BlockPtr ;

	SparseBlockIndex(
		const std::vector< Index > &innerOffsets_,
		const std::vector< Index > &outerOffsets_ )
		: innerOffsets( innerOffsets_ ), outerOffsets( outerOffsets_ )
	{}

	typedef const std::vector< BlockT >& Blocks ;
	typedef std::vector < std::pair< Index, BlockPtr > > Inner ;
	typedef std::vector < Inner > Outer ;

	Inner inner ;
	Outer outer ;
	const std::vector< Index > &innerOffsets ;
	const std::vector< Index > &outerOffsets ;

	void resizeOuter( Index size )
	{
		outer.resize( size ) ;
	}

	void insertBack( Index outIdx, Index inIdx, BlockPtr ptr )
	{
		outer[ outIdx ].push_back( std::make_pair( inIdx, ptr ) ) ;
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
} ;

template < typename BlockT >
struct SparseBlockIndex< BlockT, true >
{
	typedef typename BlockMatrixTraits< SparseBlockIndex< BlockT, true > >::Index Index ;
	typedef Index BlockPtr ;

	SparseBlockIndex(
		const std::vector< Index > &innerOffsets_,
		const std::vector< Index > &outerOffsets_ )
		: innerOffsets( innerOffsets_ ), outerOffsets( outerOffsets_ )
	{}

	typedef std::vector< Index > Inner ;
	typedef std::vector< Index > Outer ;

	Inner inner ;
	Outer outer ;
	const std::vector< Index > &innerOffsets ;
	const std::vector< Index > &outerOffsets ;

	void resizeOuter( Index size )
	{
		outer.resize( size+1 ) ;
	}

	void insertBack( Index outIdx, Index inIdx, BlockPtr ptr )
	{
		assert( ptr == inner.size() ) ;
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

	struct InnerIterator
	{
		InnerIterator( const SparseBlockIndex& index, Index outer )
			: m_it( index.outer[ outer ] ), m_end( index.outer[ outer + 1] ),
			  m_inner( index.inner )
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
		BlockPtr ptr() const { return m_it ; }

	private:
		Index m_it ;
		Index m_end ;
		const Inner& m_inner ;
	} ;
} ;


}

#endif // SPARSEBLOCKINDEX_HPP
