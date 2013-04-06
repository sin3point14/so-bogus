#ifndef BOGUS_BLOCKMATRIX_HPP
#define BOGUS_BLOCKMATRIX_HPP

#include <vector>

namespace bogus
{

typedef unsigned Index ;


template < typename Derived >
struct BlockObjectBase
{
	typedef Derived PlainObjectType ;

	const Derived& derived() const ;
	Derived& derived() ;
};

template< typename Derived >
struct BlockMatrixTraits
{ } ;

template< typename Derived >
struct BlockMatrixTraits< BlockObjectBase< Derived > > {
	typedef unsigned Index ;
	typedef Derived PlainObjectType ;
} ;

struct BlockMatrixFlags
{
	enum {
		NONE = 0,
		COMPRESSED = 0x1,
		COL_MAJOR = 0x2,
		SYMMETRIC = 0x4
	} ;
} ;


template < typename Derived >
class BlockMatrixBase : public BlockObjectBase< Derived >
{
public:
	typedef typename BlockMatrixTraits< Derived >::BlockType BlockType ;
	typedef typename BlockMatrixTraits< Derived >::Index Index ;

public:
	BlockMatrixBase() : m_rows(0), m_cols(0)
	{}

	virtual ~BlockMatrixBase()
	{}

	template < typename RhsT, typename ResT >
	void multiply( const RhsT& rhs, ResT& res, bool transposed = false ) const
	{
		this->derived().multiply( rhs, res, transposed ) ;
	}

	template < typename RhsT, typename ResT >
	void splitRowMultiply( const Index row, const RhsT& rhs, ResT& res ) const
	{
		this->derived().splitRowMultiply( row, rhs, res ) ;
	}


	const BlockType& diagonal( const Index row ) const
	{
		return this->derived().diagonal( row );
	}

	Index rows() const { return m_rows ; }
	Index cols() const { return m_cols ; }
	Index rowsOfBlocks() const { return this->derived().rowsOfBlocks() ; }
	Index colsOfBlocks() const { return this->derived().colsOfBlocks() ; }

	const std::vector< BlockType >& blocks() const { return  m_blocks ; }

protected:
	Index m_rows ;
	Index m_cols ;

	std::vector< BlockType > m_blocks ;
} ;


}

#endif // BLOCKMATRIX_HPP
