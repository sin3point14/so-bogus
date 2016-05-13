#ifndef BOGUS_BLOCK_ZERO_HPP
#define BOGUS_BLOCK_ZERO_HPP

#include "BlockObjectBase.hpp"

namespace bogus {

template< typename Scalar >
class Zero ;

template < typename Scalar_ >
struct BlockMatrixTraits< Zero< Scalar_ > >
		: public BlockMatrixTraits< BlockObjectBase< Zero< Scalar_ > > >
{
	typedef Scalar_ Scalar ;

	enum {
		is_symmetric  = 1,
	} ;

	typedef Zero<Scalar>           PlainObjectType ;
	typedef const PlainObjectType& ConstTransposeReturnType ;
	typedef PlainObjectType        TransposeObjectType ;
} ;

//! Representation of the null matrix.
/*! Useful as argument to functions expecting a BlockObjectBase */
template< typename Scalar >
class Zero : public BlockObjectBase<Zero<Scalar> > {
public:
	typedef BlockObjectBase<Zero<Scalar> > Base ;

	typedef typename Base::Index Index ;
	typedef typename Base::ConstTransposeReturnType ConstTransposeReturnType ;

	explicit Zero( Index rows = 0 )
		: m_rows(rows)
	{
		m_offsets[0] = 0 ;
		m_offsets[1] = rows ;
	}

	Index rows() const { return m_rows ; }
	Index cols() const { return m_rows ; }

	Index blockRows( Index ) const { return rows() ; }
	Index blockCols( Index ) const { return cols() ; }

	Index rowsOfBlocks() const { return 1 ; }
	Index colsOfBlocks() const { return 1 ; }

	const Index *rowOffsets( ) const { return &m_offsets ; }
	const Index *colOffsets( ) const { return &m_offsets ; }

	ConstTransposeReturnType transpose() const { *this; }

	template < bool DoTranspose, typename RhsT, typename ResT >
	void multiply( const RhsT& , ResT& res, Scalar = 1, Scalar beta = 0 ) const
	{
		if( ( Scalar ) 0 == beta )
			res.setZero() ;
		if( ( Scalar ) 1 != beta )
			res *= beta ;
	}

private:
	const Index m_rows ;
		  Index m_offsets[2] ;

};

} //bogus

#endif
