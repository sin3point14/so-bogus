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
	typedef BlockMatrixTraits< BlockObjectBase< Zero< Scalar_ > > > BaseTraits ;
	typedef typename BaseTraits::Index      Index;
	typedef typename BaseTraits::BlockPtr   BlockPtr;

	typedef Scalar_ Scalar ;

	enum {
		is_symmetric  = 1,
		is_transposed = 0,
		is_temporary  = 0
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

	Zero( Index rows, Index cols )
		: m_rows(rows), m_cols(cols), m_offset(0)
	{}

	Index rows() const { return rows() ; }
	Index cols() const { return cols() ; }

	Index blockRows( Index row ) const { return rows() ; }
	Index blockCols( Index col ) const { return cols() ; }

	Index rowsOfBlocks() const { return 1 ; }
	Index colsOfBlocks() const { return 1 ; }

	const Index *rowOffsets( ) const { return &m_offset ; }
	const Index *colOffsets( ) const { return &m_offset ; }

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
	const Index m_cols ;
	const Index m_offset ;

};

} //bogus

#endif
