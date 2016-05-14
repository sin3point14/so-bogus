#ifndef BOGUS_COMPOUND_BLOCK_MATRIX_HPP
#define BOGUS_COMPOUND_BLOCK_MATRIX_HPP

#include "IterableBlockObject.hpp"

#include "Access.hpp"

namespace bogus {

template< bool ColWise, typename MatrixT1, typename MatrixT2  >
class CompoundBlockMatrix ;


// side-by-side (ColWise = true)
template< bool ColWise, typename MatrixT1, typename MatrixT2  >
class CompoundBlockMatrix : public IterableBlockObject<CompoundBlockMatrix< ColWise, MatrixT1, MatrixT2 > > {
public:
	typedef IterableBlockObject< CompoundBlockMatrix< ColWise, MatrixT1, MatrixT2 > > Base ;

	typedef typename Base::Index  Index ;
	typedef typename Base::Scalar Scalar ;
	typedef typename Base::ConstTransposeReturnType ConstTransposeReturnType ;

	CompoundBlockMatrix(const IterableBlockObject< MatrixT1 >& first, const IterableBlockObject< MatrixT2 >& second ) ;

	// BlockObjectBase

	Index rows() const { return m_first.rows() ; }
	Index cols() const { return m_first.cols()+m_second.cols() ; }

	Index blockRows( Index row ) const { return m_first.blockRows( row ) ; }
	Index blockCols( Index col ) const {
		const Index off = col - secondBegin() ;
		return off < 0 ? m_first.blockCols( col ) : m_second.blockCols( off ) ;
	}

	Index rowsOfBlocks() const { return m_first.rowsOfBlocks() ; }
	Index colsOfBlocks() const { return m_first.colsOfBlocks() + m_second.colsOfBlocks() ; }

	const Index *rowOffsets( ) const { return m_first.rowOffsets() ; }
	const Index *colOffsets( ) const { return data_pointer(m_offsets) ; }

	ConstTransposeReturnType transpose() const { return Transpose< CompoundBlockMatrix >( *this ) ; }

	template < bool DoTranspose, typename RhsT, typename ResT >
	void multiply( const RhsT& rhs, ResT& res, Scalar alpha = 1, Scalar beta = 0 ) const ;

	// Iterable Block Object

	Index size() const { return m_first.size() + m_second.size() ; }

	template <typename Func>
	void eachBlockOfRow( const Index row, Func func ) const {
		m_first .derived().eachBlockOfRow( row, func ) ;
		m_second.derived().eachBlockOfRow( row, func ) ;
	}
	template <typename Func>
	void eachBlockOfCol( const Index col, Func func ) const {
		const Index off = col - secondBegin() ;
		if( off < 0 )
			m_first .derived().eachBlockOfCol( col, func ) ;
		else
			m_second.derived().eachBlockOfCol( off, func ) ;
	}

	template < bool DoTranspose, typename RhsT, typename ResT, typename PreOp >
	void rowMultiplyPrecompose( const Index row, const RhsT& rhs, ResT& res, const PreOp &op ) const ;

	template < bool DoTranspose, typename RhsT, typename ResT, typename PostOp >
	void colMultiplyPostcompose( const Index col, const RhsT& rhs, ResT& res, const PostOp &op ) const ;


	// Accessors

	const IterableBlockObject< MatrixT1 >& first() const { return m_first; }
	const IterableBlockObject< MatrixT2 >& second() const { return m_second; }
	const Index* compoundOffsets() const { return m_compoundOffsets ; }
	Index secondBegin() const { return m_first.colsOfBlocks() ; }

private:
	const IterableBlockObject< MatrixT1 >& m_first;
	const IterableBlockObject< MatrixT2 >& m_second;
	Index m_compoundOffsets[3] ;
	std::vector< Index > m_offsets ;

};

template< bool ColWise, typename MatrixT1, typename MatrixT2  >
struct BlockMatrixTraits< CompoundBlockMatrix< ColWise, MatrixT1, MatrixT2 > >
		: public BlockMatrixTraits< BlockObjectBase< CompoundBlockMatrix< ColWise, MatrixT1, MatrixT2 > > >
{
	typedef BlockMatrixTraits< MatrixT1 > OrigTraits;
	typedef BlockMatrixTraits< MatrixT2 > OthrTraits;

	typedef typename OrigTraits::Scalar   Scalar;
	typedef typename OrigTraits::Index    Index;

	enum {
		RowsPerBlock = SwapIf<
			 ColWise || ((int)OrigTraits::RowsPerBlock) == (int)OthrTraits::RowsPerBlock,
			internal::DYNAMIC, OrigTraits::RowsPerBlock >::First,
		ColsPerBlock = SwapIf<
			!ColWise || ((int)OrigTraits::ColsPerBlock) == (int)OthrTraits::ColsPerBlock,
			internal::DYNAMIC, OrigTraits::ColsPerBlock >::First
	};


	typedef typename MatrixT1::BlockType  BlockType ; //FIXME: BlockSolverBase / PlainObjectType
} ;

template< typename MatrixT1, typename MatrixT2  >
class CompoundBlockMatrix<false,MatrixT1, MatrixT2>
		: public IterableBlockObject<CompoundBlockMatrix< false, MatrixT1, MatrixT2 > > {
	//TODO
} ;

} //bogus

#endif
