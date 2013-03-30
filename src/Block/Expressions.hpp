#ifndef BOGUS_BLOCK_EXPRESSIONS_HPP
#define BOGUS_BLOCK_EXPRESSIONS_HPP

namespace bogus
{

template <typename MatrixT>
struct Transpose
{
	MatrixT &matrix ;
	Transpose( MatrixT& m ) : matrix( m ) {}
} ;

template < typename MatrixT >
struct TransposeTraits
{
	enum { do_transpose = false } ;
} ;
template < typename MatrixT >
struct TransposeTraits< Transpose< MatrixT > >
{
	enum { do_transpose = true } ;
} ;

template <typename LhsMatrixT, typename RhsMatrixT>
struct Product
{
	LhsMatrixT &lhs ;
	RhsMatrixT &rhs ;
	enum { transposeLhs = TransposeTraits< LhsMatrixT >::do_transpose };
	enum { transposeRhs = TransposeTraits< RhsMatrixT >::do_transpose };

	Product( LhsMatrixT& l, LhsMatrixT &r ) : lhs( l ), rhs ( r ) {}
} ;

}

#endif // EXPRESSIONS_HPP
