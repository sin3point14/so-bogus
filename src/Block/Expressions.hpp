#ifndef BOGUS_BLOCK_EXPRESSIONS_HPP
#define BOGUS_BLOCK_EXPRESSIONS_HPP

#include "BlockMatrix.hpp"

namespace bogus
{

template <typename MatrixT>
struct Transpose : public BlockObjectBase< Transpose< MatrixT > >
{
	typedef typename MatrixT::PlainObjectType PlainObjectType ;
	const PlainObjectType &matrix ;

	Transpose( const PlainObjectType &m ) : matrix( m.derived() ) {}

	const PlainObjectType& transpose() { return matrix ; }
} ;

template < typename Derived >
class SparseBlockMatrixBase ;

template <typename MatrixT>
struct Transpose< SparseBlockMatrixBase< MatrixT > > : public Transpose< BlockMatrixBase< MatrixT > >
{
	Transpose( const SparseBlockMatrixBase< MatrixT > &m )
		: Transpose< BlockMatrixBase< MatrixT > > ( m.derived() )
	{}
} ;

template < typename MatrixT >
struct TransposeTraits
{
	enum { do_transpose = false } ;
	typedef MatrixT MatrixType ;

	static const MatrixType &get( const MatrixT &m ) { return m ; }
} ;
template < typename MatrixT >
struct TransposeTraits< Transpose< MatrixT > >
{
	enum { do_transpose = true } ;
	typedef typename MatrixT::PlainObjectType MatrixType ;

	static const MatrixType &get( const Transpose< MatrixT > &m ) { return m.matrix.derived() ; }
} ;

template <typename LhsMatrixT, typename RhsMatrixT>
struct Product
{
	typedef TransposeTraits< LhsMatrixT > LhsTraits ;
	typedef TransposeTraits< RhsMatrixT > RhsTraits ;
	const typename LhsTraits::MatrixType &lhs ;
	const typename RhsTraits::MatrixType &rhs ;
	enum { transposeLhs = LhsTraits::do_transpose };
	enum { transposeRhs = RhsTraits::do_transpose };

	Product( const BlockObjectBase< LhsMatrixT >& l, const BlockObjectBase< RhsMatrixT > &r )
		: lhs( LhsTraits::get( l.derived() ) ), rhs ( RhsTraits::get( r.derived() ) )
	{}
} ;

template < bool Symmetric, bool DoTranspose >
struct BlockTranspose{
	template < typename BlockT >
	static const BlockT& get( const BlockT& src, bool )
	{ return src ; }
} ;
template <  >
struct BlockTranspose< false, true > {
	template < typename BlockT >
	static typename BlockT::ConstTransposeReturnType get( const BlockT& src, bool )
	{ return src.transpose() ; }
} ;
template < >
struct BlockTranspose< true, true > {
	template < typename BlockT >
	static BlockT get( const BlockT& src, bool afterDiag )
	{ return afterDiag ? src : src.transpose() ; }
} ;
template < >
struct BlockTranspose< true, false > {
	template < typename BlockT >
	static BlockT get( const BlockT& src, bool afterDiag )
	{ return afterDiag ? src.transpose() : src ; }
} ;

}

#endif // EXPRESSIONS_HPP
