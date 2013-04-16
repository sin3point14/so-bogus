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

// SFINAE transpose and matrix/vector product return types

template< typename BlockT >
struct SelfTransposeTraits {} ;

template< typename BlockT >
struct BlockTransposeTraits {} ;

template< > struct SelfTransposeTraits< double   > { typedef double   ReturnType ; } ;
template< > struct SelfTransposeTraits< float    > { typedef float    ReturnType ; } ;
template< > struct SelfTransposeTraits< int      > { typedef int      ReturnType ; } ;
template< > struct SelfTransposeTraits< unsigned > { typedef unsigned ReturnType ; } ;

template < typename SelfTransposeT >
const typename SelfTransposeTraits< SelfTransposeT >::ReturnType& transpose_block( const SelfTransposeT &block  ) { return  block ; }


template< typename Derived >
struct BlockVectorProductTraits
{} ;

}


#endif // EXPRESSIONS_HPP
