#ifndef BLOCK_EIGENBINDINGS_HPP
#define BLOCK_EIGENBINDINGS_HPP

#include <Eigen/Core>

#include "BlockMatrix.hpp"
#include "Expressions.hpp"

namespace bogus
{

template< typename EigenDerived >
struct BlockProductTraits< Eigen::MatrixBase< EigenDerived > >
{
	typedef Eigen::Matrix< typename Eigen::internal::traits<EigenDerived>::Scalar, EigenDerived::RowsAtCompileTime, 1 > ResVec ;
	typedef Eigen::Matrix< typename Eigen::internal::traits<EigenDerived>::Scalar, 1, EigenDerived::ColsAtCompileTime > RowResVec ;
} ;

template< typename Derived >
typename BlockProductTraits< Eigen::MatrixBase< Derived > >::ResVec getBlockProductResVec( const Eigen::MatrixBase< Derived > & )
{
	return typename BlockProductTraits< Eigen::MatrixBase< Derived > >::ResVec() ;
}

}

template < typename Derived, typename EigenDerived >
typename bogus::BlockProductTraits< Eigen::MatrixBase< EigenDerived > >::ResVec operator* ( const bogus::BlockMatrixBase< Derived >& lhs,
						 const Eigen::MatrixBase< EigenDerived > &rhs )
{
	assert( rhs.cols() == 1 ) ;
	assert( rhs.rows() == lhs.cols() ) ;

	typedef typename bogus::BlockProductTraits< Eigen::MatrixBase< EigenDerived > >::ResVec ResVec ;

	ResVec res ( lhs.rows() ) ;
	res.setZero() ;

	lhs.multiply( rhs, res ) ;
	return res ;
}

template < typename Derived, typename EigenDerived >
typename bogus::BlockProductTraits< Eigen::MatrixBase< EigenDerived > >::ResVec operator* ( const bogus::Transpose< bogus::BlockMatrixBase< Derived > >& lhs,
						 const Eigen::MatrixBase< EigenDerived > &rhs )
{
	assert( rhs.cols() == 1 ) ;
	assert( rhs.rows() == lhs.matrix.rows() ) ;

	typedef typename bogus::BlockProductTraits< Eigen::MatrixBase< EigenDerived > >::ResVec ResVec ;
	ResVec res ( lhs.matrix.cols() ) ;
	res.setZero() ;

	lhs.matrix.multiply( rhs, res, true ) ;
	return res ;
}

template < typename Derived, typename EigenDerived >
typename bogus::BlockProductTraits< Eigen::MatrixBase< EigenDerived > >::RowResVec operator* ( const Eigen::MatrixBase< EigenDerived > &lhs,
						 const bogus::BlockMatrixBase< Derived >& rhs )
{
	assert( lhs.rows() == 1 ) ;
	assert( lhs.cols() == rhs.rows() ) ;

	typedef typename bogus::BlockProductTraits< Eigen::MatrixBase< EigenDerived > >::RowResVec ResVec ;
	ResVec res ( rhs.cols() ) ;
	res.setZero() ;

	Eigen::Transpose< ResVec > resTrans ( res ) ;
	rhs.multiply( lhs.transpose(), resTrans, true ) ;
	return res ;
}

template < typename Derived, typename EigenDerived >
typename bogus::BlockProductTraits< Eigen::MatrixBase< EigenDerived > >::RowResVec operator* ( const Eigen::MatrixBase< EigenDerived > &lhs,
						 const bogus::Transpose< bogus::BlockMatrixBase< Derived > >& rhs )
{
	assert( lhs.rows() == 1 ) ;
	assert( lhs.cols() == rhs.matrix.cols() ) ;

	typedef typename bogus::BlockProductTraits< Eigen::MatrixBase< EigenDerived > >::RowResVec ResVec ;
	ResVec res ( rhs.matrix.rows() ) ;
	res.setZero() ;

	Eigen::Transpose< ResVec > resTrans ( res ) ;
	rhs.matrix.multiply( lhs.transpose(), resTrans ) ;
	return res ;
}

#endif // EIGENBINDINGS_HPP
