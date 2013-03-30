#ifndef BLOCK_EIGENBINDINGS_HPP
#define BLOCK_EIGENBINDINGS_HPP

#include <Eigen/Core>

#include "BlockMatrix.hpp"
#include "Expressions.hpp"

namespace bogus
{

template< typename EigenDerived >
struct BlockProductTraits
{
	typedef Eigen::Matrix< typename Eigen::internal::traits<EigenDerived>::Scalar, Eigen::Dynamic, 1 > Vec ;
	typedef Eigen::Matrix< typename Eigen::internal::traits<EigenDerived>::Scalar, 1, Eigen::Dynamic > RowVec ;
} ;

}

template < typename Derived, typename EigenDerived >
typename bogus::BlockProductTraits< EigenDerived >::Vec operator* ( const bogus::BlockMatrixBase< Derived >& lhs,
						 const Eigen::MatrixBase< EigenDerived > &rhs )
{
	assert( rhs.cols() == 1 ) ;
	assert( rhs.rows() == lhs.cols() ) ;
	typename bogus::BlockProductTraits< EigenDerived >::Vec res (
				bogus::BlockProductTraits< EigenDerived >::Vec::Zero( lhs.rows() ) ) ;
	lhs.multiply( rhs, res ) ;
	return res ;
}

template < typename Derived, typename EigenDerived >
typename bogus::BlockProductTraits< EigenDerived >::Vec operator* ( const bogus::Transpose< bogus::BlockMatrixBase< Derived > >& lhs,
						 const Eigen::MatrixBase< EigenDerived > &rhs )
{
	assert( rhs.cols() == 1 ) ;
	assert( rhs.rows() == lhs.matrix.rows() ) ;
	typename bogus::BlockProductTraits< EigenDerived >::Vec res (
				bogus::BlockProductTraits< EigenDerived >::Vec::Zero( lhs.matrix.cols() ) ) ;
	lhs.matrix.multiply( rhs, res, true ) ;
	return res ;
}

template < typename Derived, typename EigenDerived >
typename bogus::BlockProductTraits< EigenDerived >::RowVec operator* ( const Eigen::MatrixBase< EigenDerived > &lhs,
						 const bogus::BlockMatrixBase< Derived >& rhs )
{
	assert( lhs.rows() == 1 ) ;
	assert( lhs.cols() == rhs.rows() ) ;
	typename bogus::BlockProductTraits< EigenDerived >::RowVec res (
				bogus::BlockProductTraits< EigenDerived >::RowVec::Zero( rhs.cols() ) ) ;
	Eigen::Transpose< typename bogus::BlockProductTraits< EigenDerived >::RowVec > resTrans ( res ) ;
	rhs.multiply( lhs.transpose(), resTrans, true ) ;
	return res ;
}

template < typename Derived, typename EigenDerived >
typename bogus::BlockProductTraits< EigenDerived >::RowVec operator* ( const Eigen::MatrixBase< EigenDerived > &lhs,
						 const bogus::Transpose< bogus::BlockMatrixBase< Derived > >& rhs )
{
	assert( lhs.rows() == 1 ) ;
	assert( lhs.cols() == rhs.matrix.cols() ) ;
	typename bogus::BlockProductTraits< EigenDerived >::RowVec res (
				bogus::BlockProductTraits< EigenDerived >::RowVec::Zero( rhs.matrix.rows() ) ) ;
	Eigen::Transpose< typename bogus::BlockProductTraits< EigenDerived >::RowVec > resTrans ( res ) ;
	rhs.matrix.multiply( lhs.transpose(), resTrans ) ;
	return res ;
}

#endif // EIGENBINDINGS_HPP
