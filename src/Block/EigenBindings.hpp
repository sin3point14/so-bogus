#ifndef BLOCK_EIGENBINDINGS_HPP
#define BLOCK_EIGENBINDINGS_HPP

#include <Eigen/Core>

#include "BlockMatrix.hpp"

namespace bogus
{

template< typename EigenDerived >
struct BlockProductTraits
{
	typedef Eigen::Matrix< typename Eigen::internal::traits<EigenDerived>::Scalar, Eigen::Dynamic, 1 > Vec ;
	typedef Eigen::Matrix< typename Eigen::internal::traits<EigenDerived>::Scalar, 1, Eigen::Dynamic > RowVec ;
} ;

template < typename Derived, typename EigenDerived >
typename BlockProductTraits< EigenDerived >::Vec operator* ( const BlockMatrixBase< Derived >& lhs,
						 const Eigen::MatrixBase< EigenDerived > &rhs )
{
	assert( rhs.cols() == 1 ) ;
	assert( rhs.rows() == lhs.cols() ) ;
	typename BlockProductTraits< EigenDerived >::Vec res ( BlockProductTraits< EigenDerived >::Vec::Zero( lhs.rows() ) ) ;
	lhs.multiply( rhs, res ) ;
	return res ;
}

template < typename Derived, typename EigenDerived >
typename BlockProductTraits< EigenDerived >::RowVec operator* ( const Eigen::MatrixBase< EigenDerived > &lhs,
						 const BlockMatrixBase< Derived >& rhs )
{
	assert( lhs.rows() == 1 ) ;
	assert( lhs.cols() == rhs.rows() ) ;
	typename BlockProductTraits< EigenDerived >::RowVec res ( BlockProductTraits< EigenDerived >::RowVec::Zero( rhs.cols() ) ) ;
	Eigen::Transpose< typename BlockProductTraits< EigenDerived >::RowVec > resTrans ( res ) ;
	rhs.multiply( lhs.transpose(), resTrans, true ) ;
	return res ;
}

}

#endif // EIGENBINDINGS_HPP
