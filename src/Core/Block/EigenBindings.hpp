#ifndef BLOCK_EIGENBINDINGS_HPP
#define BLOCK_EIGENBINDINGS_HPP

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "BlockMatrix.hpp"
#include "Expressions.hpp"

#ifndef BOGUS_BLOCK_WITHOUT_LINEAR_SOLVERS
#include "../Utils/EigenLinearSolvers.hpp"
#include "../Utils/EigenSparseLinearSolvers.hpp"
#endif

namespace bogus
{

// Transpose

template< typename EigenDerived >
typename EigenDerived::ConstTransposeReturnType
transpose_block ( const Eigen::MatrixBase< EigenDerived >& block )
{
	return block.transpose() ;
}
template< typename EigenDerived >
const Eigen::Transpose< const EigenDerived >
transpose_block( const Eigen::SparseMatrixBase< EigenDerived >& block )
{
	return block.transpose() ;
}

template< typename Derived >
const LinearSolverBase< Derived >& transpose_block( const LinearSolverBase< Derived >& b )
{
	assert( 0 && "transpose_block should not be called on LinearSolvers ! ") ;
	return b ;
}


// Matrix vector product

template< typename EigenDerived >
struct BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >
{
	typedef Eigen::Matrix< typename Eigen::internal::traits<EigenDerived>::Scalar, EigenDerived::RowsAtCompileTime, 1 > ResVec ;
	typedef Eigen::Matrix< typename Eigen::internal::traits<EigenDerived>::Scalar, 1, EigenDerived::ColsAtCompileTime > RowResVec ;
} ;

template< typename Derived >
typename BlockVectorProductTraits< Eigen::MatrixBase< Derived > >::ResVec getBlockProductResVec( const Eigen::MatrixBase< Derived > & )
{
	return typename BlockVectorProductTraits< Eigen::MatrixBase< Derived > >::ResVec() ;
}

}

template < typename Derived, typename EigenDerived >
typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::ResVec operator* ( const bogus::BlockMatrixBase< Derived >& lhs,
						 const Eigen::MatrixBase< EigenDerived > &rhs )
{
	assert( rhs.cols() == 1 ) ;
	assert( rhs.rows() == lhs.cols() ) ;

	typedef typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::ResVec ResVec ;

	ResVec res ( lhs.rows() ) ;
	res.setZero() ;

	lhs.multiply( rhs, res ) ;
	return res ;
}

template < typename Derived, typename EigenDerived >
typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::ResVec operator* ( const bogus::Transpose< bogus::BlockMatrixBase< Derived > >& lhs,
						 const Eigen::MatrixBase< EigenDerived > &rhs )
{
	assert( rhs.cols() == 1 ) ;
	assert( rhs.rows() == lhs.matrix.rows() ) ;

	typedef typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::ResVec ResVec ;
	ResVec res ( lhs.matrix.cols() ) ;
	res.setZero() ;

	lhs.matrix.multiply( rhs, res, true ) ;
	return res ;
}

template < typename Derived, typename EigenDerived >
typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::RowResVec operator* ( const Eigen::MatrixBase< EigenDerived > &lhs,
						 const bogus::BlockMatrixBase< Derived >& rhs )
{
	assert( lhs.rows() == 1 ) ;
	assert( lhs.cols() == rhs.rows() ) ;

	typedef typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::RowResVec ResVec ;
	ResVec res ( rhs.cols() ) ;
	res.setZero() ;

	Eigen::Transpose< ResVec > resTrans ( res ) ;
	rhs.multiply( lhs.transpose(), resTrans, true ) ;
	return res ;
}

template < typename Derived, typename EigenDerived >
typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::RowResVec operator* ( const Eigen::MatrixBase< EigenDerived > &lhs,
						 const bogus::Transpose< bogus::BlockMatrixBase< Derived > >& rhs )
{
	assert( lhs.rows() == 1 ) ;
	assert( lhs.cols() == rhs.matrix.cols() ) ;

	typedef typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::RowResVec ResVec ;
	ResVec res ( rhs.matrix.rows() ) ;
	res.setZero() ;

	Eigen::Transpose< ResVec > resTrans ( res ) ;
	rhs.matrix.multiply( lhs.transpose(), resTrans ) ;
	return res ;
}

#endif // EIGENBINDINGS_HPP
