/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver          
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>                       
 *
 * This Source Code Form is subject to the terms of the Mozilla Public 
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BLOCK_EIGENBINDINGS_HPP
#define BLOCK_EIGENBINDINGS_HPP

#include <Eigen/Core>

#ifndef BOGUS_BLOCK_WITHOUT_EIGEN_SPARSE
#if EIGEN_VERSION_AT_LEAST(3,1,0)
#include <Eigen/SparseCore>
#else
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Sparse>
#endif
#endif

#include "BlockMatrix.hpp"
#include "Expressions.hpp"

#ifndef BOGUS_BLOCK_WITHOUT_LINEAR_SOLVERS
#include "../Utils/EigenLinearSolvers.hpp"
#ifndef BOGUS_BLOCK_WITHOUT_EIGEN_SPARSE
#include "../Utils/EigenSparseLinearSolvers.hpp"
#endif
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

#ifndef BOGUS_BLOCK_WITHOUT_EIGEN_SPARSE
template < typename BlockT >
struct BlockTransposeTraits< Eigen::SparseMatrixBase < BlockT > > {
	typedef const Eigen::Transpose< const BlockT > ReturnType ;
} ;

template< typename EigenDerived >
const Eigen::Transpose< const EigenDerived >
transpose_block( const Eigen::SparseMatrixBase< EigenDerived >& block )
{
	return block.transpose() ;
}
#endif

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

	lhs.template multiply< false >( rhs, res ) ;
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

	lhs.matrix.template multiply< true >( rhs, res ) ;
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
	rhs.template multiply< true >( lhs.transpose(), resTrans ) ;
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
	rhs.matrix.template multiply< false >( lhs.transpose(), resTrans ) ;
	return res ;
}

#endif // EIGENBINDINGS_HPP
