/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*! \file
	Necessary bindings to use Eigen matrices as block types, and
	\c operator* specialization for matrix/vector products
*/


#ifndef BLOCK_EIGENBINDINGS_HPP
#define BLOCK_EIGENBINDINGS_HPP

#include <Eigen/Core>

#ifndef BOGUS_BLOCK_WITHOUT_EIGEN_SPARSE
#include "SparseHeader.hpp"
#endif

#include "../Block/BlockMatrix.hpp"
#include "../Block/Expressions.hpp"

#ifndef BOGUS_BLOCK_WITHOUT_LINEAR_SOLVERS
#include "EigenLinearSolvers.hpp"
#ifndef BOGUS_BLOCK_WITHOUT_EIGEN_SPARSE
#include "EigenSparseLinearSolvers.hpp"
#endif
#endif

namespace bogus
{

// transpose_block, is_zero

template< typename EigenDerived >
typename EigenDerived::ConstTransposeReturnType
transpose_block ( const Eigen::MatrixBase< EigenDerived >& block )
{
	return block.transpose() ;
}

template< typename EigenDerived >
bool is_zero ( const Eigen::MatrixBase< EigenDerived >& block,
			   typename EigenDerived::Scalar precision )
{
	return block.isZero( precision ) ;
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

template< typename EigenDerived >
bool is_zero ( const Eigen::SparseMatrixBase< EigenDerived >& block,
			   typename EigenDerived::Scalar precision )
{
	return block.isZero( precision ) ;
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

} //namespace bogus


template < typename Derived, typename EigenDerived >
typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::ResVec operator* ( const bogus::BlockObjectBase< Derived >& lhs,
						 const Eigen::MatrixBase< EigenDerived > &rhs )
{
	assert( rhs.cols() == 1 ) ;
	assert( rhs.rows() == lhs.cols() ) ;

	typedef typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::ResVec ResVec ;
	ResVec res ( lhs.rows() ) ;

	lhs.eval().template multiply< Derived::is_transposed >( rhs, res, 1, 0 ) ;
	return res ;
}

template < typename Derived, typename EigenDerived >
typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::ResVec operator* ( const bogus::Scaling< Derived >& lhs,
						 const Eigen::MatrixBase< EigenDerived > &rhs )
{
	assert( rhs.cols() == 1 ) ;
	assert( rhs.rows() == lhs.cols() ) ;

	typedef typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::ResVec ResVec ;
	ResVec res ( lhs.rows() ) ;

	lhs.operand.object.eval().template multiply< Derived::is_transposed >( rhs, res, lhs.operand.scaling, 0 ) ;
	return res ;
}

template < typename Derived, typename EigenDerived >
typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::RowResVec operator* ( const Eigen::MatrixBase< EigenDerived > &lhs,
						 const bogus::BlockObjectBase< Derived >& rhs )
{
	assert( lhs.rows() == 1 ) ;
	assert( lhs.cols() == rhs.rows() ) ;

	typedef typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::RowResVec ResVec ;
	ResVec res ( rhs.cols() ) ;

	Eigen::Transpose< ResVec > resTrans ( res ) ;
	rhs.eval().template multiply< !Derived::is_transposed >( lhs.transpose(), resTrans, 1, 0 ) ;
	return res ;
}

template < typename Derived, typename EigenDerived >
typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::RowResVec operator* ( const Eigen::MatrixBase< EigenDerived > &lhs,
						 const bogus::Scaling< Derived >& rhs )
{
	assert( lhs.rows() == 1 ) ;
	assert( lhs.cols() == rhs.rows() ) ;

	typedef typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::RowResVec ResVec ;
	ResVec res ( rhs.cols() ) ;

	Eigen::Transpose< ResVec > resTrans ( res ) ;
	rhs.operand.object.eval().template multiply< !Derived::is_transposed >( lhs.transpose(), resTrans, rhs.operand.scaling, 0 ) ;
	return res ;
}



#endif // EIGENBINDINGS_HPP
