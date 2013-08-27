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

#include "../Block/BlockMatrixBase.hpp"
#include "../Block/Expressions.hpp"

#ifndef BOGUS_BLOCK_WITHOUT_LINEAR_SOLVERS
#include "EigenLinearSolvers.hpp"
#ifndef BOGUS_BLOCK_WITHOUT_EIGEN_SPARSE
#include "EigenSparseLinearSolvers.hpp"
#endif
#endif

#include "../Utils/CppTools.hpp"

namespace bogus
{

// transpose_block, is_zero

template< typename EigenDerived >
inline bool is_zero ( const Eigen::MatrixBase< EigenDerived >& block,
			   typename EigenDerived::Scalar precision )
{
	return block.isZero( precision ) ;
}

template< typename EigenDerived >
inline const typename EigenDerived::Scalar* data_pointer ( const Eigen::MatrixBase< EigenDerived >& block )
{
	return block.derived().data() ;
}

#ifndef BOGUS_BLOCK_WITHOUT_EIGEN_SPARSE
template < typename BlockT >
struct BlockTransposeTraits< Eigen::SparseMatrixBase < BlockT > > {
	typedef const Eigen::Transpose< const BlockT > ReturnType ;
} ;

template< typename EigenDerived >
inline const Eigen::Transpose< const EigenDerived >
transpose_block( const Eigen::SparseMatrixBase< EigenDerived >& block )
{
	return block.transpose() ;
}

template< typename EigenDerived >
inline bool is_zero ( const Eigen::SparseMatrixBase< EigenDerived >& block,
			   typename EigenDerived::Scalar precision )
{
	return block.isZero( precision ) ;
}

#endif

// Block traits for Eigen::Matrix

template< typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols >
struct BlockTraits < Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >
{
	typedef _Scalar Scalar ;
	enum {
		RowsAtCompileTime = _Rows,
		ColsAtCompileTime = _Cols,
		uses_plain_array_storage = 1,
		is_row_major = _Options & Eigen::RowMajorBit,
		is_self_transpose = ( _Rows == _Cols ) && ( _Rows == 1 )
	} ;

} ;

// Block/block product return type

template<
	typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols,
	typename _Scalar2, int _Rows2, int _Cols2, int _Options2, int _MaxRows2, int _MaxCols2,
	bool TransposeLhs, bool TransposeRhs >
struct BlockBlockProductTraits <
	   Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>,
	   Eigen::Matrix<_Scalar2, _Rows2, _Cols2, _Options2, _MaxRows2, _MaxCols2>,
		TransposeLhs, TransposeRhs >
{
	typedef Eigen::Matrix< _Scalar,
		SwapIf< TransposeLhs, _Rows, _Cols >::First,
		SwapIf< TransposeRhs, _Rows2, _Cols2 >::Second,
		_Options,
		SwapIf< TransposeLhs, _MaxRows, _MaxCols >::First,
		SwapIf< TransposeRhs, _MaxRows2, _MaxCols2 >::Second >
	ReturnType ;
} ;

// Matrix vector product traits and operator*

template< typename EigenDerived >
struct BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >
{
	typedef Eigen::Matrix< typename Eigen::internal::traits<EigenDerived>::Scalar, EigenDerived::RowsAtCompileTime, EigenDerived::ColsAtCompileTime > ResVec ;
} ;

template< typename Derived >
inline typename BlockVectorProductTraits< Eigen::MatrixBase< Derived > >::ResVec
get_mutable_vector( const Eigen::MatrixBase< Derived > & )
{
	return typename BlockVectorProductTraits< Eigen::MatrixBase< Derived > >::ResVec() ;
}

} //namespace bogus


template < typename Derived, typename EigenDerived >
typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::ResVec operator* ( const bogus::BlockObjectBase< Derived >& lhs,
						 const Eigen::MatrixBase< EigenDerived > &rhs )
{
	assert( rhs.rows() == lhs.cols() ) ;

	typedef typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::ResVec ResVec ;
	ResVec res ( lhs.rows(), rhs.cols() ) ;

	lhs.eval()->template multiply< Derived::is_transposed >( rhs.derived(), res, 1, 0 ) ;
	return res ;
}

template < typename Derived, typename EigenDerived >
typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::ResVec operator* ( const bogus::Scaling< Derived >& lhs,
						 const Eigen::MatrixBase< EigenDerived > &rhs )
{
	assert( rhs.rows() == lhs.cols() ) ;

	typedef typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::ResVec ResVec ;
	ResVec res ( lhs.rows(), rhs.cols() ) ;

	lhs.operand.object.eval()->template multiply< Derived::is_transposed >( rhs.derived(), res, lhs.operand.scaling, 0 ) ;
	return res ;
}

template < typename Derived, typename EigenDerived >
typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::ResVec operator* ( const Eigen::MatrixBase< EigenDerived > &lhs,
						 const bogus::BlockObjectBase< Derived >& rhs )
{
	assert( lhs.cols() == rhs.rows() ) ;

	typedef typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::ResVec ResVec ;
	ResVec res ( lhs.rows(), rhs.cols() ) ;

	Eigen::Transpose< ResVec > resTrans ( res ) ;
	rhs.eval()->template multiply< !Derived::is_transposed >( lhs.transpose(), resTrans, 1, 0 ) ;
	return res ;
}

template < typename Derived, typename EigenDerived >
typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::ResVec operator* ( const Eigen::MatrixBase< EigenDerived > &lhs,
						 const bogus::Scaling< Derived >& rhs )
{
	assert( lhs.cols() == rhs.rows() ) ;

	typedef typename bogus::BlockVectorProductTraits< Eigen::MatrixBase< EigenDerived > >::ResVec ResVec ;
	ResVec res ( lhs.rows(), rhs.cols() ) ;

	Eigen::Transpose< ResVec > resTrans ( res ) ;
	rhs.operand.object.eval()->template multiply< !Derived::is_transposed >( lhs.transpose(), resTrans, rhs.operand.scaling, 0 ) ;
	return res ;
}



#endif // EIGENBINDINGS_HPP
