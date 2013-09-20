/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

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

// transpose_block, is_zero, resize, set_identity

template< typename EigenDerived >
inline bool is_zero ( const Eigen::MatrixBase< EigenDerived >& block,
				 typename EigenDerived::Scalar precision )
{
	return block.isZero( precision ) ;
}

template< typename EigenDerived >
inline void set_identity ( Eigen::MatrixBase< EigenDerived >& block )
{
	block.derived().setIdentity( ) ;
}

template< typename EigenDerived >
inline void resize ( Eigen::MatrixBase< EigenDerived >& block, int rows, int cols )
{
	block.derived().resize( rows, cols ) ;
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

template < typename Scalar, int Options, typename Index >
inline void set_identity ( Eigen::SparseMatrix< Scalar, Options, Index >& block )
{
	return block.setIdentity( ) ;
}

template < typename Scalar, int Options, typename Index >
inline void resize ( Eigen::SparseMatrix< Scalar, Options, Index >& block, Index rows, Index cols )
{
	block.resize( rows, cols ) ;
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

	typedef Eigen::Matrix< _Scalar, _Cols, _Rows, _Options, _MaxCols, _MaxRows >
	TransposeStorageType ;

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

template< typename Derived >
inline typename BlockVectorProductTraits< Eigen::MatrixBase< Derived > >::ResVec
get_mutable_vector( const Eigen::MatrixBase< Derived > & )
{
	return typename Eigen::internal::plain_matrix_type<Derived>::type() ;
}

// Matrix vector product return types and operator*

template< typename Derived, typename EigenDerived, bool MatrixOnTheLeft >
struct BlockMatrixEigenProduct
		: public Eigen::ReturnByValue< BlockMatrixEigenProduct< Derived, EigenDerived, MatrixOnTheLeft > >
{
	typedef typename Eigen::ReturnByValue< BlockMatrixEigenProduct >::Index Index ;
	typedef typename Eigen::ReturnByValue< BlockMatrixEigenProduct >::Scalar Scalar ;

	BlockMatrixEigenProduct( const Derived &matrix, const EigenDerived &vector, Scalar scaling = 1 )
		: m_matrix( matrix) , m_vector( vector ), m_scaling ( scaling )
	{}

	template<typename Dest> inline void evalTo(Dest& dst) const
	{
		m_matrix.eval()->template multiply< Derived::is_transposed >( m_vector.derived(), dst, m_scaling, 0 ) ;
	}

	Index rows() const { return m_matrix.rows() ; }
	Index cols() const { return m_vector.cols() ; }

private:
	const Derived& m_matrix ;
	const EigenDerived &m_vector ;
	const Scalar m_scaling ;
};
template< typename Derived, typename EigenDerived >
struct BlockMatrixEigenProduct< Derived, EigenDerived, false >
		: public Eigen::ReturnByValue< BlockMatrixEigenProduct< Derived, EigenDerived, false > >
{
	typedef typename Eigen::ReturnByValue< BlockMatrixEigenProduct >::Index Index ;
	typedef typename Eigen::ReturnByValue< BlockMatrixEigenProduct >::Scalar Scalar ;

	BlockMatrixEigenProduct( const Derived &matrix, const EigenDerived &vector, Scalar scaling = 1 )
		: m_matrix( matrix) , m_vector( vector ), m_scaling ( scaling )
	{}

	template<typename Dest> inline void evalTo(Dest& dst) const
	{
		Eigen::Transpose< Dest > transposed( dst.transpose() ) ;
		m_matrix.eval()->template multiply< !Derived::is_transposed >( m_vector.transpose(), transposed, m_scaling, 0 ) ;
	}

	Index rows() const { return m_vector.rows() ; }
	Index cols() const { return m_matrix.cols() ; }

private:
	const Derived& m_matrix ;
	const EigenDerived &m_vector ;
	const Scalar m_scaling ;
};
} //namespace bogus

namespace Eigen {
namespace internal {
template < typename Derived, typename EigenDerived, bool MatrixOnTheLeft >
struct traits< bogus::BlockMatrixEigenProduct< Derived, EigenDerived, MatrixOnTheLeft > >
{
	typedef Eigen::Matrix<
		typename EigenDerived::Scalar,
		Eigen::Dynamic, EigenDerived::ColsAtCompileTime,
		0,
		Eigen::Dynamic, EigenDerived::MaxColsAtCompileTime >
	ReturnType ;
} ;
template < typename Derived, typename EigenDerived >
struct traits< bogus::BlockMatrixEigenProduct< Derived, EigenDerived, false > >
{
	typedef Eigen::Matrix<
		typename EigenDerived::Scalar,
		EigenDerived::RowsAtCompileTime, Eigen::Dynamic,
		Eigen::RowMajor,
		EigenDerived::MaxRowsAtCompileTime, Eigen::Dynamic >
	ReturnType ;
} ;
}
} //namespace Eigen

namespace bogus{

template < typename Derived, typename EigenDerived >
BlockMatrixEigenProduct< Derived, EigenDerived, true > operator* (
		const BlockObjectBase< Derived >& lhs,
		const Eigen::MatrixBase< EigenDerived > &rhs )
{
	assert( rhs.rows() == lhs.cols() ) ;
	return BlockMatrixEigenProduct< Derived, EigenDerived, true > ( lhs.derived(), rhs.derived() ) ;
}

template < typename Derived, typename EigenDerived >
BlockMatrixEigenProduct< Derived, EigenDerived, true > operator* (
		const Scaling< Derived >& lhs,
		const Eigen::MatrixBase< EigenDerived > &rhs )
{
	assert( rhs.rows() == lhs.cols() ) ;
	return BlockMatrixEigenProduct< Derived, EigenDerived, true > ( lhs.operand.object, rhs.derived(), lhs.operand.scaling ) ;
}

template < typename Derived, typename EigenDerived >
BlockMatrixEigenProduct< Derived, EigenDerived, false > operator* (
		const Eigen::MatrixBase< EigenDerived > &lhs,
		const BlockObjectBase< Derived >& rhs )
{
	assert( lhs.cols() == rhs.rows() ) ;
	return BlockMatrixEigenProduct< Derived, EigenDerived, false > ( rhs.derived(), lhs.derived() ) ;
}

template < typename Derived, typename EigenDerived >
BlockMatrixEigenProduct< Derived, EigenDerived, false > operator* (
		const Eigen::MatrixBase< EigenDerived > &lhs,
		const Scaling< Derived >& rhs )
{
	assert( lhs.cols() == rhs.rows() ) ;
	return BlockMatrixEigenProduct< Derived, EigenDerived, false > ( rhs.operand.object, lhs.derived(), rhs.operand.scaling ) ;
}


} //namespace bogus

#endif // EIGENBINDINGS_HPP
