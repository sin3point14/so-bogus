/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_BLOCK_EXPRESSIONS_HPP
#define BOGUS_BLOCK_EXPRESSIONS_HPP

#include "BlockMatrix.hpp"

namespace bogus
{

//! Base class for Transpose views of a BlockObjectBase
template <typename MatrixT>
struct Transpose : public BlockObjectBase< Transpose< MatrixT > >
{
	typedef BlockMatrixTraits< Transpose< MatrixT > > Traits ;
	typedef typename Traits::PlainObjectType PlainObjectType ;
	typedef typename Traits::Index Index ;

	const PlainObjectType &matrix ;

	Transpose( const PlainObjectType &m ) : matrix( m.derived() ) {}

	const PlainObjectType& transpose() const { return matrix ; }
	const PlainObjectType& eval() const { return matrix.eval() ; }

	Index rows() const { return matrix.cols() ; }
	Index cols() const { return matrix.rows() ; }
} ;


template <typename MatrixT>
struct BlockMatrixTraits< Transpose< MatrixT > >
{
	enum { is_transposed = 1 } ;
	enum { is_temporary = 1 } ;

	typedef BlockMatrixTraits< MatrixT > OrigTraits;
	typedef typename OrigTraits::Index Index ;
	typedef typename OrigTraits::BlockPtr BlockPtr ;
	typedef typename OrigTraits::PlainObjectType PlainObjectType ;
	typedef typename OrigTraits::Scalar Scalar ;

	typedef PlainObjectType ConstTransposeReturnType ;
} ;

template < typename ObjectT, bool IsTemporary >
struct BlockStorage
{
	typedef const ObjectT& ConstValue ;
} ;
template < typename ObjectT >
struct BlockStorage< ObjectT, true >
{
	typedef const ObjectT  ConstValue ;
} ;


template < typename ObjectT >
struct BlockOperand
{
	typedef ObjectT ObjectType ;
	typedef typename ObjectT::PlainObjectType PlainObjectType ;

	typedef BlockMatrixTraits< ObjectT > Traits ;
	enum { do_transpose = Traits::is_transposed } ;
	typedef typename Traits::Scalar Scalar ;

	typename BlockStorage< ObjectT, Traits::is_temporary >::ConstValue object ;
	Scalar scaling ;

	BlockOperand( const ObjectT & o, Scalar s = 1 )
		: object(o), scaling(s)
	{}
} ;

template < template < typename LhsT, typename RhsT > class BlockOp, typename LhsMatrixT, typename RhsMatrixT>
struct BinaryBlockOp : public BlockObjectBase< BlockOp< LhsMatrixT, RhsMatrixT > >
{
	typedef BlockOperand< LhsMatrixT > Lhs ;
	typedef BlockOperand< RhsMatrixT > Rhs ;

	typedef typename Lhs::PlainObjectType PlainLhsMatrixType ;
	typedef typename Rhs::PlainObjectType PlainRhsMatrixType ;

	const Lhs lhs ;
	const Rhs rhs ;
	enum { transposeLhs = Lhs::do_transpose };
	enum { transposeRhs = Rhs::do_transpose };

	BinaryBlockOp( const LhsMatrixT& l,const RhsMatrixT& r,
				   typename Lhs::Scalar lscaling = 1, typename Lhs::Scalar rscaling = 1 )
		: lhs( l, lscaling ), rhs ( r, rscaling )
	{}
} ;

template <typename LhsMatrixT, typename RhsMatrixT>
struct Product : public BinaryBlockOp< Product, LhsMatrixT, RhsMatrixT >
{
	typedef BinaryBlockOp< bogus::Product, LhsMatrixT, RhsMatrixT > Base ;

	Product( const LhsMatrixT& l, const RhsMatrixT &r,
			  typename Base::Lhs::Scalar lscaling = 1, typename Base::Lhs::Scalar rscaling = 1 )
		: Base( l, r, lscaling, rscaling )
	{}

	typename Base::ConstTransposeReturnType transpose()
	{
		return typename Base::ConstTransposeReturnType (
					Base::rhs.object.transpose(), Base::lhs.object.transpose(),
					Base::rhs.scaling, Base::lhs.scaling ) ;
	}
} ;

template <typename LhsMatrixT, typename RhsMatrixT>
struct BlockMatrixTraits< Product< LhsMatrixT, RhsMatrixT > >
{
	enum { is_transposed = 0 } ;
	enum { is_temporary = 1 } ;

	typedef BlockMatrixTraits< LhsMatrixT > OrigTraits;
	typedef typename OrigTraits::Index Index ;
	typedef typename OrigTraits::BlockPtr BlockPtr ;
	typedef typename OrigTraits::PlainObjectType PlainObjectType ;
	typedef typename OrigTraits::Scalar Scalar ;

	typedef Product< typename RhsMatrixT::ConstTransposeReturnType,
					typename LhsMatrixT::ConstTransposeReturnType >
	ConstTransposeReturnType ;
} ;

template <typename LhsMatrixT, typename RhsMatrixT>
struct Addition : public BinaryBlockOp< Addition, LhsMatrixT, RhsMatrixT >
{
	typedef BinaryBlockOp< bogus::Addition, LhsMatrixT, RhsMatrixT > Base ;
	Addition( const LhsMatrixT& l, const RhsMatrixT &r,
			  typename Base::Lhs::Scalar lscaling = 1, typename Base::Lhs::Scalar rscaling = 1 )
		: Base( l, r, lscaling, rscaling )
	{}

	typename Base::ConstTransposeReturnType transpose()
	{
		return typename Base::ConstTransposeReturnType (
					Base::lhs.object.transpose(), Base::rhs.object.transpose(),
					Base::lhs.scaling, Base::rhs.scaling ) ;
	}
} ;

template <typename LhsMatrixT, typename RhsMatrixT>
struct BlockMatrixTraits< Addition< LhsMatrixT, RhsMatrixT > >
{
	enum { is_transposed = 0 } ;
	enum { is_temporary = 1 } ;

	typedef BlockMatrixTraits< LhsMatrixT > OrigTraits;
	typedef typename OrigTraits::Index Index ;
	typedef typename OrigTraits::BlockPtr BlockPtr ;
	typedef typename OrigTraits::PlainObjectType PlainObjectType ;
	typedef typename OrigTraits::Scalar Scalar ;

	typedef Addition< typename LhsMatrixT::ConstTransposeReturnType,
					 typename RhsMatrixT::ConstTransposeReturnType >
	ConstTransposeReturnType ;
} ;

template <typename MatrixT>
struct Scaling : public BlockObjectBase< Scaling< MatrixT > >
{
	typedef BlockOperand< MatrixT > Operand ;
	typedef typename Operand::PlainObjectType PlainOperandMatrixType ;
	Operand operand ;

	enum { transposeOperand = Operand::do_transpose };

	typedef BlockObjectBase< Scaling< MatrixT > > Base ;

	Scaling( const MatrixT &object, const typename MatrixT::Scalar scaling )
		: operand( object, scaling )
	{}

	typename Base::ConstTransposeReturnType transpose()
	{
		return typename Base::ConstTransposeReturnType (
					operand.object.transpose(), operand.scaling ) ;
	}

	typename Base::Index rows() const { return operand.object.rows() ; }
	typename Base::Index cols() const { return operand.object.rows() ; }
} ;

template <typename MatrixT>
struct BlockMatrixTraits< Scaling< MatrixT > >
{
	enum { is_transposed = 0 } ;
	enum { is_temporary = 1 } ;

	typedef BlockMatrixTraits< MatrixT > OrigTraits;
	typedef typename OrigTraits::Index Index ;
	typedef typename OrigTraits::BlockPtr BlockPtr ;
	typedef typename OrigTraits::PlainObjectType PlainObjectType ;
	typedef typename OrigTraits::Scalar Scalar ;

	typedef Scaling< typename MatrixT::ConstTransposeReturnType >
	ConstTransposeReturnType ;
} ;

template < typename ObjectT >
struct BlockOperand< Scaling< ObjectT > > : public BlockOperand< ObjectT >
{
	typedef BlockOperand< ObjectT > Base ;

	BlockOperand( const Scaling< ObjectT > & o, typename Base::Scalar s = 1 )
		: Base(o.operand.object, s*o.operand.scaling)
	{}
	BlockOperand( const typename Base::ObjectType & o, typename Base::Scalar s = 1 )
		: Base(o, s)
	{}
} ;

// Transpose and matrix/vector product return types
// Specialization of these structures should define a ReturnType if the operation is allowed

template< typename BlockT >
struct BlockTransposeTraits {} ;

template< typename BlockT >
struct SelfTransposeTraits {} ;

template< typename Derived >
struct BlockVectorProductTraits {} ;

}


#endif // EXPRESSIONS_HPP
