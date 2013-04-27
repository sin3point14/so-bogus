/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_EIGEN_PROBLEM_TRAITS_HPP
#define BOGUS_EIGEN_PROBLEM_TRAITS_HPP

#include "EigenMatrixTraits.hpp"

namespace bogus
{

template < typename LocalMatrixType >
struct ProblemTraits : public MatrixTraits< LocalMatrixType >
{
	typedef MatrixTraits< LocalMatrixType > Base ;
	typedef typename Base::Scalar Scalar ;
	enum{ dimension = Base::dimension } ;

	typedef Eigen::Matrix< Scalar, Eigen::Dynamic, 1 > DynVector ;

	template< typename VectorType >
	static typename VectorType::template FixedSegmentReturnType< dimension >::Type segment( const unsigned i, VectorType& v )
	{ return v.template segment< dimension > ( i * dimension ) ; }
	template< typename VectorType >
	static typename VectorType::template ConstFixedSegmentReturnType< dimension >::Type segment( const unsigned i, const VectorType& v )
	{ return v.template segment< dimension > ( i * dimension ) ; }

} ;

template< unsigned Dimension, typename Scalar >
struct LocalProblemTraits : public ProblemTraits< Eigen::Matrix< Scalar, Dimension, Dimension > >
{

	typedef Eigen::Matrix< Scalar, Dimension, 1 > Vector ;
	typedef Eigen::Matrix< Scalar, Dimension, Dimension > Matrix ;

	static Scalar np( const Vector & v )
	{ return v[0] ; }
	static Scalar& np( Vector & v )
	{ return v[0] ; }

	static typename Vector::template ConstFixedSegmentReturnType< Dimension - 1 >::Type
	tp( const Vector & v )  { return v.template segment< Dimension - 1 >( 1 ) ; }
	static typename Vector::template FixedSegmentReturnType< Dimension - 1 >::Type
	tp(       Vector & v )  { return v.template segment< Dimension - 1 >( 1 ) ; }
} ;


} //namespace bogus

#endif
