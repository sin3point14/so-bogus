/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver          
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>                       
 *
 * This Source Code Form is subject to the terms of the Mozilla Public 
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_EIGEN_MATRIX_TRAITS_HPP
#define BOGUS_EIGEN_MATRIX_TRAITS_HPP

#include "LinearSolver.hpp"

#include <Eigen/Core>

namespace bogus
{

template< unsigned Dimension, typename ScalarType >
struct MatrixTraits
{
	typedef ScalarType Scalar ;
	typedef Eigen::Matrix< Scalar, Dimension, 1 > Vector ;
	typedef Eigen::Matrix< Scalar, Dimension, Dimension > Matrix ;

	typedef LU  < Eigen::MatrixBase< Matrix > > LUType ;
	typedef LDLT< Eigen::MatrixBase< Matrix > > LDLTType ;

	enum{ dimension = Dimension } ;

	static ScalarType np( const Vector & v )
	{ return v[0] ; }
	static ScalarType& np( Vector & v )
	{ return v[0] ; }

	static typename Vector::template ConstFixedSegmentReturnType< Dimension - 1 >::Type
	tp( const Vector & v )  { return v.template segment< Dimension - 1 >( 1 ) ; }
	static typename Vector::template FixedSegmentReturnType< Dimension - 1 >::Type
	tp(       Vector & v )  { return v.template segment< Dimension - 1 >( 1 ) ; }
} ;

}

#endif
