/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver          
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>                       
 *
 * This Source Code Form is subject to the terms of the Mozilla Public 
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_EIGEN_MATRIX_TRAITS_HPP
#define BOGUS_EIGEN_MATRIX_TRAITS_HPP

#include "../Utils/LinearSolverBase.hpp"

#include <Eigen/Core>

namespace bogus
{

template < typename _MatrixType >
struct MatrixTraits
{
	typedef _MatrixType MatrixType ;
	enum { dimension = MatrixType::RowsAtCompileTime  } ;
	typedef typename MatrixType::Scalar Scalar ;

	typedef LU  < Eigen::MatrixBase< MatrixType > > LUType ;
	typedef LDLT< Eigen::MatrixBase< MatrixType > > LDLTType ;

} ;

}

#endif
