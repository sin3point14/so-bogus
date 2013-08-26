/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef BOGUS_EIGEN_BLOCK_CONTAINERS_HPP
#define BOGUS_EIGEN_BLOCK_CONTAINERS_HPP

#define BOGUS_USE_ALLIGNED_ALLOCATOR( MatrixType ) \
	template<> struct ResizableSequenceContainer < MatrixType > { \
		typedef std::vector<MatrixType,Eigen::aligned_allocator<MatrixType> > Type ;\
	}

#include <Eigen/StdVector>

namespace bogus {

namespace internal {
	typedef Eigen::Matrix<double, 2, 2, Eigen::RowMajor> RowMatrix2d ;
	typedef Eigen::Matrix<float , 2, 2, Eigen::RowMajor> RowMatrix2f ;
	typedef Eigen::Matrix<double, 4, 4, Eigen::RowMajor> RowMatrix4d ;
	typedef Eigen::Matrix<float , 4, 4, Eigen::RowMajor> RowMatrix4f ;
} 

BOGUS_USE_ALLIGNED_ALLOCATOR( Eigen::Vector2d );
BOGUS_USE_ALLIGNED_ALLOCATOR( Eigen::Vector2f );
BOGUS_USE_ALLIGNED_ALLOCATOR( Eigen::Vector4d );
BOGUS_USE_ALLIGNED_ALLOCATOR( Eigen::Vector4f );
BOGUS_USE_ALLIGNED_ALLOCATOR( Eigen::Matrix2d );
BOGUS_USE_ALLIGNED_ALLOCATOR( Eigen::Matrix2f );
BOGUS_USE_ALLIGNED_ALLOCATOR( Eigen::Matrix4d );
BOGUS_USE_ALLIGNED_ALLOCATOR( Eigen::Matrix4f );
BOGUS_USE_ALLIGNED_ALLOCATOR( internal::RowMatrix2d );
BOGUS_USE_ALLIGNED_ALLOCATOR( internal::RowMatrix2f );
BOGUS_USE_ALLIGNED_ALLOCATOR( internal::RowMatrix4d );
BOGUS_USE_ALLIGNED_ALLOCATOR( internal::RowMatrix4f );

} //namespace bogus

#endif
