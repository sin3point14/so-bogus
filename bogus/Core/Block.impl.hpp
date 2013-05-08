/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_BLOCK_IMPL_HPP
#define BOGUS_BLOCK_IMPL_HPP

#include "Block.hpp"

#ifndef BOGUS_WITHOUT_EIGEN
#include "Eigen/BlockBindings.hpp"
#endif

#include "Block/BlockMatrix.impl.hpp"
#include "Block/SparseBlockMatrix.impl.hpp"
#include "Block/SparseTranspose.impl.hpp"
#include "Block/SparseMatrixVectorProduct.impl.hpp"
#include "Block/SparseMatrixMatrixProduct.impl.hpp"

#endif
