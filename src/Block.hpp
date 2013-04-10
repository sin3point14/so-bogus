#ifndef BOGUS_BLOCK_HPP
#define BOGUS_BLOCK_HPP

#include "Block/BlockMatrix.hpp"
#include "Block/Expressions.hpp"
#include "Block/SparseBlockMatrix.hpp"

#include "Block/BlockMatrix.impl.hpp"
#include "Block/SparseBlockMatrix.impl.hpp"
#include "Block/SparseTranspose.impl.hpp"
#include "Block/SparseProduct.impl.hpp"

#ifndef BOGUS_BLOCK_WITHOUT_EIGEN_BINDINGS
#include "Block/EigenBindings.hpp"
#endif

#endif
