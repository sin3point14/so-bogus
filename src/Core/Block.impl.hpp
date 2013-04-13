#ifndef BOGUS_BLOCK_IMPL_HPP
#define BOGUS_BLOCK_IMPL_HPP

#include "Block.hpp"

#ifndef BOGUS_BLOCK_WITHOUT_EIGEN_BINDINGS
#include "Block/EigenBindings.hpp"
#endif
#ifndef BOGUS_BLOCK_WITHOUT_STREAMS
#include "Block/Streams.hpp"
#endif

#include "Block/BlockMatrix.impl.hpp"
#include "Block/SparseBlockMatrix.impl.hpp"
#include "Block/SparseTranspose.impl.hpp"
#include "Block/SparseProduct.impl.hpp"

#endif
