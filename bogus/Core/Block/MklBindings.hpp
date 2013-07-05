#ifndef BOGUS_BLOCK_MKL_BINDINGS_HPP
#define BOGUS_BLOCK_MKL_BINDINGS_HPP

#include "SparseBlockMatrix.hpp"

#include <mkl.h>

// Creates a compile error with boost
#ifdef P4
#undef P4
#endif

#include <iostream>

namespace bogus
{

template < bool NativeOrder, bool Transpose >
struct SparseBlockMatrixVectorMultiplier< true, false, NativeOrder, Transpose >
{
    template < typename Derived, typename RhsT, typename ResT, typename ScalarT >
    static void multiply( const SparseBlockMatrixBase< Derived >& matrix, const RhsT& rhs, ResT& res, ScalarT alpha )
    {
        std::cout << "HEY " << rhs.transpose()  << std::endl ;
    }
} ;


} //namespace bogus

#endif
