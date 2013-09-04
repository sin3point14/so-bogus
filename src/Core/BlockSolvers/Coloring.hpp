/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_COLORING_HPP
#define BOGUS_COLORING_HPP

#include "../Block.fwd.hpp"

namespace bogus {

//! Coloring
struct Coloring {

    std::vector< std::size_t > 	  permutation ;
    std::vector< std::ptrdiff_t > colors ;

    Coloring() : upToDate( false )
    {}

    void invalidate() { upToDate = false ; }

    template < typename Derived >
    void update( const bool enable, const BlockMatrixBase< Derived >& matrix ) ;


private:

    bool upToDate ;

    void reset( std::size_t n )
    {
        permutation.resize( n ) ;
        colors.clear() ;
        colors.push_back( 0 ) ;
        colors.push_back( n ) ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
        for( std::ptrdiff_t i = 0 ; i < (std::ptrdiff_t) n ; ++ i ) { permutation[i] = i ; }
    }

    template < typename Derived >
    void compute( const SparseBlockMatrixBase< Derived >& matrix ) ;

    template < typename Derived >
    void compute( const BlockMatrixBase< Derived >& matrix ) ;
} ;

}

#endif

