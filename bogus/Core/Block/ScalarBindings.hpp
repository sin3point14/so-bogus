/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver          
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>                       
 *
 * This Source Code Form is subject to the terms of the Mozilla Public 
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef BOGUS_BLOCK_SCALAR_BINDGINS
#define BOGUS_BLOCK_SCALAR_BINDGINS

#include "Expressions.hpp"

#define BOGUS_BLOCK_SCALAR_TYPES \
	BOGUS_PROCESS_SCALAR( double   ) \
	BOGUS_PROCESS_SCALAR( float    ) \
	BOGUS_PROCESS_SCALAR( int      ) \

namespace bogus {

#define BOGUS_PROCESS_SCALAR( Scalar ) \
	template< > struct SelfTransposeTraits< Scalar > { typedef Scalar ReturnType ; } ;
BOGUS_BLOCK_SCALAR_TYPES
#undef BOGUS_PROCESS_SCALAR

#define BOGUS_PROCESS_SCALAR( Scalar ) \
	inline bool is_zero( Scalar s, Scalar precision ) { return std::abs( s ) <= precision ; }
BOGUS_BLOCK_SCALAR_TYPES
#undef BOGUS_PROCESS_SCALAR

#define BOGUS_PROCESS_SCALAR( Scalar_ ) \
    template< > struct BlockTraits< Scalar_ > { typedef Scalar_ Scalar ; } ;
BOGUS_BLOCK_SCALAR_TYPES
#undef BOGUS_PROCESS_SCALAR

}


#endif
