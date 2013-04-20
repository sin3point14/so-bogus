/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver          
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>                       
 *
 * This Source Code Form is subject to the terms of the Mozilla Public 
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_BLOCKMATRIX_IMPL_HPP
#define BOGUS_BLOCKMATRIX_IMPL_HPP

#include "BlockMatrix.hpp"

namespace bogus
{

template < typename Derived >
const Derived& BlockObjectBase< Derived >::derived() const
{
	return static_cast< const Derived& >( *this ) ;
}

template < typename Derived >
Derived& BlockObjectBase< Derived >::derived() {
	return static_cast< Derived& >( *this ) ;
}

}

#endif
