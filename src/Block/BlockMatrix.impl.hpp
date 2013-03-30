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
