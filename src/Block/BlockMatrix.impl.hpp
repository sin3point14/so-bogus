#ifndef BLOCKMATRIX_IMPL_HPP
#define BLOCKMATRIX_IMPL_HPP

#include "BlockMatrix.hpp"

namespace bogus
{

template < typename Derived >
const Derived& BlockMatrixBase< Derived >::derived() const
{
	return static_cast< const Derived& >( *this ) ;
}

template < typename Derived >
Derived& BlockMatrixBase< Derived >::derived() {
	return static_cast< Derived& >( *this ) ;
}

}

#endif
