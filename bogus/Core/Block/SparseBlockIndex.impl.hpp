#ifndef BOGUS_SPARSE_BLOCK_INDEX_IMPL_HPP
#define BOGUS_SPARSE_BLOCK_INDEX_IMPL_HPP

#include "SparseBlockIndex.hpp"

#include "CompoundSparseBlockIndex.hpp"

namespace bogus {

template< typename Derived >
Derived & SparseBlockIndexBase< Derived >::derived()
{
    return static_cast< Derived& >( *this ) ;
}

template< typename Derived >
const Derived & SparseBlockIndexBase< Derived >::derived() const
{
    return static_cast< const Derived& >( *this ) ;
}

template< typename Derived >
typename SparseBlockIndexBase< Derived >::Index SparseBlockIndexBase< Derived >::innerSize() const
{
    return innerOffsetsArray().size() - 1  ;
}

template< typename Derived >
typename SparseBlockIndexBase< Derived >::Index SparseBlockIndexBase< Derived >::outerSize() const
{
    return derived().outerSize() ;
}

template< typename Derived >
const typename SparseBlockIndexBase< Derived >::InnerOffsetsType & SparseBlockIndexBase< Derived >::innerOffsetsArray() const
{
    return derived().innerOffsetsArray() ;
}

template< typename Derived >
bool SparseBlockIndexBase< Derived >::hasInnerOffsets() const
{
    return !innerOffsetsArray().empty();
}

template < bool Compressed, typename Index, typename BlockPtr  >
template < typename SourceDerived >
SparseBlockIndex< Compressed, Index, BlockPtr > & SparseBlockIndex< Compressed, Index, BlockPtr >::operator=(
        const SparseBlockIndexBase< SourceDerived > &source )
{
    clear() ;
    resizeOuter( source.outerSize() ) ;

    for( typename SourceDerived::Index i = 0 ; i < source.outerSize() ; ++i )
    {
        for( typename SourceDerived::InnerIterator it( source.derived(), i ) ;
             it ; ++ it )
        {
            insertBack( i, it.inner(), it.ptr() ) ;
        }
    }

    finalize() ;
    valid = source.valid ;
    if( source.hasInnerOffsets() ) {
        innerOffsets.resize( source.innerOffsetsArray().size() ) ;
        std::copy( source.innerOffsetsArray().begin(), source.innerOffsetsArray().end(), innerOffsets.begin() ) ;
    }

    return *this ;
}


} // namespace bogus

#endif

