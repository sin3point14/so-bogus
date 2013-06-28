/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver          
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>                       
 *
 * This Source Code Form is subject to the terms of the Mozilla Public 
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef BOGUS_SPARSE_SCALEADD_HPP
#define BOGUS_SPARSE_SCALEADD_HPP

#include "SparseBlockMatrix.hpp"
#include "BlockTranspose.hpp"
#include "SparseBlockIndexComputer.hpp"
namespace bogus {

template < typename Derived >
Derived& SparseBlockMatrixBase< Derived >::scale( const Scalar alpha )
{

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
    for( int i = 0 ; i < (int) blocks().size() ;  ++i )
    {
        block( i ) *= alpha ;
    }

    return derived() ;
}

template < typename Derived >
template < bool Transpose, typename OtherDerived >
Derived& SparseBlockMatrixBase<Derived>::add( const SparseBlockMatrixBase< OtherDerived > &rhs, Scalar alpha )
{
	typedef typename SparseBlockMatrixBase< OtherDerived >::Traits OtherTraits ;
	typedef std::pair< BlockPtr, typename OtherTraits::BlockPtr > PtrPair ;
	typedef std::pair< Index, PtrPair > NonZero ;

	std::vector< std::vector< NonZero > > nonZeros ;

	// I - Compute non-zeros
	{

		SparseBlockIndexComputer< OtherDerived, OtherTraits::is_symmetric, Traits::is_col_major, Transpose >
				indexComputer( rhs ) ;
		typedef typename SparseBlockIndexComputer< OtherDerived, OtherTraits::is_symmetric, Traits::is_col_major, Transpose >::ReturnType
				SourceIndexType ;
		const SourceIndexType &rhsIndex = indexComputer.get() ;

		const MajorIndexType &lhsIndex = majorIndex() ;

		assert( rhsIndex.outerSize() == lhsIndex.outerSize() ) ;
		assert( rhsIndex.innerSize() == lhsIndex.innerSize() ) ;

		nonZeros.resize( lhsIndex.outerSize() ) ;

	#ifndef BOGUS_DONT_PARALLELIZE
	#pragma omp parallel for
	#endif
		for ( int i = 0 ; i < (int) lhsIndex.outerSize() ; ++i )
		{
			typename MajorIndexType::InnerIterator lhs_it ( lhsIndex, i ) ;
			typename SourceIndexType::InnerIterator rhs_it ( rhsIndex, i ) ;

			NonZero nz ;
			while( lhs_it || rhs_it )
			{
				if( lhs_it && ( !rhs_it || lhs_it.inner() < rhs_it.inner() ) )
				{
					nz.first = lhs_it.inner() ;
					nz.second.first = lhs_it.ptr() ;
					nz.second.second = OtherDerived::InvalidBlockPtr ;
					++ lhs_it ;
				} else if ( rhs_it && ( !lhs_it || rhs_it.inner() < lhs_it.inner() ) ) {
					nz.first = rhs_it.inner() ;
					nz.second.first = InvalidBlockPtr ;
					nz.second.second = rhs_it.ptr() ;
					++ rhs_it ;
				} else {
					nz.first = lhs_it.inner() ;
					nz.second.first = lhs_it.ptr() ;
					nz.second.second = rhs_it.ptr() ;
					++lhs_it ;
					++rhs_it ;
				}

				if( Traits::is_symmetric && nz.first > i )
					break ;

				nonZeros[i].push_back( nz ) ;
			}
		}
	}


	std::vector< BlockPtr > offsets( nonZeros.size() + 1 ) ;

	MajorIndexType resIndex ;
	resIndex.resizeOuter( nonZeros.size() ) ;
    offsets[0] = 0 ;
    for( unsigned i = 0 ; i < nonZeros.size() ; ++i )
    {
        offsets[i+1] = offsets[i] + nonZeros[i].size() ;
    }

	resIndex.reserve( offsets.back() ) ;
	for( unsigned i = 0 ; i < nonZeros.size() ; ++i )
	{
		for( unsigned j = 0 ; j < nonZeros[i].size() ; ++j )
		{
			const BlockPtr ptr = offsets[i]+j ;
			const NonZero &nz = nonZeros[i][j] ;
			resIndex.insertBack( i, nz.first, ptr ) ;
		}
	}

    typename BlockContainerTraits< BlockType >::Type resBlocks( offsets.back() ) ;

	BlockTransposeOption< OtherTraits::is_symmetric, Transpose > rhsGetter ;

	// II - Proper addition

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( int i = 0 ; i < (int) nonZeros.size() ; ++i )
	{
		for( unsigned j = 0 ; j < nonZeros[i].size() ; ++j )
		{
			const BlockPtr ptr = offsets[i]+j ;
			const NonZero &nz = nonZeros[i][j] ;

			BlockType &res = resBlocks[ptr] ;
			const bool afterDiag = ( (bool) Traits::is_col_major )  == ( (bool) OtherTraits::is_col_major )
									? (nz.first > i) : (nz.first < i) ;
			if( nz.second.first == InvalidBlockPtr )
			{
				res = alpha * rhsGetter.get( rhs.block( nz.second.second ), afterDiag ) ;
			} else if( nz.second.second == OtherDerived::InvalidBlockPtr )
			{
				res = block( nz.second.first ) ;
			} else {
				res = block( nz.second.first) + alpha * rhsGetter.get( rhs.block( nz.second.second ), afterDiag ) ;
			}
		}
	}
	resIndex.finalize() ;

	clear() ;
    m_majorIndex.move( resIndex );
	resBlocks.swap( m_blocks ) ;
	m_nBlocks = m_blocks.size() ;
	m_minorIndex.valid = false ;

	Finalizer::finalize( *this ) ;

	return derived() ;
}

template < typename Derived >
template < typename LhsT, typename RhsT >
Derived& SparseBlockMatrixBase<Derived>::operator=( const Addition< LhsT, RhsT > &addition )
{
	typedef Addition< LhsT, RhsT> Add ;

	// WARNING -- Not safe w.r.t aliasing

	Scaling< typename Add::Lhs::ObjectType> lhs ( addition.lhs.object, addition.lhs.scaling ) ;
	*this = lhs ;

	typename Add::Rhs::EvalType rhs = addition.rhs.object.eval() ;
	return add< Add::transposeRhs >( *rhs, addition.rhs.scaling ) ;
}

} //namespace bogus

// operators +, -, *, /

template < typename LhsT, typename RhsT >
bogus::Addition< LhsT, RhsT > operator+ ( const bogus::BlockObjectBase< LhsT >& lhs,
			 const bogus::BlockObjectBase< RhsT > &rhs )
{
	return bogus::Addition< LhsT, RhsT >( lhs.derived(), rhs.derived() ) ;
}

template < typename LhsT, typename RhsT >
bogus::Addition< LhsT, RhsT > operator- ( const bogus::BlockObjectBase< LhsT >& lhs,
			 const bogus::BlockObjectBase< RhsT > &rhs )
{
	return bogus::Addition<  LhsT, RhsT >( lhs.derived(), rhs.derived(), 1, -1 ) ;
}

template < typename Derived >
bogus::Scaling< Derived > operator* ( const bogus::BlockObjectBase< Derived >& lhs,
			 typename Derived::Scalar rhs )
{
	return bogus::Scaling< Derived >( lhs.derived(), rhs ) ;
}

template < typename Derived >
bogus::Scaling< Derived > operator* ( typename Derived::Scalar lhs ,
			 const bogus::BlockObjectBase< Derived >& rhs)
{
	return bogus::Scaling< Derived >( rhs.derived(), lhs ) ;
}

template < typename Derived >
bogus::Scaling< Derived > operator/ ( const bogus::BlockObjectBase< Derived >& lhs,
			 typename Derived::Scalar rhs )
{
	return bogus::Scaling< Derived >( lhs.derived(), 1/rhs ) ;
}

template < typename Derived >
bogus::Scaling< Derived > operator/ ( typename Derived::Scalar lhs ,
			 const bogus::BlockObjectBase< Derived >& rhs)
{
	return bogus::Scaling< Derived >( rhs.derived(), 1/lhs ) ;
}


#endif

