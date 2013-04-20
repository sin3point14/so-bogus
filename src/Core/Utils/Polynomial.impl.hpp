/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver          
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>                       
 *
 * This Source Code Form is subject to the terms of the Mozilla Public 
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_POLYNOMIAL_IMPL_HPP
#define BOGUS_POLYNOMIAL_IMPL_HPP

#include <Eigen/Eigenvalues>

#include "Polynomial.hpp"
#include "NumTraits.hpp"


namespace bogus
{

namespace polynomial {


template< unsigned Dimension, typename Scalar >
struct CompanionMatrix
{
	typedef Eigen::Matrix< Scalar, Dimension, Dimension > BaseType ;
	typedef typename BaseType::MapType ReturnType ;

	static ReturnType get()
	{
		static __thread double s_matrix_data[ Dimension*Dimension ] ;
		static __thread bool s_matrix_initialized = false ;

		ReturnType matrix( s_matrix_data ) ;

		if( !s_matrix_initialized )
		{
			matrix.template block< 1, Dimension-1> ( 0, 0 ).setZero() ;
			matrix.template block< Dimension - 1, Dimension -1 >( 1, 0 ).setIdentity() ;
			s_matrix_initialized = true ;
		}

		return matrix ;
	}
} ;

template< unsigned Dimension, typename Scalar >
unsigned RootsFinder< Dimension, Scalar>::getRealRoots(const Scalar *coeffs, Scalar *realRoots,
		RealRootsFilter filter )
{
	typedef CompanionMatrix< Dimension, Scalar > CM ;
	typename CM::ReturnType matrix = CM::get() ;

	matrix.template block< Dimension, 1 >( 0, Dimension - 1 ) = -Eigen::Matrix< Scalar, Dimension, 1 >::Map( coeffs ) ;
	const typename Eigen::EigenSolver< typename CM::BaseType >::EigenvalueType& ev = matrix.eigenvalues() ;

	unsigned count = 0 ;
	for( unsigned i = 0 ; i < Dimension ; ++i )
	{
		if( NumTraits< Scalar >::isZero( std::imag( ev[i] ) ) )
		{
			const bool discard =
					( filter == StrictlyPositiveRoots && std::real( ev[i] ) <= 0 ) ||
					( filter == StrictlyNegativeRoots && std::real( ev[i] ) >= 0 ) ;
			if( !discard ) realRoots[ count++ ] = std::real( ev[i] ) ;
		}
	}
	return count ;
}

template< typename Scalar >
struct PossiblyDegenerateRootsFinder< 0, Scalar >
{
	static unsigned getRealRoots( const Scalar* coeffs,
								  Scalar* realRoots,
								  RealRootsFilter filter = AllRoots )
	{
		realRoots[0] = 0. ;
		return filter == AllRoots && NumTraits< Scalar >::isZero( coeffs[0] ) ;
	}
} ;

template< unsigned Dimension, typename Scalar >
unsigned PossiblyDegenerateRootsFinder< Dimension, Scalar>::getRealRoots( Scalar *coeffs, Scalar *realRoots,
		RealRootsFilter filter )
{
	if( NumTraits< Scalar >::isZero( coeffs[ Dimension ] ) )
	{
		return PossiblyDegenerateRootsFinder< Dimension - 1, Scalar >::getRealRoots( coeffs, realRoots, filter ) ;
	}
	const Scalar inv = 1./coeffs[Dimension] ;
	for( unsigned k = 0 ; k < Dimension ; ++k )
	{
		coeffs[k] *= inv ;
	}
	return RootsFinder< Dimension, Scalar >::getRealRoots( coeffs, realRoots, filter ) ;
}

} //namespace polynomial

} //namespace bogus

#endif
