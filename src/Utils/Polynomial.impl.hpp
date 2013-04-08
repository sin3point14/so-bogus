#ifndef BOGUS_POLYNOMIAL_IMPL_HPP
#define BOGUS_POLYNOMIAL_IMPL_HPP

#include "Polynomial.hpp"
#include "NumTraits.hpp"

#include <Eigen/EigenValues>

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
unsigned RootsFinder< Dimension, Scalar>::getRealRoots(
		const Scalar coeffs[Dimension], Scalar realRoots[Dimension],
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

}


}

#endif
