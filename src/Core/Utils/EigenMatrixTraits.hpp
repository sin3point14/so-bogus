#ifndef BOGUS_EIGEN_MATRIX_TRAITS_HPP
#define BOGUS_EIGEN_MATRIX_TRAITS_HPP

#include "LinearSolver.hpp"

#include <Eigen/Core>

namespace bogus
{

template< unsigned Dimension, typename ScalarType >
struct MatrixTraits
{
	typedef ScalarType Scalar ;
	typedef Eigen::Matrix< Scalar, Dimension, 1 > Vector ;
	typedef Eigen::Matrix< Scalar, Dimension, Dimension > Matrix ;

	typedef LU  < Eigen::MatrixBase< Matrix > > LUType ;
	typedef LDLT< Eigen::MatrixBase< Matrix > > LDLTType ;

	enum{ dimension = Dimension } ;

	static ScalarType np( const Vector & v )
	{ return v[0] ; }
	static ScalarType& np( Vector & v )
	{ return v[0] ; }

	static typename Vector::template ConstFixedSegmentReturnType< Dimension - 1 >::Type
	tp( const Vector & v )  { return v.template segment< Dimension - 1 >( 1 ) ; }
	static typename Vector::template FixedSegmentReturnType< Dimension - 1 >::Type
	tp(       Vector & v )  { return v.template segment< Dimension - 1 >( 1 ) ; }
} ;

}

#endif
