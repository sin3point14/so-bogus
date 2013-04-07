#ifndef BOGUS_NUMTRAITS_HPP
#define BOGUS_NUMTRAITS_HPP

#include <limits>
#include <cmath>

#include <Eigen/Core>

namespace bogus
{

template <typename Scalar>
struct NumTraits
{
	static Scalar epsilon()
	{ return std::numeric_limits< Scalar >::epsilon() ; }
	static bool isZero( Scalar s )
	{ return std::fabs(s) < std::numeric_limits< Scalar >::epsilon() ; }
} ;

template< unsigned Dimension, typename ScalarType >
struct MatrixTraits
{
        typedef ScalarType Scalar ;
	typedef Eigen::Matrix< Scalar, Dimension, 1 > Vector ;
	typedef Eigen::Matrix< Scalar, Dimension, Dimension > Matrix ;

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
