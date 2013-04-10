#ifndef BOGUS_NUMTRAITS_HPP
#define BOGUS_NUMTRAITS_HPP

#include <limits>
#include <cmath>

namespace bogus
{

template <typename Scalar>
struct NumTraits
{
	static Scalar epsilon()
	{ return std::numeric_limits< Scalar >::epsilon() ; }
	static bool isZero( Scalar s )
	{ return std::fabs(s) < std::numeric_limits< Scalar >::epsilon() ; }
	static bool isSquareZero( Scalar s )
	{ return s*s < std::numeric_limits< Scalar >::epsilon() ; }
} ;

}

#endif
