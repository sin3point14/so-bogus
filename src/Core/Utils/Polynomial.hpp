#ifndef BOGUS_POLYNOMIAL_HPP
#define BOGUS_POLYNOMIAL_HPP

namespace bogus
{

namespace polynomial {

enum RealRootsFilter
{
	StrictlyPositiveRoots,
	StrictlyNegativeRoots,
	AllRoots
} ;

template< unsigned Dimension, typename Scalar >
struct RootsFinder
{
	static unsigned getRealRoots( const Scalar coeffs[Dimension],
								  Scalar realRoots[Dimension],
								  RealRootsFilter filter = AllRoots ) ;
} ;


template< unsigned Dimension, typename Scalar >
unsigned getRealRoots( const Scalar (&coeffs)[Dimension],
					   Scalar (&realRoots)[Dimension],
					   RealRootsFilter filter = AllRoots )
{
	return RootsFinder< Dimension, Scalar >::getRealRoots( coeffs, realRoots, filter ) ;
}


}

}

#endif
