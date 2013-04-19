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
	static unsigned getRealRoots( const Scalar* coeffs,
								  Scalar* realRoots,
								  RealRootsFilter filter = AllRoots ) ;
} ;

template< unsigned Dimension, typename Scalar >
struct PossiblyDegenerateRootsFinder
{
	static unsigned getRealRoots( Scalar* coeffs,
								  Scalar* realRoots,
								  RealRootsFilter filter = AllRoots ) ;
} ;

template< unsigned Dimension, typename Scalar >
unsigned getRealRoots( const Scalar (&coeffs)[Dimension],
					   Scalar (&realRoots)[Dimension],
					   RealRootsFilter filter = AllRoots )
{
	return RootsFinder< Dimension, Scalar >::getRealRoots( coeffs, realRoots, filter ) ;
}

template< unsigned Dimension, typename Scalar >
unsigned getRealRoots( Scalar (&coeffs)[Dimension+1],
					   Scalar (&realRoots)[Dimension],
					   RealRootsFilter filter = AllRoots )
{
	return PossiblyDegenerateRootsFinder< Dimension, Scalar >::getRealRoots( coeffs, realRoots, filter ) ;
}

} //namespace polynomial

} //namespace bogus

#endif
