
#include "Utils/Polynomial.hpp"
#include "Utils/Polynomial.impl.hpp"

#include <gtest/gtest.h>

TEST( Polynomial, Quadratic )
{
	double c[2] = {-1, 0} ;
	double x[2] ;

	unsigned nRoots = bogus::polynomial::getRealRoots( c, x ) ;
	EXPECT_EQ( 2u, nRoots ) ;
	nRoots = bogus::polynomial::getRealRoots( c, x, bogus::polynomial::StrictlyPositiveRoots ) ;
	EXPECT_EQ( 1u, nRoots ) ;
	EXPECT_DOUBLE_EQ(  1., x[0] ) ;
	nRoots = bogus::polynomial::getRealRoots( c, x, bogus::polynomial::StrictlyNegativeRoots ) ;
	EXPECT_EQ( 1u, nRoots ) ;
	EXPECT_DOUBLE_EQ(  -1., x[0] ) ;

	c[0] = 1 ;
	nRoots = bogus::polynomial::getRealRoots( c, x ) ;
	EXPECT_EQ( 0u, nRoots ) ;
	c[1] = -2 ;
	nRoots = bogus::polynomial::getRealRoots( c, x ) ;
	EXPECT_EQ( 2u, nRoots ) ;
	EXPECT_DOUBLE_EQ(  1., x[0] ) ;
	EXPECT_DOUBLE_EQ(  1., x[1] ) ;
}
