
#include "Core/Utils/Polynomial.hpp"
#include "Core/Utils/Polynomial.impl.hpp"

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

TEST( Polynomial, SOQPQuartic )
{

	Eigen::Matrix3d W ;
	W <<  0.01344, -9.421e-07, 0.001486,
		  -9.421e-07, 0.1061, 0.0001733,
		  0.001486, 0.0001733, 0.001442 ;
	Eigen::Vector3d b ;
	b << -0.1458, -0.2484, -0.1515 ;

	Eigen::Vector3d r ;
	Eigen::Vector3d u ;
	const double mu = 0.6 ;

	 typedef double Scalar ;
	 typedef Eigen::Matrix< Scalar, 2, 1 > Vec2 ;
	 typedef Eigen::Matrix< Scalar, 2, 2 > Mat2 ;

	 const Scalar wN = W(0,0) ;
	 ASSERT_FALSE( wN < bogus::NumTraits< Scalar >::epsilon() ) ;

	 const Scalar b0=b[0], b1=b[1], b2=b[2],
			 w00=W(0,0), w01=W(0,1), w02=W(0,2), w11=W(1,1), w12=W(1,2),
			 w22=W(2,2) ;

	 Scalar coeffs[5] ;
	 coeffs[0] =  (b1*w12 - b2*w11)*mu*mu + (b0*w12 + b1*w02 - 2*b2*w01)*mu + b0*w02 - b2*w00 ;
	 coeffs[1] =  2*(b1*w22 - b2*w12)*mu*mu - 2*(b0*w11 - b0*w22 - b1*w01 + b2*w02)*mu - 2*b0*w01 + 2*b1*w00 ;
	 coeffs[2] =  -6*(b0*w12 - b1*w02)*mu ;
	 coeffs[3] =  2*(b1*w22 - b2*w12)*mu*mu + 2*(b0*w11 - b0*w22 - b1*w01 + b2*w02)*mu - 2*b0*w01 + 2*b1*w00 ;
	 coeffs[4] =  -(b1*w12 - b2*w11)*mu*mu + (b0*w12 + b1*w02 - 2*b2*w01)*mu - b0*w02 + b2*w00 ;

	 Scalar roots[4] ;
	 const unsigned nRoots =
			 bogus::polynomial::getRealRoots( coeffs, roots, bogus::polynomial::AllRoots ) ;

	 ASSERT_EQ( 2u, nRoots ) ;

	 for ( unsigned i = 0 ; i != nRoots ; ++ i )
	 {
		 Scalar t = roots[i] ;

		 const Scalar CT = ( 1 - t*t ) / ( 1 + t*t ) ;
		 const Scalar ST = 2*t / ( 1 + t*t ) ;

		 const Eigen::Vector3d dir ( 1, mu*CT, mu*ST ) ;

		 const Scalar den = ( mu * W.col(1) + CT * W.col( 0 )).dot( dir ) ;
		 if( bogus::NumTraits< Scalar >::isZero( den ) )
			 continue ;

		 const Scalar rN = -(CT*b0 + b1*mu)/den ;

		 r = rN * dir ;
		 u = W*r + b ;

		 if( r[0] > 0 && u[0] > 0 )
		 {
			 std::cout << "Found " << t << std::endl ;

			 std::cout << r.transpose() << std::endl ;
			 std::cout <<  u.transpose()  << std::endl ;
			 std::cout <<  u.dot( r ) << std::endl ;

			 break ;
		 }
	 }

	 ASSERT_FLOAT_EQ( u[0], u.segment<2>(1) .norm() * mu ) ;
	 ASSERT_TRUE( u.dot(r) < 1.e-12 ) ;


}

