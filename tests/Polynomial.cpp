
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
	const double mu = 0.6 ;

	 typedef double Scalar ;
	 typedef Eigen::Matrix< Scalar, 2, 1 > Vec2 ;
	 typedef Eigen::Matrix< Scalar, 2, 2 > Mat2 ;

	 const Scalar wN = W(0,0) ;
	 if( wN < bogus::NumTraits< Scalar >::epsilon() )
		 return ; // Could we do something better ?
/*
	 const Vec2 wT = W.block< 2, 1 >( 1, 0 ) ;
	 const Mat2 WT = W.block< 2, 2 >( 1, 1 ) ;

	 const Scalar bN = b[0];
	 const Vec2 bT = b.segment<2>(1) ;

	 const Vec2 wT_wN = wT/wN ;
	 const Mat2 Wbar = WT - wT_wN * wT.transpose() ;
	 const Vec2 bbar = bT/bN - wT_wN ;

	 const Scalar imu2 = 1. / ( mu * mu ) ;

	 const Scalar A = Wbar.trace() - wT.dot( bbar ) ;
	 const Vec2   B ( Wbar(1,1)*bbar[0] - Wbar(1,0)*bbar[1],
					  Wbar(0,0)*bbar[1] - Wbar(0,1)*bbar[0] ) ;
	 const Scalar C = Wbar.determinant() - wT.dot( B ) ;
	 const Scalar D = wN*wN * imu2 ;

	 Scalar coeffs[4] = {
		 C*C - D * B.squaredNorm(),
		 2*( C*A - D * bbar.dot( B ) ),
		 2*C + A*A - D * bbar.squaredNorm(),
		 2*A
	 } ;


	 const Scalar E = - 2 * wN * imu2 ;
	 const Scalar F = bbar.squaredNorm() ;
	 const Scalar G = 2 * B.dot( bbar );
	 const Scalar H = B.squaredNorm() ;

	 coeffs[3] += G * imu2 + E * F ;
	 coeffs[2] += H * imu2 + E * G ;
	 coeffs[1] += E * H ;

	 const Scalar coeff4 = 1 + F * imu2 ;
	 for( unsigned i = 0 ; i < 4 ; ++i )
	 {
		 coeffs[i] /= coeff4 ;
	 }

	 Scalar roots[4] ;
	 const unsigned nRoots =
			 bogus::polynomial::getRealRoots( coeffs, roots, bogus::polynomial::StrictlyPositiveRoots ) ;

	 if( 0 == nRoots ) return ;

	 Scalar alpha = roots[0] ;
	 //Get the minimal one, this is as good an heuristic as any
	 for ( unsigned i = 1 ; i != nRoots ; ++ i )
	 {
		 if( roots[i] < alpha ) alpha = roots[i] ;
	 }*/

	 /*
	 const double alpha = W(0,0) ;



	 std::cout << "Found " << alpha << std::endl ;

	 const Mat2 M = Wbar + alpha * Mat2::Identity() ;
	 r.segment<2>(1) = - bN * M.fullPivLu().solve( bbar ) ;
	 r[0] = r.segment<2>(1).norm() / mu ;*/

	 Eigen::Matrix3d M = W ;
	 M(0,0) = 0 ;
	 M(1,1) += wN ;
	 M(2,2) += wN ;
	 r = - M.fullPivLu().solve( b ) ;

	 std::cout << r.transpose() << std::endl ;
	 Eigen::Vector3d u = W*r + b ;
	 std::cout <<  u[0]  << std::endl ;
	 std::cout <<  u.segment<2>(1) .norm() << std::endl ;
	 std::cout <<  u.segment<2>(1) .norm() * mu << std::endl ;
	 std::cout <<  u.dot( r ) << std::endl ;


}

