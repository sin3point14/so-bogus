
/*
 * Any copyright is dedicated to the Public Domain.
 * http://creativecommons.org/publicdomain/zero/1.0/
*/

#include "Core/Block.impl.hpp"
#include "Core/BlockSolvers/ProjectedGradient.impl.hpp"
#include "Extra/SecondOrder.impl.hpp"

#include "ResidualInfo.hpp"
#include "SmallFrictionPb.hpp"

#include <gtest/gtest.h>

TEST_F( SmallFrictionPb, ProjectedGradient )
{
	ResidualInfo ri ;

	Eigen::VectorXd b = w - H * ( InvMassMat * f );

	typedef bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::SYMMETRIC > WType ;
	WType W ;
	W = H * InvMassMat * H.transpose() ;

	Eigen::VectorXd x( W.rows() ) ;
	double res = -1 ;

	bogus::ProjectedGradient< WType > pg( W ) ;
	ri.bindTo( pg.callback() );
	pg.setTol( 1.e-8 );

	x.setOnes() ;
	ri.setMethodName( "PG_DEF" );
	res = pg.solve( bogus::SOC3D( 2, mu.data() ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;

	x.setOnes() ;
	ri.setMethodName( "PG_DESC" );
	res = pg.solve< bogus::projected_gradient::Descent >( bogus::SOC3D( 2, mu.data() ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;

	x.setOnes() ;
	ri.setMethodName( "PG_STD" );
	res = pg.solve< bogus::projected_gradient::Standard >( bogus::SOC3D( 2, mu.data() ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;

	x.setOnes() ;
	ri.setMethodName( "PG_CONJ" );
	res = pg.solve< bogus::projected_gradient::Conjugated >( bogus::SOC3D( 2, mu.data() ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;

	x.setOnes() ;
	ri.setMethodName( "PG_APGD" );
	res = pg.solve< bogus::projected_gradient::APGD >( bogus::SOC3D( 2, mu.data() ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;

	x.setOnes() ;
	ri.setMethodName( "PG_SPG" );
	res = pg.solve< bogus::projected_gradient::SPG >( bogus::SOC3D( 2, mu.data() ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;
}

TEST_F( SmallFrictionPb, MatrixFreeProjectedGradient )
{
	ResidualInfo ri ;

	Eigen::VectorXd b = w - H * ( InvMassMat * f );

	Eigen::VectorXd x( b.rows() ) ;
	double res = -1 ;

	// Without assembling W
	typedef bogus::Product< bogus::Product< HType, MType >, bogus::Transpose< HType > > Prod ;
	Prod prod = H * InvMassMat * H.transpose() ;

	x.setOnes() ;
	bogus::ProjectedGradient< Prod > pgProd( prod ) ;
	ri.bindTo( pgProd.callback() );
	pgProd.setTol( 1.e-8 );

	ri.setMethodName( "PG_PROD" );
	res = pgProd.solve( bogus::SOC3D( 2, mu.data() ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;
}


TEST( ProjectedGradient, Projection )
{

	double mu[2] = { 0.5, 3 } ;
	bogus::SOC3D law( 2, mu ) ;

	Eigen::Vector3d x ;

	x << 2, 0, 1 ;
	law.projectOnConstraint( 0, x );
	EXPECT_EQ( Eigen::Vector3d( 2, 0, 1 ), x ) ;

	x << 1, 3, 0 ;
	law.projectOnConstraint( 1, x );
	EXPECT_EQ( Eigen::Vector3d( 1, 3, 0 ), x ) ;
	law.projectOnConstraint( 0, x );
	EXPECT_TRUE( Eigen::Vector3d( 2, 1, 0 ).isApprox( x ) ) ;

	x << -1, 7, 0 ;
	law.projectOnConstraint( 0, x );
	EXPECT_TRUE( Eigen::Vector3d( 2, 1, 0 ).isApprox( x ) ) ;

	x << -1, .3, 0 ;
	law.projectOnConstraint( 1, x );
	EXPECT_TRUE( x.isZero() ) ;
}

