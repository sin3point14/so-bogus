
/*
 * Any copyright is dedicated to the Public Domain.
 * http://creativecommons.org/publicdomain/zero/1.0/
*/

#include "Core/Block.impl.hpp"
#include "Core/BlockSolvers/ProductGaussSeidel.impl.hpp"

#include "Extra/SecondOrder.impl.hpp"

#include "ResidualInfo.hpp"
#include "SmallFrictionPb.hpp"

#include <gtest/gtest.h>

TEST_F( SmallFrictionPb, ProductGaussSeidel )
{
	ResidualInfo ri ;

	// End problem definition

	Eigen::VectorXd b = w - H * ( InvMassMat * f );

	Eigen::VectorXd x( b.rows() ) ;
	double res = -1 ;

	bogus::Coulomb3D law( 2, mu.data()  ) ;

	{
		bogus::ProductGaussSeidel< HType > pgs( H ) ;
		ri.bindTo( pgs.callback() );

		x.setOnes() ;
		ri.setMethodName( "Product_GS" );
		res = pgs.solve( law, b, x ) ;
		ASSERT_LT( res, 1.e-8 ) ;
		Eigen::VectorXd y = H * H.transpose() * x + b ;
		ASSERT_LT( pgs.eval(law, y, x), 1.e-8 ) ;
	}
	{
		bogus::ProductGaussSeidel< HType, MType > pgs( H, InvMassMat ) ;

		x.setOnes() ;
		ri.setMethodName( "Product_GS" );
		res = pgs.solve( law, b, x ) ;
		ASSERT_LT( res, 1.e-8 ) ;
		ASSERT_TRUE( sol.isApprox( x, 1.e-4 ) ) ;
	}
	///////////


	Eigen::VectorXd exp = H*f ;
	Eigen::VectorXd tst  = Eigen::VectorXd::Zero( exp.rows() ) ;

	bogus::Segmenter< 3, Eigen::VectorXd, HType::Index >
			tstSeg( tst, H.rowOffsets() ) ;

	for( int i = 0 ; i < H.rowsOfBlocks() ; ++i ) {
		Eigen::Vector3d lx = tstSeg[i] ;
		H.rowMultiply< false >( i, f, lx ) ;
		tstSeg[i] = lx ;
	}

	EXPECT_EQ( exp, tst ) ;
	Eigen::VectorXd v = H.transpose() * tst ;
	Eigen::VectorXd u = Eigen::VectorXd::Zero( v.rows() ) ;

	for( int i = 0 ; i < H.rowsOfBlocks() ; ++i ) {
		H.colMultiply< true >( i, tstSeg[i], u ) ;
	}

	EXPECT_EQ( v, u ) ;

}
