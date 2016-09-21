
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
//	ri.setVerbose(true);

	// End problem definition

	Eigen::VectorXd b = w - H * ( InvMassMat * f );

	Eigen::VectorXd x( b.rows() ) ;
	double res = -1 ;

	bogus::Coulomb3D law( 2, mu.data()  ) ;

	{
		ASSERT_TRUE( bogus::ProductGaussSeidel< HType >::has_trivial_diagonal ) ;
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
		ASSERT_FALSE( (bogus::ProductGaussSeidel< HType, MType >::has_trivial_diagonal) ) ;
		bogus::ProductGaussSeidel< HType, MType, false > pgs( H, InvMassMat ) ;
		ri.bindTo( pgs.callback() );

		x.setOnes() ;
		ri.setMethodName( "Product_GS_D" );
		res = pgs.solve( law, b, x ) ;
		ASSERT_LT( res, 1.e-8 ) ;
		ASSERT_TRUE( sol.isApprox( x, 1.e-4 ) ) ;
	}
	{
		bogus::ProductGaussSeidel< HType, MType > pgs( H, InvMassMat ) ;
		ri.bindTo( pgs.callback() );

		x.setOnes() ;
		ri.setMethodName( "Product_GS_D_precomp" );
		res = pgs.solve( law, b, x ) ;
		ASSERT_LT( res, 1.e-8 ) ;
		ASSERT_TRUE( sol.isApprox( x, 1.e-4 ) ) ;
	}

}
