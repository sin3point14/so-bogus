
/*
 * Any copyright is dedicated to the Public Domain.
 * http://creativecommons.org/publicdomain/zero/1.0/
*/

#include "Core/Block.impl.hpp"
#include "Core/BlockSolvers/GaussSeidel.impl.hpp"
#include "Core/BlockSolvers/LCPLaw.impl.hpp"
#include "Core/BlockSolvers/PyramidLaw.impl.hpp"

#include "Extra/SecondOrder.impl.hpp"

#include "ResidualInfo.hpp"
#include "SmallFrictionPb.hpp"

#include <gtest/gtest.h>

TEST_F( SmallFrictionPb, GaussSeidel )
{
	ResidualInfo ri ;

	// End problem definition

	Eigen::VectorXd b = w - H * ( InvMassMat * f );


	typedef bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::SYMMETRIC > WType ;
	WType W ;
	W = H * InvMassMat * H.transpose() ;


	Eigen::VectorXd x( W.rows() ) ;
	double res = -1 ;

	bogus::GaussSeidel< WType > gs( W ) ;
	ri.bindTo( gs.callback() );

	x.setOnes() ;
	ri.setMethodName( "GS_Hyb" );
	res = gs.solve( bogus::SOCLaw< 3u, double, true, bogus::local_soc_solver::Hybrid >( 2, mu.data() ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;
	ASSERT_TRUE( sol.isApprox( x, 1.e-4 ) ) ;

	x.setOnes() ;
	ri.setMethodName( "GS_PureN" );
	res = gs.solve( bogus::SOCLaw< 3u, double, true, bogus::local_soc_solver::PureNewton >( 2, mu.data() ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;
	ASSERT_TRUE( sol.isApprox( x, 1.e-4 ) ) ;

	x.setOnes() ;
	ri.setMethodName( "GS_PureE" );
	res = gs.solve( bogus::SOCLaw< 3u, double, true, bogus::local_soc_solver::PureEnumerative >( 2, mu.data() ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;
	ASSERT_TRUE( sol.isApprox( x, 1.e-4 ) ) ;

	x.setOnes() ;
	ri.setMethodName( "GS_RevHyb" );
	res = gs.solve( bogus::SOCLaw< 3u, double, true, bogus::local_soc_solver::RevHybrid >( 2, mu.data() ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;
	ASSERT_TRUE( sol.isApprox( x, 1.e-4 ) ) ;

	x.setOnes() ;
	ri.setMethodName( "GS_Def" );
	gs.setTol( 1.e-8 );
	res = gs.solve( bogus::SOC3D( 2, mu.data() ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;

}

TEST( GaussSeidel, LCP )
{
	ResidualInfo ri ;

	typedef Eigen::Matrix< double, 1, 3 > GradBlockT ;
	bogus::SparseBlockMatrix< GradBlockT > H ;

	H.setRows( 4 ) ;
	H.setCols( 3 ) ;

	H.insertBack( 0, 0 ).setOnes() ;
	H.insertBack( 1, 2 ) = Eigen::Vector3d(1,2,3) ;
	H.insertBack( 2, 1 ) = Eigen::Vector3d(-1,2,-3) ;

	H.finalize() ;

	Eigen::VectorXd k ( H.cols() ) ;
	k << -1, -2, 3, 4, -5, 3, -6, 7, -8 ;

	ri.setMethodName( "LCP_GS" ) ;

	Eigen::VectorXd b = H * k ;

	{
		typedef Eigen::Matrix< double, 1, 1 > WBlockT ;
		typedef bogus::SparseBlockMatrix< WBlockT, bogus::SYMMETRIC > WType ;

		WType W = H * H.transpose() ;

		bogus::GaussSeidel< WType > gs( W ) ;
		ri.bindTo( gs.callback() );

		Eigen::VectorXd x ;
		x.setZero( b.rows() ) ;

		gs.setTol(1.e-16) ;
		double res = gs.solve( bogus::LCPLaw< double >(), b, x ) ;

		Eigen::VectorXd y = W*x + b ;

		ASSERT_GT(1.e-16, res) ;
		ASSERT_GT(1.e-8, (x.array()*y.array()).matrix().squaredNorm() ) ;
		ASSERT_LT(-1.e-8, x.minCoeff() ) ;
		ASSERT_LT(-1.e-16, y.minCoeff() ) ;
	}

	{
		typedef double WBlockT ;
		typedef bogus::SparseBlockMatrix< WBlockT, bogus::SYMMETRIC > WType ;

		WType W = H * H.transpose() ;

		bogus::GaussSeidel< WType > gs( W ) ;
		ri.bindTo( gs.callback() );

		Eigen::VectorXd x ;
		x.setZero( b.rows() ) ;

		gs.setTol(1.e-16) ;
		double res = gs.solve( bogus::LCPLaw< double >(), b, x ) ;

		Eigen::VectorXd y = W*x + b ;

		ASSERT_GT(1.e-16, res) ;
		ASSERT_GT(1.e-8, (x.array()*y.array()).matrix().squaredNorm() ) ;
		ASSERT_LT(-1.e-8, x.minCoeff() ) ;
		ASSERT_LT(-1.e-16, y.minCoeff() ) ;
	}
}


TEST( GaussSeidel, Pyramid )
{
	double mu = 0.4 ;

	bogus::PyramidLaw< 3u, double > Plaw( 1, &mu ) ;

	Eigen::Vector3d r  ;
	r.setOnes() ;
	Eigen::Vector3d u  ;
	u.setZero() ;

	ASSERT_LT( 1, Plaw.eval( 0, r, u) ) ;
	Plaw.projectOnConstraint(0,r);
	ASSERT_NEAR( 0, Plaw.eval( 0, r, u), 1.e-16 ) ;

	Eigen::Vector3d r2 = r ;
	Plaw.projectOnConstraint(0,r2);
	ASSERT_TRUE( r2.isApprox(r) ) ;

	Eigen::Vector3d rperp  ;
	rperp[0] = 0 ;
	rperp.tail<2>() = -r.tail<2>() ;

	ASSERT_NEAR( 0, Plaw.eval( 0, r, rperp), 1.e-16 ) ;
	r2 *= -1 ;

	r = r2 ;
	r[0] = 0 ;
	ASSERT_NEAR( 0, Plaw.eval( 0, u, r), 1.e-16 ) ;

	Plaw.projectOnConstraint(0,r2);
	ASSERT_TRUE( r2.isZero() ) ;

	r2[0] = 1 ;
	r = r2 ;
	Plaw.projectOnConstraint(0,r2);
	ASSERT_TRUE( r2.isApprox(r) ) ;

	Eigen::Matrix3d A ;
	A << 3,1,1,1,2,1,1,1,3 ;

	Eigen::Vector3d b ;

	b.setOnes() ;
	Plaw.solveLocal(0,A,b,r,1) ;
	u = A*r + b ;
	ASSERT_NEAR( 0, Plaw.eval( 0, r, u), 1.e-16 ) ;

	b.setConstant(-1) ;
	Plaw.solveLocal(0,A,b,r,1) ;
	u = A*r + b ;
	ASSERT_NEAR( 0, Plaw.eval( 0, r, u), 1.e-6 ) ;

	b[1] = 1 ;
	Plaw.solveLocal(0,A,b,r,1) ;
	u = A*r + b ;
	ASSERT_NEAR( 0, Plaw.eval( 0, r, u), 1.e-6 ) ;
}
