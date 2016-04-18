
/*
 * Any copyright is dedicated to the Public Domain.
 * http://creativecommons.org/publicdomain/zero/1.0/
*/

#include "Core/Block.impl.hpp"
#include "Core/BlockSolvers/GaussSeidel.impl.hpp"
#include "Core/BlockSolvers/ProjectedGradient.impl.hpp"
#include "Extra/SecondOrder.impl.hpp"
#include "Core/BlockSolvers/LCPLaw.impl.hpp"
#include "Core/BlockSolvers/PyramidLaw.impl.hpp"

#include <Eigen/LU>
#include <Eigen/Cholesky>

#include <gtest/gtest.h>

struct ResidualInfo {

	static const char* s_meth ;

	static void ack( unsigned iter, double err ) {
		std::cout << s_meth << ": " << iter << " ==> " << err << std::endl ;
	}
};
const char* ResidualInfo::s_meth = "" ;

TEST( GaussSeidel, Small )
{

	bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::UNCOMPRESSED > MassMat ;
	typedef bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::UNCOMPRESSED > MType ;
	MType InvMassMat ;

	const unsigned dofs[2] = { 4, 2 } ;

	MassMat.setRows( 2, dofs ) ;
	MassMat.setCols( 2, dofs ) ;

	MassMat.insertBackAndResize( 0, 0 ) << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 ;
	MassMat.insertBackAndResize( 1 ,1 ) << 2, 0, 0, 2 ;

	MassMat.finalize() ;

	InvMassMat.cloneStructure( MassMat ) ;
	for( unsigned i = 0 ; i < MassMat.nBlocks() ; ++i )
	{
		InvMassMat.block( i ) = MassMat.block( i ).inverse() ;
	}


	typedef Eigen::Matrix< double, 3, Eigen::Dynamic > GradBlockT ;
	typedef bogus::SparseBlockMatrix< GradBlockT > HType ;
	HType H ;

	H.setCols( 2, dofs ) ;
	H.setRows( 2 ) ;
	H.insertBackAndResize( 0, 0 ) << 3, 5, 7, -1, 2, 4, 5, 4, 1, 6, 9, 1 ;
	H.insertBackAndResize( 0, 1 ) << -3, -6, -6, -5, -5, -8 ;
	H.insertBackAndResize( 1, 0 ) << 3, 7, 1, -1, 3, 6, 3, 2, 4, 4, 7, -1 ;

	GradBlockT H2B( 3, dofs[0] ) ;
	H2B << 6, 7, 10, 9, 5, 8, 11, 8, 4, 9, 12, 7 ;
	H.block( 2 ) -= H2B ;

	Eigen::Matrix3d E ;
	E << 0, 0, 1, 1, 0, 0, 0, 1, 0 ;

	H.block( 2 ) = E.transpose() * H.block( 2 ) ;
	H.finalize() ;

	//std::cout << MassMat << std::endl ;
	//	std::cout << InvMassMat << std::endl ;
	//	std::cout << H << std::endl ;

	typedef bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::SYMMETRIC > WType ;
	WType W ;
	W = H * InvMassMat * H.transpose() ;

	//	std::cout << W << std::endl ;
	//	W.cacheTranspose();

	Eigen::VectorXd f( 6 ) ;
	f << 1, 2, 3, 4, 5, 6 ;
	Eigen::VectorXd w( 6 ) ;
	w << 2, 1, 2, 1, 3, 3 ;

	Eigen::VectorXd b = w - H * ( InvMassMat * f );

	double mu[2] = { 0.5, 0.7 } ;

	Eigen::VectorXd sol( 6 ) ;
	sol << 0.0152695, 0.0073010, 0.0022325, 0.0, 0.0, 0.0 ;

	Eigen::VectorXd x( W.rows() ) ;
	double res = -1 ;

	bogus::GaussSeidel< WType > gs( W ) ;
	gs.callback().connect( &ResidualInfo::ack );

	x.setOnes() ;
	ResidualInfo::s_meth = "GS_Hyb" ;
	res = gs.solve( bogus::SOCLaw< 3u, double, true, bogus::local_soc_solver::Hybrid >( 2, mu ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;
	ASSERT_TRUE( sol.isApprox( x, 1.e-4 ) ) ;

	x.setOnes() ;
	ResidualInfo::s_meth = "GS_PureN" ;
	res = gs.solve( bogus::SOCLaw< 3u, double, true, bogus::local_soc_solver::PureNewton >( 2, mu ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;
	ASSERT_TRUE( sol.isApprox( x, 1.e-4 ) ) ;

	x.setOnes() ;
	ResidualInfo::s_meth = "GS_PureE" ;
	res = gs.solve( bogus::SOCLaw< 3u, double, true, bogus::local_soc_solver::PureEnumerative >( 2, mu ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;
	ASSERT_TRUE( sol.isApprox( x, 1.e-4 ) ) ;

	x.setOnes() ;
	ResidualInfo::s_meth = "GS_RevHyb" ;
	res = gs.solve( bogus::SOCLaw< 3u, double, true, bogus::local_soc_solver::RevHybrid >( 2, mu ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;
	ASSERT_TRUE( sol.isApprox( x, 1.e-4 ) ) ;


	x.setOnes() ;
	gs.setTol( 1.e-8 );
	res = gs.solve( bogus::SOC3D( 2, mu ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;

	bogus::ProjectedGradient< WType > pg( W ) ;
	pg.callback().connect( &ResidualInfo::ack );
	pg.setTol( 1.e-8 );

	x.setOnes() ;
	ResidualInfo::s_meth = "PG_DEF" ;
	res = pg.solve( bogus::SOC3D( 2, mu ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;

	x.setOnes() ;
	ResidualInfo::s_meth = "PG_DESC" ;
	res = pg.solve< bogus::projected_gradient::Descent >( bogus::SOC3D( 2, mu ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;

	x.setOnes() ;
	ResidualInfo::s_meth = "PG_STD" ;
	res = pg.solve< bogus::projected_gradient::Standard >( bogus::SOC3D( 2, mu ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;

	x.setOnes() ;
	ResidualInfo::s_meth = "PG_CONJ" ;
	res = pg.solve< bogus::projected_gradient::Conjugated >( bogus::SOC3D( 2, mu ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;

	x.setOnes() ;
	ResidualInfo::s_meth = "PG_APGD" ;
	res = pg.solve< bogus::projected_gradient::APGD >( bogus::SOC3D( 2, mu ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;

	x.setOnes() ;
	ResidualInfo::s_meth = "PG_SPG" ;
	res = pg.solve< bogus::projected_gradient::SPG >( bogus::SOC3D( 2, mu ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;

	// Without assembling W
	typedef bogus::Product< bogus::Product< HType, MType >, bogus::Transpose< HType > > Prod ;
	Prod prod = H * InvMassMat * H.transpose() ;

	x.setOnes() ;
	bogus::ProjectedGradient< Prod > pgProd( prod ) ;
	ResidualInfo::s_meth = "PG_PROD" ;
	pgProd.callback().connect( &ResidualInfo::ack );
	pgProd.setTol( 1.e-8 );

	res = pgProd.solve( bogus::SOC3D( 2, mu ), b, x ) ;
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

TEST( GaussSeidel, LCP )
{
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

	ResidualInfo::s_meth = "LCP_GS" ;

	Eigen::VectorXd b = H * k ;

	{
		typedef Eigen::Matrix< double, 1, 1 > WBlockT ;
		typedef bogus::SparseBlockMatrix< WBlockT, bogus::SYMMETRIC > WType ;

		WType W = H * H.transpose() ;

		bogus::GaussSeidel< WType > gs( W ) ;
		gs.callback().connect( &ResidualInfo::ack );

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
		gs.callback().connect( &ResidualInfo::ack );

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
