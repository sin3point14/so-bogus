
#include "Core/Block.impl.hpp"
#include "Core/BlockSolvers/GaussSeidel.impl.hpp"
#include "Core/BlockSolvers/ProjectedGradient.impl.hpp"

#include "Core/BlockSolvers/ADMM.impl.hpp"

#include "Extra/SecondOrder.impl.hpp"

#include <Eigen/Cholesky>
#include <gtest/gtest.h>

static void ack( unsigned iter, double err ) {
	std::cout << "A : " << iter << " ==> " << err << std::endl ;
}


//#define AMA
//#define ACC
//#define CB


TEST( ADMM, Small )
{

	typedef bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::SYMMETRIC > MType ;
	MType MassMat, InvMassMat ;

	const unsigned n = 2 ;
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

	Eigen::VectorXd f( 6 ) ;
	f << 1, 2, 3, 4, 5, 6 ;
	Eigen::VectorXd w( 6 ) ;
	w << 2, 1, 2, 1, 3, 3 ;
//	w.setZero() ;

	Eigen::VectorXd b = w - H * ( InvMassMat * f );

	double mu[2] = { 0.5, 0.7 } ;

	Eigen::VectorXd sol( 6 ) ;
	sol << 0.0152695, 0.0073010, 0.0022325, 0.0, 0.0, 0.0 ;

	Eigen::VectorXd x( W.rows() ) ;
	double res = -1 ;

//	bogus::GaussSeidel< WType > gs( W ) ;
//	gs.callback().connect( &ack );

//	x.setOnes() ;
//	res = gs.solve( bogus::SOCLaw< 3u, double, true, bogus::local_soc_solver::Hybrid >( 2, mu ), b, x ) ;
//	ASSERT_LT( res, 1.e-8 ) ;
//	ASSERT_TRUE( sol.isApprox( x, 1.e-4 ) ) ;

	bogus::ProjectedGradient< WType > pg( W ) ;
	pg.callback().connect( &ack );
	pg.setTol( 1.e-8 );

	x.setOnes() ;
	res = pg.solve( bogus::SOC3D( 2, mu ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;

	Eigen::VectorXd v = InvMassMat * ( H.transpose() * x - f ) ;
	std::cout << "Optimal v " << v.transpose() << std::endl ;
	std::cout << "Optimal r " << x.transpose() << std::endl ;
	std::cout << " J = " << v.transpose() * ( .5 * MassMat * v + f ) << std::endl ;

	// Test ADMM

	// J = .5 * vMv + f'v
	// g = I(C)(v), C = (Hv \in Kimu)

	const double lambda = 1.e-1 ;
		  double gamma  = 1.e-4 ;
		  double rho    = 1.    ;
		  double tau    = 1.    ;

	v.setZero() ;
	Eigen::VectorXd r = Eigen::VectorXd::Zero( 6 ) ;
	Eigen::VectorXd s = Eigen::VectorXd::Zero( 6 ) ;
	Eigen::VectorXd ut = H*v + w ;
	Eigen::VectorXd z  = ut ;

	Eigen::VectorXd tmp ;
	Eigen::ArrayXd imu = 1./Eigen::ArrayXd::Map( mu, 2 ) ;
	bogus::SOC3D Pkmu ( imu.size(), mu ) ;
	bogus::SOC3D Pkimu( imu.size(), imu.data() ) ;

	MType InvMLambda ;
	InvMLambda.cloneStructure( MassMat ) ;
	for( unsigned i = 0 ; i < MassMat.nBlocks() ; ++i )
	{
		Eigen::MatrixXd Mi =
				MassMat.block( i ) + 1./lambda * Eigen::MatrixXd::Identity( dofs[i], dofs[i] ) ;

		InvMLambda.block( i ) = Mi.inverse() ;
	}

	std::cout << " ============= " << std::endl ;

	r.setOnes() ;
	v = InvMassMat * ( H.transpose() * r - f ) ;

//	bogus::QuadraticProxOp< MType > prox( InvMassMat, 0, f ) ;
	bogus::QuadraticProxOp< MType > prox( InvMLambda, 1.e-1, f ) ;
	bogus::ADMM< HType > admm( H ) ;

	admm.setMaxIters(1000);
	admm.setStepSize(1.e-2);
	admm.callback().connect( &ack );

	admm.solve< bogus::admm::Accelerated >( Pkimu, prox, w, v, r ) ;

}
