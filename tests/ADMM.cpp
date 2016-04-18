
#include "Core/Block.impl.hpp"
#include "Core/BlockSolvers/GaussSeidel.impl.hpp"
#include "Core/BlockSolvers/ProjectedGradient.impl.hpp"

#include "Extra/SecondOrder.impl.hpp"

#include <Eigen/Cholesky>
#include <gtest/gtest.h>

static void ack( unsigned iter, double err ) {
	std::cout << "A : " << iter << " ==> " << err << std::endl ;
}


#define AMA
#define ACC
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

	bogus::ProjectedGradient< MType > ppg( MassMat ) ;


	MType InvMLambda ;
	InvMLambda.cloneStructure( MassMat ) ;
	for( unsigned i = 0 ; i < MassMat.nBlocks() ; ++i )
	{
		Eigen::MatrixXd Mi =
				lambda * MassMat.block( i ) + Eigen::MatrixXd::Identity( dofs[i], dofs[i] ) ;

		InvMLambda.block( i ) = Mi.inverse() ;
	}

	for( unsigned fp = 0 ; fp < 1 ; ++fp )
	{

//		ut = H*v + w ;
//		for( unsigned i = 0 ; i < n ; ++ i ) {
//			s[3*i] = (1 - rho) * s[3*i] + rho * mu[i] * ut.segment<2>(3*i+1).norm() ;
//		}
//		std::cout << "s  " << s.transpose() << std::endl;

//		v.setZero() ;
//		r.setZero() ;
//		ut = H*v+w +s ;
//		z = ut ;

		Eigen::VectorXd lbda = r, lbdap = r, zp = z,
				vh = v, vp = v;
		double thetap = 1 ;

		double gp = 0 ;

		for( unsigned k = 0 ; k < 1000 ; ++ k ) {

#ifdef CB
			z = -r - tau * ( ut + H*(vh-v) ) ;
			pg.projectOnConstraints( Pkmu, z ) ;

			const double g = (z+r).squaredNorm() ;

			r = -z ;

#endif

#ifdef AMA
			//AMA ( = ADMM w/ lambda = 0)
			// v = argmin ( J(v) - < Hv, r > )
			// v = M^{-1}( H^T r - f )
			v = - InvMassMat * ( gamma/(lambda*lambda)*H.transpose()*r + f ) ;
#else
			// ADMM
			// v = prox_{lambda J }( v - 1/\lambda H^T( Hv - z + u )  )
			tmp = v - gamma/lambda * H.transpose()  * ( ut - z + r  ) ;
			// v = prox_{lambda J }( tmp )
			//   = arginf_( .5 vMv + vf + 1./(2 \lambda) (v-tmp)(v-tmp) )
			//   = arginf_( .5 v(M + 1./(2 lambda) I )v + v'(f - 1/lambda tmp)
			//   = [ M + I/(2lambda) ]^{-1} (tmp/lambda - f)
			v = InvMLambda * ( tmp - lambda * f ) ;
#endif

#ifdef CB
#ifdef ACC
			thetap = 1./sqrt( 1 + 2 * 1.e-8 * gamma ) ;
			gamma = thetap * gamma ;
			tau = tau / thetap ;

			vh = v + thetap * ( v - vp ) ;
			vp = v ;
#else
			vh = v ;
#endif
#endif
			ut = H*v + w ;
			for( unsigned i = 0 ; i < n ; ++ i ) {
//				s[3*i] = (1 - rho) * s[3*i] + rho * mu[i] * ut.segment<2>(3*i+1).norm() ;
			}
			ut += s ;

#ifndef CB

			// z =  P_{Kimu}( H v + u )
			z =  ut + r ;
			ppg.projectOnConstraints( Pkimu, z ) ;

			const double g = std::max( (ut-z).squaredNorm(),
									   (gamma/(lambda*lambda)*H.transpose()*( z - zp)).squaredNorm() )  ;
			if( g < 1.e-24 ) break ;

			zp = z ;

#ifdef ACC
			if( gp < g ) thetap = 1 ;
			gp = g ;

			lbda = r + ut - z ;
			const double beta = bogus::pg_impl::nesterov_inertia( thetap, 0. ) ;
//			double thetan =  (1 + std::sqrt(1 + 4*thetap*thetap))/2 ;
			r = lbda  + beta * ( lbda - lbdap )  ;

			lbdap = lbda ;
//			thetap = theta ;
//			theta  = thetan ;
#else
			// u += Hv - z
			r += ut - z ;
#endif
#endif



	//		std::cout << k << " " << u.transpose() << std::endl ;
			std::cout << k << " g " << g << " \t" << thetap << std::endl ;
	//		std::cout << " J = " << v.transpose() * ( .5 * MassMat * v + f ) << std::endl ;
	//		x = MassMat*v + f ;
	//		double err = ppg.eval( Pkimu, z, x) ;
	//		std::cout << " err " << err << std::endl ;
		}

	std::cout << " g " << (H*v+w+s-z).squaredNorm() << std::endl ;
	std::cout << "Optimal v " << v.transpose() << std::endl ;
	std::cout << "Optimal r " << (-r.transpose() * gamma/lambda/lambda ) << std::endl ;
	std::cout << " J = " << v.transpose() * ( .5 * MassMat * v + f ) << std::endl ;

	}

}
