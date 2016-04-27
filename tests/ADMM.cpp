
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

		  double lambda = 1.e-1 ;
		  double gamma  = 1.e-4 ;
//		  double rho    = 1.    ;
//		  double tau    = 1.    ;

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
	admm.setTol(1.e-8);
	admm.callback().connect( &ack );

	res = admm.solve< bogus::admm::Accelerated >( Pkimu, prox, w, v, r ) ;
	ASSERT_LT( res, 1.e-8 ) ;

	std::cout << " ============= " << std::endl ;
	std::cout << " TEST DUAL " << std::endl ;

	// Rationale
	/*
	 * min J(r) + IKmu(r)
	 *
	 * J(r) =.5  r W r + ( w - HM^{-1} f ) ' r
	 *
	 * J_0 = .5 r' H M^-1 H' r - r' H M^{-1} f = .5 x' M^-1 x - x' M^-1 f
	 * G   = w'r + IKmu(r)
	 *
	 * min J(x) + G(r)  with x = H'r
	 *
	 * prox_{l, G}( y ) = argmin_r  I_c(r) + r'w + 1/2l || r - y ||^2
	 *                  = argmin_r  I_c(r) + 1/2l r^2 - 1/l r'( y - l w )
	 *                  = argmin_r  I_c(r) + 1/2l || r - ( y - l w ) || ^2
	 *                  = prox_{l, I_c} ( y - l w )
	 *
	 * AMA -- Goldstein
	 * H = J_0  ; G = G ; A = Id ; B = -H'
	 *
	 * x = argmin J_0(x) + <v, -x>  = argmin .5 x' M^-1 x - x'(M^-1 f + v) = f + M v
	 * r
	 *   = argmin  G(r) + <v, H'r> + gamma/2 || H'r - x  ||
	 *   = argmin Ic(r) + <H v + w, r> + gamma/2 || H'r - x  ||
	 *   = argmin Ic(r) + <H v + w, r> + gamma < H H'r^k - H x, r > + 1/2l || r - r^k ||
	 *   = argmin Ic(r) + <H v + w + gamma ( H H'r^k - H x ), r > + 1/2l || r - r^k ||
	 *   = prox_{Ic,l}(  r^k - l (H v + w) - l gamma ( H H'r^k - H x ) )
	 *   = Pi_Kmu (  r - l gamma H ( H'r - x ) - l ( Hv + w ) )
	 *
	 * v = v + gamma ( H' r - x )
	 *
	 * ADMM -- Boyd
	 * min G(r) + J_0( H'r )
	 * f = G ; g = J ; A = H'
	 *
	 * r = prox_g,G ( r - g/l H ( H'r - z + u ) )
	 *   = Pi_Kmu (  r - g/l H ( H'r - z + u ) - g w )
	 * z = prox_{l,J_0} ( u - H'r )
	 *   = argmin .5 z M^-1 z - z' M^-1 f + 1/2l || z - (u - H'r) ||^2
	 *   = argmin .5 z ( M^-1 +I/2l ) z - z' ( M^-1 f + (u - H'r) )
	 *   = ( M^-1 + 1/2l )^-1  ( M^-1 f + 1/l (u - H'r) )
	 * u = u + H'r - z
	 *
	 * avec u/l = v, u = l v
	 *
	 * r = prox_g,G ( r - g H ( (H'r - z)/l + v ) )
	 *   = Pi_Kmu (  r - 1/l H ( H'r - z ) - g( Hv + w ) )
	 *
	 * z = ( M^-1 + 1/2l )^-1  ( M^-1 f + v - H'r/l ) )
	 *   = f + Mv ( pour 1/lambda = 0 )
	 *
	 * v = v + 1/l ( H'r - z )
	 *
	 * Chambolle Pock
	 * min F(H'r) + G(r)
	 * F = J_0 ; G = G ; K = H'
	 *
	 * v = prox_F*,s ( v + s H' r_ )
	 * r = prox_G,l  ( r - l H v ) = prox_Ic,l ( r - l H v - l w )
	 *   = PKmu ( r - (Hv + w) )
	 * r_ = r + t ( r - r_{k-1} )
	 *
	 *
	 * prox_F*,s = v - s prox_{F, 1/s}( v/s )
	 *           = v - s argmin( J_0(x) + s/2 || x - v/s  || )
	 *           = v - s argmin( .5 x'M^-1 x - x'(M^-1 f) + s/2 || x - v/s  || )
	 *           = v - s argmin( .5 x' ( M^-1 + s I ) x - x'(M^-1 f + v) )
	 *           = v - s ( M^-1 + s I )^-1 ( x'(M^-1 f + v) )
	 *
	 */

	 std::cout << v.transpose() << std::endl ;
#ifdef TEST_INLINE
	 r.setZero() ;
	v = InvMassMat * ( H.transpose() * r - f ) ;

	 lambda = 1.e-1 ;
	 gamma  = 2.e-1 ;

	 for( unsigned k = 0 ; k < 300 ; ++k ) {
		 x = MassMat * v + f ;
		 ut = H*v + w ;

		 r = r - lambda * ( gamma * H * ( H.transpose() * r - x ) + ut ) ;
		 pg.projectOnConstraints( Pkmu, r ) ;

		 std::cout << (H.transpose() * r -x).squaredNorm() << std::endl ;
//		 std::cout << k << "   " << pg.eval( Pkmu, ut, r ) << std::endl;

		 v += gamma * ( H.transpose() * r - x ) ;

	 }

	 std::cout << " ============= " << std::endl ;
#endif
	r.setOnes() ;
	v = InvMassMat * ( H.transpose() * r - f ) ;

	bogus::DualAMA< HType > dama( H ) ;
	dama.callback().connect( &ack );

	dama.setFpStepSize(1.e-1);
	dama.setProjStepSize(1.e-1);
	dama.setTol(1.e-8);

	res = dama.solve< bogus::admm::Accelerated >( Pkmu, MassMat, f, w, v, r ) ;
	ASSERT_LT( res, 1.e-8 ) ;

	 std::cout << v.transpose() << std::endl ;
}
