
#include "Core/Block.impl.hpp"

#include "Core/BlockSolvers/ADMM.impl.hpp"

#include "Extra/SecondOrder.impl.hpp"

#include "ResidualInfo.hpp"
#include "SmallFrictionPb.hpp"

#include <gtest/gtest.h>

#define TEST_PRIMAL
#define TEST_DUAL


#ifdef TEST_PRIMAL
TEST_F( SmallFrictionPb, ADMM )
{
	ResidualInfo ri ;

	// J = .5 * vMv + f'v
	// g = I(C)(v), C = (Hv \in Kimu)

	Eigen::VectorXd r( H.rows() ), v ( MassMat.rows() ) ;
	double res = -1 ;

	Eigen::ArrayXd imu = 1./mu.array() ;
	bogus::SOC3D Pkimu( imu.size(), imu.data() ) ;

	double lambda = 1.e-1 ;

	MType InvMLambda ;
	InvMLambda.cloneStructure( MassMat ) ;
	for( unsigned i = 0 ; i < MassMat.nBlocks() ; ++i )
	{
		const unsigned dofs = MassMat.blockRows( i ) ;
		Eigen::MatrixXd Mi =
				MassMat.block( i ) + 1./lambda * Eigen::MatrixXd::Identity( dofs, dofs ) ;

		InvMLambda.block( i ) = Mi.inverse() ;
	}

	bogus::ADMM< HType > admm( H ) ;
	admm.setMaxIters(1000);
	admm.setStepSize(1.e-2);
	admm.setTol(1.e-8);

	ri.bindTo( admm.callback() ) ;

	r.setOnes() ;
	v = InvMassMat * ( H.transpose() * r - f ) ;

	ri.setMethodName("AMA");
	bogus::QuadraticProxOp< MType > prox0( InvMassMat, 0, f ) ;
	res = admm.solve< bogus::admm::Accelerated >( Pkimu, prox0, w, v, r ) ;
	ASSERT_LT( res, 1.e-8 ) ;

	r.setOnes() ;
	v = InvMassMat * ( H.transpose() * r - f ) ;

	ri.setMethodName("ADMM");
	bogus::QuadraticProxOp< MType > prox( InvMLambda, lambda, f ) ;
	res = admm.solve< bogus::admm::Accelerated >( Pkimu, prox, w, v, r ) ;
	ASSERT_LT( res, 1.e-8 ) ;
}
#endif

#ifdef TEST_DUAL
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
TEST_F( SmallFrictionPb, DualAMA )
{
	ResidualInfo ri ;


	Eigen::VectorXd r( H.rows() ), v ( MassMat.rows() ) ;
	double res = -1 ;

	bogus::SOC3D Pkmu ( mu.size(), mu.data() ) ;
	bogus::Coulomb3D Cmu ( mu.size(), mu.data() ) ;


	bogus::DualAMA< HType > dama( H ) ;
	ri.bindTo( dama.callback() ) ;

//	dama.callback().connect( &ack );

	dama.setFpStepSize(1.e-1);
	dama.setProjStepSize(1.e0);
	dama.setTol(1.e-16);

	dama.setLineSearchIterations(4);
	dama.setLineSearchOptimisticFactor(1.25);

	ri.setMethodName("DualAMA-Std-Kmu");
	r.setZero() ;
	v = InvMassMat * ( H.transpose() * r - f ) ;
	res = dama.solve< bogus::admm::Standard >( Pkmu, MassMat, f, w, v, r ) ;
	ASSERT_LT( res, 1.e-8 ) ;

	ri.setMethodName("DualAMA-Std-Pkmu");
	r.setZero() ;
	v = InvMassMat * ( H.transpose() * r - f ) ;
	res = dama.solve< bogus::admm::Standard >( Cmu, MassMat, f, w, v, r ) ;
	ASSERT_LT( res, 1.e-8 ) ;
	ASSERT_TRUE( sol.isApprox( r, 1.e-4 ) ) ;

	dama.setLineSearchIterations(0);
	dama.setProjStepSize(1.e-1);

	ri.setMethodName("DualAMA-Acc-Kmu");
	r.setZero() ;
	v = InvMassMat * ( H.transpose() * r - f ) ;
	res = dama.solve< bogus::admm::Accelerated >( Pkmu, MassMat, f, w, v, r ) ;
	ASSERT_LT( res, 1.e-8 ) ;


	dama.setProjStepSize(1.e-2);

	ri.setMethodName("DualAMA-Acc-Cmu");
	r.setZero() ;
	v = InvMassMat * ( H.transpose() * r - f ) ;
	res = dama.solve< bogus::admm::Accelerated >( Cmu, MassMat, f, w, v, r ) ;
	ASSERT_LT( res, 1.e-8 ) ;
	ASSERT_TRUE( sol.isApprox( r, 1.e-4 ) ) ;

//	 std::cout << v.transpose() << std::endl ;
//	 std::cout << r.transpose() << std::endl ;
}
#endif

