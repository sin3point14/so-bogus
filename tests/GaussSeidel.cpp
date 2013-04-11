
#include "Core/Block.impl.hpp"
#include "Core/GaussSeidel.impl.hpp"
#include "Core/SecondOrder.impl.hpp"

#include <Eigen/LU>
#include <Eigen/Cholesky>

#include <gtest/gtest.h>

TEST( GaussSeidel, Small )
{

	bogus::SparseBlockMatrix< Eigen::MatrixXd > MassMat ;
	bogus::SparseBlockMatrix< Eigen::MatrixXd > InvMassMat ;

	std::vector < unsigned > dofs ;
	dofs.push_back( 4 ) ;
	dofs.push_back( 2 ) ;

	MassMat.setRows( dofs ) ;
	MassMat.setCols( dofs ) ;

	MassMat.insertBackAndResize( 0, 0 ) << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 ;
	MassMat.insertBackAndResize( 1 ,1 ) << 2, 0, 0, 2 ;

	MassMat.finalize() ;

	InvMassMat.cloneStructure( MassMat ) ;
	for( unsigned i = 0 ; i < MassMat.nBlocks() ; ++i )
	{
		InvMassMat.block( i ) = MassMat.block( i ).inverse() ;
	}


	typedef Eigen::Matrix< double, 3, Eigen::Dynamic > GradBlockT ;
	bogus::SparseBlockMatrix< GradBlockT > H ;

	H.setCols( dofs ) ;
	H.setRows( 2, 3 ) ;
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

	typedef Eigen::Matrix< double, Eigen::Dynamic, 3 > GradBlockTransT ;
	bogus::SparseBlockMatrix< GradBlockTransT, bogus::flags::COL_MAJOR >
			MInvHt = InvMassMat * H.transpose() ;

//	std::cout << MInvHt << std::endl ;

	typedef bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::SYMMETRIC | bogus::flags::COMPRESSED > WType ;
	WType W = H * MInvHt ;

//	std::cout << W << std::endl ;
	W.cacheTranspose();

	Eigen::VectorXd f( 6 ) ;
	f << 1, 2, 3, 4, 5, 6 ;
	Eigen::VectorXd w( 6 ) ;
	w << 2, 1, 2, 1, 3, 3 ;

	Eigen::VectorXd b = w - H * ( InvMassMat * f );

	bogus::GaussSeidel< WType > gs( W ) ;

	std::vector< double > mu ;
	mu.resize(2) ;
	mu[0] = 0.5 ;
	mu[1] = 0.7 ;

	Eigen::VectorXd sol( 6 ) ;
	sol << 0.0152695, 0.0073010, 0.0022325, 0.0, 0.0, 0.0 ;

	Eigen::VectorXd x( W.rows() ) ;

	x.setOnes() ;
	double res = gs.solve( bogus::Coulomb3D( mu ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;
	ASSERT_TRUE( sol.isApprox( x, 1.e-4 ) ) ;

	x.setOnes() ;
	res = gs.solve( bogus::SOCLaw< Eigen::Matrix3d, true, bogus::local_soc_solver::PureNewton >( mu ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;
	ASSERT_TRUE( sol.isApprox( x, 1.e-4 ) ) ;

	x.setOnes() ;
	res = gs.solve( bogus::SOCLaw< Eigen::Matrix3d, true, bogus::local_soc_solver::PureEnumerative >( mu ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;
	ASSERT_TRUE( sol.isApprox( x, 1.e-4 ) ) ;

	x.setOnes() ;
	res = gs.solve( bogus::SOCLaw< Eigen::Matrix3d, true, bogus::local_soc_solver::RevHybrid >( mu ), b, x ) ;
	ASSERT_LT( res, 1.e-8 ) ;
	ASSERT_TRUE( sol.isApprox( x, 1.e-4 ) ) ;
}
