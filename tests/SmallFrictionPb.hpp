/*
 * Any copyright is dedicated to the Public Domain.
 * http://creativecommons.org/publicdomain/zero/1.0/
*/

#ifndef D6_TEST_SMALL_FRICTION_PB
#define D6_TEST_SMALL_FRICTION_PB

#include <gtest/gtest.h>

#include "Core/Block.impl.hpp"

#include <Eigen/LU>
#include <Eigen/Cholesky>

class SmallFrictionPb : public ::testing::Test {
protected:

	virtual void SetUp() {

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

		//	std::cout << W << std::endl ;
		//	W.cacheTranspose();

		f.resize( 6 ) ;
		f << 1, 2, 3, 4, 5, 6 ;
		w.resize( 6 ) ;
		w << 2, 1, 2, 1, 3, 3 ;

		mu.resize(2) ;
		mu <<  0.5, 0.7 ;

		sol.resize( 6 ) ;
		sol << 0.0152695, 0.0073010, 0.0022325, 0.0, 0.0, 0.0 ;
	}

	typedef bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::UNCOMPRESSED > MType ;
	typedef Eigen::Matrix< double, 3, Eigen::Dynamic > GradBlockT ;
	typedef bogus::SparseBlockMatrix< GradBlockT > HType ;

	MType MassMat ;
	MType InvMassMat ;

	HType H ;
	Eigen::VectorXd f ;
	Eigen::VectorXd w ;

	Eigen::VectorXd mu ;

	Eigen::VectorXd sol ;

	// virtual void TearDown() {}

};


#endif
