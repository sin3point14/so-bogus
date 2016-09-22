/*
 * Any copyright is dedicated to the Public Domain.
 * http://creativecommons.org/publicdomain/zero/1.0/
*/

#include <iostream>

#include "Core/Block.impl.hpp"
#include "Core/Block.io.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <gtest/gtest.h>

class FlatSBM : public ::testing::Test {
protected:

	virtual void SetUp() {

		//Eigen::initParallel()
		Eigen::MatrixXf A = Eigen::MatrixXf::Zero(1,1); A = A*A;

		EXPECT_TRUE( bogus::IsTransposable< BlockT >::Value ) ;
		sbm.setRows( 5, 3 ) ;
		sbm.setCols( 2, 4 ) ;

		sbm.insertBackAndResize( 0, 1 ).setOnes()  ;
		sbm.insertBack( 3, 0 ) = 2 * BlockT::Ones(3,4) ;
		sbm.insertBack( 3, 1 ).setConstant(3) ;

		sbm.finalize();

		rhs.resize( sbm.cols() ) ;
		rhs.setOnes() ;

		expected_1.resize(15) ;  // sbm * rhs
		expected_2.resize(8) ;   //
		expected_3.resize(15) ;

		expected_1 << 4, 4, 4, 0, 0, 0, 0, 0, 0, 20, 20, 20, 0, 0, 0 ;
		expected_2 << 120, 120, 120, 120, 192, 192, 192, 192 ;
		expected_3 << 768, 768, 768, 0, 0, 0, 0, 0, 0, 3264, 3264, 3264, 0, 0, 0 ;

		ASSERT_EQ( 2u, sbm.blockPtr( 3, 1 ) ) ;
		ASSERT_EQ( sbm.InvalidBlockPtr, sbm.blockPtr( 0, 0 ) ) ;
		ASSERT_EQ( sbm.InvalidBlockPtr, sbm.blockPtr( 1, 1 ) ) ;
	}

	typedef Eigen::MatrixXd BlockT ;
	typedef bogus::FlatSparseBlockMatrix< BlockT > SBMT ;
	SBMT sbm ;
	Eigen::VectorXd rhs  ;

	Eigen::VectorXd expected_1 ; // sbm * rhs
	Eigen::VectorXd expected_2 ; // sbm.transpose() * (sbm * rhs)
	Eigen::VectorXd expected_3 ; // sbm * ( sbm.transpose() * (sbm * rhs) )

} ;

TEST_F(FlatSBM, Construction)
{
	Eigen::VectorXd res ( sbm.rows() ) ;
	res.setZero() ;
	sbm.multiply< false >( rhs,res ) ;

	EXPECT_EQ( expected_1, res ) ;

	bogus::FlatSparseBlockMatrix< BlockT > ss = 2*sbm ;
	ss.multiply< false >( rhs,res ) ;
	EXPECT_EQ( 2*expected_1, res ) ;

	bogus::FlatSparseBlockMatrix< BlockT > aa = ss + sbm ;
	aa.multiply< false >( rhs,res ) ;
	EXPECT_EQ( 3*expected_1, res ) ;

	bogus::FlatSparseBlockMatrix< BlockT > mm = aa.transpose()*sbm ;
	res.resize( mm.rows() ) ;
	mm.multiply< false >( rhs,res ) ;
	EXPECT_EQ( aa.transpose()*expected_1, res ) ;
}

TEST_F( FlatSBM, Convert )
{
	Eigen::VectorXd res ( sbm.rows() ) ;
	res.setZero() ;

	Eigen::SparseMatrix< double, Eigen::RowMajor > csr ;
	bogus::convert( sbm, csr ) ;

	SBMT sbm2 ;
	bogus::convert( csr, sbm2, 3, 4 ) ;

	sbm2.multiply< false >( rhs,res ) ;
	EXPECT_EQ( expected_1, res ) ;
}
