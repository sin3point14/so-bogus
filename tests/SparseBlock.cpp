#include <iostream>

#include "Block.hpp"

#include <Eigen/Core>

#include <gtest/gtest.h>


TEST( SparseBlock, MatrixVector )
{
	Eigen::VectorXd expected_1(15), expected_2(8), expected_3(15) ;
	expected_1 << 4, 4, 4, 0, 0, 0, 0, 0, 0, 20, 20, 20, 0, 0, 0 ;
	expected_2 << 120, 120, 120, 120, 192, 192, 192, 192 ;
	expected_3 << 768, 768, 768, 0, 0, 0, 0, 0, 0, 3264, 3264, 3264, 0, 0, 0 ;

	typedef Eigen::MatrixXd BlockT ;
	bogus::SparseBlockMatrix< BlockT > sbm ;
	sbm.setRows( 5, 3 ) ;
	sbm.setCols( 2, 4 ) ;

	sbm.insertBack( 0, 1 ) = BlockT::Ones( 3,4 ) ;
	sbm.insertBack( 3, 0 ) = 2 * BlockT::Ones( 3,4 ) ;
	sbm.insertBack( 3, 1 ) = 3 * BlockT::Ones( 3,4 ) ;

	sbm.finalize();
	//std::cout << sbm << std::endl ;

	Eigen::VectorXd rhs ( sbm.cols() ) ;
	Eigen::VectorXd res ( sbm.rows() ) ;

	rhs.setOnes() ;
	res.setZero() ;
	sbm.multiply( rhs, res ) ;

	EXPECT_EQ( expected_1, res ) ;
	EXPECT_EQ( expected_1, sbm*rhs ) ;

	rhs.setZero() ;
	sbm.multiply( res, rhs, true ) ;

	EXPECT_EQ( expected_2, rhs ) ;
	EXPECT_EQ( expected_2, ( res.transpose() * sbm ).transpose() ) ;
	EXPECT_EQ( expected_2, ( sbm.transpose() * res )  ) ;
	EXPECT_EQ( expected_3, ( rhs.transpose() * sbm.transpose() ).transpose()  ) ;

	sbm.cacheTranspose();
	EXPECT_EQ( expected_2, ( res.transpose() * sbm ).transpose() ) ;
	EXPECT_EQ( expected_2, ( sbm.transpose() * res )  ) ;

	bogus::SparseBlockMatrix< Eigen::MatrixXd > copyofsbm ;
	copyofsbm.cloneStructure( sbm ) ;
	copyofsbm = sbm ;
	sbm = copyofsbm.transpose() ;
	//std::cout << sbm << std::endl;

	EXPECT_EQ( expected_2, ( sbm * res )  ) ;

	bogus::SparseBlockMatrix< Eigen::MatrixXd > prod( sbm.transpose() * copyofsbm );
}

TEST( SparseBlock, Symmetric )
{
	Eigen::VectorXd expected_1(9), expected_2(9) ;
	expected_1 << 9, 7, 5, 2, 4, 6, 9, 9, 9 ;
	expected_2 << 6, 4, 2, 2, 4, 6, 0, 0, 0 ;

	Eigen::VectorXd res, rhs ;

	bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::BlockMatrixFlags::SYMMETRIC | bogus::BlockMatrixFlags::COMPRESSED > ssbm ;
	ssbm.setRows( 3 ) ;

	ssbm.insertBack( 0, 0 ) = Eigen::Matrix3d::Ones() ;
	Eigen::Matrix3d secBlock ;
	secBlock << 2, 0, 0, 2, 2, 0, 2, 2, 2 ;
	ssbm.insertBack( 1, 0 ) = secBlock ;
	ssbm.insertBack( 2, 2 ) = 3*Eigen::Matrix3d::Ones() ;

	ssbm.finalize() ;

	rhs.resize( ssbm.rows() ) ;
	rhs.setOnes() ;
	EXPECT_EQ( rhs.rows(), 9 ) ;

	EXPECT_EQ( expected_1, ssbm * rhs ) ;
	EXPECT_EQ( expected_1, ssbm.transpose() * rhs ) ;

	ssbm.cacheTranspose() ;

	EXPECT_EQ( expected_1, ssbm * rhs ) ;
	EXPECT_EQ( expected_1, ssbm.transpose() * rhs ) ;

	res.resize( ssbm.rows() ) ;
	res.setZero() ;
	EXPECT_EQ( res.rows(), 9 ) ;

	for( unsigned k = 0 ; k < 3 ; ++ k )
	{
		Eigen::VectorXd::SegmentReturnType seg ( res.segment ( 3*k, 3 ) ) ;
		ssbm.splitRowMultiply( k, rhs, seg ) ;
	}
	EXPECT_EQ( expected_2, res ) ;

	bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::BlockMatrixFlags::SYMMETRIC | bogus::BlockMatrixFlags::COMPRESSED | bogus::BlockMatrixFlags::COL_MAJOR > ssbm_col_major = ssbm ;
	//std::cout << ssbm_copy << std::endl ;

	res.setZero() ;

	for( unsigned k = 0 ; k < 3 ; ++ k )
	{
		Eigen::VectorXd::SegmentReturnType seg ( res.segment ( 3*k, 3 ) ) ;
		ssbm_col_major.splitRowMultiply( k, rhs, seg ) ;
	}
	EXPECT_EQ( expected_2, res ) ;


	bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::BlockMatrixFlags::SYMMETRIC > ussbm ( ssbm );

	rhs.resize( ussbm.rows() ) ;
	rhs.setOnes() ;

	EXPECT_EQ( true, ussbm.computeMinorIndex() ) ;

	for( unsigned k = 0 ; k < 3 ; ++ k )
	{
		Eigen::VectorXd::SegmentReturnType seg ( res.segment ( 3*k, 3 ) ) ;
		ssbm.splitRowMultiply( k, rhs, seg ) ;
	}
	EXPECT_EQ( 2*expected_2, res ) ;

	ussbm.cacheTranspose() ;

	res.setZero() ;
	for( unsigned k = 0 ; k < 3 ; ++ k )
	{
		Eigen::VectorXd::SegmentReturnType seg ( res.segment ( 3*k, 3 ) ) ;
		ssbm.splitRowMultiply( k, rhs, seg ) ;
	}
	EXPECT_EQ( expected_2, res ) ;

}


TEST( SparseBlock, ColMajor )
{
	Eigen::VectorXd expected_1(12) ;
	expected_1 << 13, 14, 15, 0, 0, 0, 0, 0, 0, 49, 54, 59 ;

	typedef Eigen::MatrixXd BlockT ;
	bogus::SparseBlockMatrix< BlockT, bogus::BlockMatrixFlags::COL_MAJOR > sbm ;
	sbm.setRows( 4, 3 ) ;
	sbm.setCols( 2, 4 ) ;

	BlockT baseBlock( 3, 4 ) ;
	baseBlock << 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1 ;

	sbm.insertBack( 3, 0 ) = 2 * baseBlock ;
	sbm.insertBack( 0, 1 ) = baseBlock ;
	sbm.insertBack( 3, 1 ) = 3 * baseBlock ;

	Eigen::VectorXd rhs( sbm.cols() ) ;
	EXPECT_EQ( 8, rhs.rows() ) ;
	rhs << 1, 2, 3, 4, 5, 6, 7, 8 ;

	EXPECT_EQ( expected_1, sbm * rhs ) ;
	std::cout << ( sbm * rhs ).transpose() << std::endl ;
	bogus::SparseBlockMatrix< BlockT, bogus::BlockMatrixFlags::COMPRESSED > rsbm ( sbm );
	EXPECT_EQ( 4u, rsbm.rowsOfBlocks() ) ;
	EXPECT_EQ( 4u, rsbm.majorIndex().outerSize() ) ;
	EXPECT_EQ( true, rsbm.majorIndex().valid ) ;
	EXPECT_EQ( 2u, rsbm.colsOfBlocks() ) ;
	std::cout << rsbm << std::endl ;
	std::cout << ( rsbm*rhs ) .transpose() << std::endl ;
	EXPECT_EQ( expected_1, rsbm * rhs ) ;
	bogus::SparseBlockMatrix< BlockT, bogus::BlockMatrixFlags::COL_MAJOR > tsbm ( rsbm.transpose() );
	EXPECT_EQ( 8u, tsbm.rows() ) ;
	EXPECT_EQ( 12u, tsbm.cols() ) ;
	EXPECT_EQ( 2u, tsbm.rowsOfBlocks() ) ;
	EXPECT_EQ( 4u, tsbm.colsOfBlocks() ) ;
	EXPECT_EQ( expected_1, tsbm.transpose() * rhs ) ;

	bogus::SparseBlockMatrix< BlockT, bogus::BlockMatrixFlags::COL_MAJOR > ttsbm ( tsbm.transpose() );
	EXPECT_EQ( 12u, ttsbm.rows() ) ;
	EXPECT_EQ( 8u, ttsbm.cols() ) ;
	EXPECT_EQ( 4u, ttsbm.rowsOfBlocks() ) ;
	EXPECT_EQ( 2u, ttsbm.colsOfBlocks() ) ;
	EXPECT_EQ( expected_1, ttsbm * rhs ) ;

}



