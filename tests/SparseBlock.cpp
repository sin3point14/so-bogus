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

	bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::BlockMatrixFlags::COMPRESSED > csbm( sbm ) ;
	EXPECT_EQ( expected_2, ( res.transpose() * csbm ).transpose() ) ;
	EXPECT_EQ( expected_2, ( csbm.transpose() * res )  ) ;
	EXPECT_EQ( expected_3, ( rhs.transpose() * csbm.transpose() ).transpose()  ) ;

	sbm.cacheTranspose();
	EXPECT_EQ( expected_2, ( res.transpose() * sbm ).transpose() ) ;
	EXPECT_EQ( expected_2, ( sbm.transpose() * res )  ) ;

	bogus::SparseBlockMatrix< Eigen::MatrixXd > structofsbm ;
	structofsbm.cloneStructure( sbm ) ;
	EXPECT_EQ( true, structofsbm.majorIndex().valid ) ;
	EXPECT_EQ( false, structofsbm.transposeCached() ) ;
	EXPECT_EQ( false, structofsbm.minorIndex().valid ) ;
	structofsbm = sbm ;

	bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::BlockMatrixFlags::COMPRESSED > copyofsbm ;
	copyofsbm = sbm ;
	EXPECT_EQ( true, copyofsbm.majorIndex().valid ) ;
	EXPECT_EQ( true, copyofsbm.transposeCached() ) ;
	EXPECT_EQ( true, copyofsbm.minorIndex().valid ) ;

	rhs.setOnes() ;
	EXPECT_EQ( expected_1, copyofsbm * rhs ) ;
	EXPECT_EQ( expected_2, ( copyofsbm.transpose() * res )  ) ;

	sbm = copyofsbm.transpose() ;
	EXPECT_EQ( true, sbm.majorIndex().valid ) ;
	EXPECT_EQ( expected_1, sbm.transpose() * rhs ) ;
	EXPECT_EQ( expected_2, ( sbm * res )  ) ;

	sbm = structofsbm.transpose() ;
	EXPECT_EQ( true, sbm.majorIndex().valid ) ;
	EXPECT_EQ( expected_1, sbm.transpose() * rhs ) ;
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
	Eigen::RowVectorXd expected_2(8) ;
	expected_2 << 6, 4, 2, 12, 21, 17, 13, 51 ;

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
	Eigen::RowVectorXd lhs( sbm.rows() ) ;
	EXPECT_EQ( 8, rhs.rows() ) ;
	EXPECT_EQ( 12, lhs.cols() ) ;
	rhs << 1, 2, 3, 4, 5, 6, 7, 8 ;
	lhs << 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 ;

	EXPECT_EQ( expected_1, sbm * rhs ) ;
//	std::cout << sbm << std::endl ;
	bogus::SparseBlockMatrix< BlockT, bogus::BlockMatrixFlags::COMPRESSED > rsbm ( sbm );
	EXPECT_EQ( 4u, rsbm.rowsOfBlocks() ) ;
	EXPECT_EQ( 4u, rsbm.majorIndex().outerSize() ) ;
	EXPECT_EQ( true, rsbm.majorIndex().valid ) ;
	EXPECT_EQ( 2u, rsbm.colsOfBlocks() ) ;
	EXPECT_EQ( expected_1, rsbm * rhs ) ;
	bogus::SparseBlockMatrix< BlockT > ursbm ( sbm );
	EXPECT_EQ( 4u, ursbm.rowsOfBlocks() ) ;
	EXPECT_EQ( 4u, ursbm.majorIndex().outerSize() ) ;
	EXPECT_EQ( true, ursbm.majorIndex().valid ) ;
	EXPECT_EQ( 2u, ursbm.colsOfBlocks() ) ;
	EXPECT_EQ( expected_1, ursbm * rhs ) ;
	bogus::SparseBlockMatrix< BlockT, bogus::BlockMatrixFlags::COL_MAJOR > tsbm ( rsbm.transpose() );
	EXPECT_EQ( 8u, tsbm.rows() ) ;
	EXPECT_EQ( 12u, tsbm.cols() ) ;
	EXPECT_EQ( 2u, tsbm.rowsOfBlocks() ) ;
	EXPECT_EQ( 4u, tsbm.colsOfBlocks() ) ;
	EXPECT_EQ( expected_1, tsbm.transpose() * rhs ) ;
	bogus::SparseBlockMatrix< BlockT, bogus::BlockMatrixFlags::COL_MAJOR | bogus::BlockMatrixFlags::COMPRESSED > ttsbm ( tsbm.transpose() );
	EXPECT_EQ( 12u, ttsbm.rows() ) ;
	EXPECT_EQ( 8u, ttsbm.cols() ) ;
	EXPECT_EQ( 4u, ttsbm.rowsOfBlocks() ) ;
	EXPECT_EQ( 2u, ttsbm.colsOfBlocks() ) ;
	EXPECT_EQ( expected_1, ttsbm * rhs ) ;
	bogus::SparseBlockMatrix< BlockT, bogus::BlockMatrixFlags::COL_MAJOR > uttsbm ( tsbm.transpose() );
	EXPECT_EQ( 12u, uttsbm.rows() ) ;
	EXPECT_EQ( 8u, uttsbm.cols() ) ;
	EXPECT_EQ( 4u, uttsbm.rowsOfBlocks() ) ;
	EXPECT_EQ( 2u, uttsbm.colsOfBlocks() ) ;
	EXPECT_EQ( expected_1, uttsbm * rhs ) ;

	EXPECT_EQ( expected_2, lhs * uttsbm ) ;
	EXPECT_EQ( expected_2, lhs * ttsbm ) ;
	EXPECT_EQ( expected_2, lhs * tsbm.transpose() ) ;
	EXPECT_EQ( expected_2, lhs * ursbm ) ;
	EXPECT_EQ( expected_2, lhs * rsbm ) ;
	EXPECT_EQ( expected_2, lhs * sbm ) ;

	uttsbm.cacheTranspose();
	ttsbm.cacheTranspose();
	tsbm.cacheTranspose();
	ursbm.cacheTranspose();
	rsbm.cacheTranspose();
	sbm.cacheTranspose();

	EXPECT_EQ( expected_2, lhs * uttsbm ) ;
	EXPECT_EQ( expected_2, lhs * ttsbm ) ;
	EXPECT_EQ( expected_2, lhs * tsbm.transpose() ) ;
	EXPECT_EQ( expected_2, lhs * ursbm ) ;
	EXPECT_EQ( expected_2, lhs * rsbm ) ;
	EXPECT_EQ( expected_2, lhs * sbm ) ;

}



