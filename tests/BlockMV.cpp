/*
 * Any copyright is dedicated to the Public Domain.
 * http://creativecommons.org/publicdomain/zero/1.0/
*/

#include "Core/Block.impl.hpp"

#include <Eigen/Core>

#include <gtest/gtest.h>

#include <iostream>

class SmallSBM : public ::testing::Test {
protected:

	virtual void SetUp() {

		//Eigen::initParallel()
		Eigen::MatrixXf A = Eigen::MatrixXf::Zero(1,1); A = A*A;

		EXPECT_TRUE( bogus::IsTransposable< BlockT >::Value ) ;
		sbm.setRows( 5, 3 ) ;
		sbm.setCols( 2, 4 ) ;

		sbm.insertBack( 0, 1 ) = BlockT::Ones( 3,4 ) ;
		sbm.insertBack( 3, 0 ) = 2 * BlockT::Ones( 3,4 ) ;
		sbm.insertBack( 3, 1 ) = 3 * BlockT::Ones( 3,4 ) ;

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
	typedef bogus::SparseBlockMatrix< BlockT, bogus::flags::UNCOMPRESSED > SBMT ;
	SBMT sbm ;
	Eigen::VectorXd rhs  ;

	Eigen::VectorXd expected_1 ; // sbm * rhs
	Eigen::VectorXd expected_2 ; // sbm.transpose() * (sbm * rhs)
	Eigen::VectorXd expected_3 ; // sbm * ( sbm.transpose() * (sbm * rhs) )

} ;

TEST_F( SmallSBM, MatrixVector )
{

	Eigen::VectorXd res ( sbm.rows() ) ;

	res.setZero() ;
	sbm.multiply< false >( rhs,res ) ;

	EXPECT_EQ( expected_1, res ) ;

	res.noalias() += sbm*rhs ;
	res -= sbm*rhs ;
	EXPECT_EQ( expected_1, sbm*rhs ) ;
	rhs.setZero() ;
	sbm.multiply< true >( res, rhs ) ;

	EXPECT_EQ( expected_2, rhs ) ;
	EXPECT_EQ( expected_2, ( res.transpose() * sbm ).transpose() ) ;
	EXPECT_EQ( expected_2, ( sbm.transpose() * res )  ) ;
	EXPECT_EQ( expected_3, ( rhs.transpose() * sbm.transpose() ).transpose()  ) ;

	bogus::SparseBlockMatrix< Eigen::MatrixXd > csbm( sbm ) ;
	EXPECT_EQ( expected_2, ( res.transpose() * csbm ).transpose() ) ;
	EXPECT_EQ( expected_2, ( csbm.transpose() * res )  ) ;
	EXPECT_EQ( expected_3, ( rhs.transpose() * csbm.transpose() ).transpose()  ) ;

	EXPECT_FALSE( sbm.minorIndex().valid );
	EXPECT_FALSE( sbm.transposeIndex().valid );
	sbm.cacheTranspose();
	ASSERT_TRUE( sbm.minorIndex().valid );
	ASSERT_TRUE( sbm.transposeIndex().valid );
	EXPECT_EQ( expected_2, ( res.transpose() * sbm ).transpose() ) ;
	EXPECT_EQ( expected_2, ( sbm.transpose() * res )  ) ;

	bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::flags::UNCOMPRESSED  > structofsbm ;
	structofsbm.cloneStructure( sbm ) ;
	EXPECT_TRUE( structofsbm.majorIndex().valid ) ;
	EXPECT_FALSE( structofsbm.transposeCached() ) ;
	EXPECT_TRUE( structofsbm.minorIndex().valid ) ;
	structofsbm = sbm ;

	bogus::SparseBlockMatrix< BlockT, bogus::flags::UNCOMPRESSED > closbm ;
	closbm.setRows( 5, 2 ) ;
	closbm.setCols( 2, 3 ) ;
	closbm.cloneIndex( sbm.majorIndex() );
	EXPECT_EQ( closbm.nBlocks(), sbm.nBlocks() ) ;
	EXPECT_EQ( closbm.rowsOfBlocks(), sbm.rowsOfBlocks() ) ;
	EXPECT_EQ( closbm.colsOfBlocks(), sbm.colsOfBlocks() ) ;
	EXPECT_EQ( closbm.rows(), 10 ) ;
	EXPECT_EQ( closbm.rowOffset( 3 ), 6 ) ;

	bogus::SparseBlockMatrix< Eigen::MatrixXd > copyofsbm ;
	copyofsbm = sbm ;
	EXPECT_TRUE( copyofsbm.majorIndex().valid ) ;
	EXPECT_TRUE( copyofsbm.transposeCached() ) ;
	EXPECT_EQ( copyofsbm.transposeBlocks().size(), sbm.nBlocks() ) ;
	EXPECT_EQ( sbm.transposeBlocks().size(), copyofsbm.nBlocks() ) ;
	EXPECT_TRUE( copyofsbm.minorIndex().valid ) ;

	EXPECT_EQ( 2u, copyofsbm.blockPtr( 3, 1 ) ) ;
	EXPECT_EQ( copyofsbm.InvalidBlockPtr, copyofsbm.blockPtr( 0, 0 ) ) ;
	EXPECT_EQ( copyofsbm.InvalidBlockPtr, copyofsbm.blockPtr( 1, 1 ) ) ;

	rhs.setOnes() ;
	EXPECT_EQ( expected_1, copyofsbm * rhs ) ;
	EXPECT_EQ( expected_2, ( copyofsbm.transpose() * res )  ) ;

	sbm = copyofsbm.transpose() ;
	EXPECT_TRUE( sbm.majorIndex().valid ) ;
	EXPECT_FALSE( sbm.minorIndex().valid ) ;
	EXPECT_TRUE( sbm.transposeIndex().valid ) ;
	EXPECT_EQ( expected_1, sbm.transpose() * rhs ) ;
	EXPECT_EQ( expected_2, ( sbm * res )  ) ;

	bogus::SparseBlockMatrix< Eigen::MatrixXd > prod( sbm * copyofsbm );

	sbm = structofsbm.transpose() ;
	EXPECT_TRUE( sbm.majorIndex().valid ) ;
	EXPECT_EQ( expected_1, sbm.transpose() * rhs ) ;
	EXPECT_EQ( expected_2, ( sbm * res )  ) ;


	Eigen::MatrixXd mrhs( res.rows(), 2 ) ;
	mrhs.col(0) = res ;
	mrhs.col(1) = 2 * res ;
	Eigen::MatrixXd mres = sbm * mrhs ;
	EXPECT_EQ( expected_2, mres.col(0) ) ;
	EXPECT_EQ( 2*expected_2, mres.col(1) ) ;

	mres = ( mrhs.transpose() * sbm.transpose() ).transpose() ;
	EXPECT_EQ( expected_2, mres.col(0) ) ;
	EXPECT_EQ( 2*expected_2, mres.col(1) ) ;

	sbm.setBlocksToZero() ;
	mres = sbm * mrhs ;
	EXPECT_TRUE( mres.isZero() ) ;
}

TEST_F( SmallSBM, ExprVector )
{
	bogus::SparseBlockMatrix< BlockT > sbm_2 = sbm/2 ;

	Eigen::VectorXd res ( sbm.rows() ) ;
	res.setZero() ;

	(2*sbm_2).multiply< false >( rhs,res ) ;
	EXPECT_EQ( expected_1, res ) ;

	res.setZero() ;

	(3*sbm + 10*sbm_2).multiply< false >( rhs, res ) ;
	res /= 8;

	EXPECT_EQ( expected_1, res ) ;

	res.noalias() += sbm*rhs ;
	res -= sbm*rhs ;
	EXPECT_EQ( expected_1, sbm*rhs ) ;
	rhs.setZero() ;

	sbm.transpose().multiply< false >( res, rhs ) ;
	EXPECT_EQ( expected_2, rhs ) ;

	bogus::SparseBlockMatrix< BlockT > sbm3 = sbm + sbm_2 ;

	Eigen::VectorXd mmv = sbm.transpose() * ( sbm3 * rhs ) ;
	Eigen::VectorXd res2( mmv.rows() ) ;

	(sbm.transpose() * sbm3).multiply< false >( rhs, res2 ) ;
	EXPECT_EQ( res2, mmv ) ;

	mmv += 6*sbm3.transpose() * ( sbm * rhs ) ;

	( (3*sbm).transpose() * (2*sbm3)).multiply< true >( rhs, res2, 1, 1 ) ;
	EXPECT_EQ( res2, mmv ) ;

}

TEST( BlockMV, Symmetric )
{
	Eigen::VectorXd expected_1(9), expected_2(9) ;
	expected_1 << 9, 7, 5, 2, 4, 6, 9, 9, 9 ;
	expected_2 << 6, 4, 2, 2, 4, 6, 0, 0, 0 ;

	Eigen::VectorXd res, rhs ;

	bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::SYMMETRIC > ssbm ;
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

	{
		Eigen::MatrixXd mrhs( rhs.rows(), 2 ) ;
		mrhs.col(0) = rhs ;
		mrhs.col(1) = 2 * rhs ;
		Eigen::MatrixXd mres = ssbm * mrhs ;
		EXPECT_EQ( expected_1, mres.col(0) ) ;
		EXPECT_EQ( 2*expected_1, mres.col(1) ) ;

		ssbm.cacheTranspose() ;

		mres = ( mrhs.transpose() * ssbm.transpose() ).transpose() ;
		EXPECT_EQ( expected_1, mres.col(0) ) ;
		EXPECT_EQ( 2*expected_1, mres.col(1) ) ;
	}


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

	bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::SYMMETRIC | bogus::flags::COL_MAJOR > ssbm_col_major = ssbm ;
	//std::cout << ssbm_col_major << std::endl ;

	res.setZero() ;

	for( unsigned k = 0 ; k < 3 ; ++ k )
	{
		Eigen::VectorXd::SegmentReturnType seg ( res.segment ( 3*k, 3 ) ) ;
		ssbm_col_major.splitRowMultiply( k, rhs, seg ) ;
	}
	EXPECT_EQ( expected_2, res ) ;


	bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::SYMMETRIC > ussbm ( ssbm );

	rhs.resize( ussbm.rows() ) ;
	rhs.setOnes() ;

	EXPECT_TRUE( ussbm.computeMinorIndex() ) ;

	for( unsigned k = 0 ; k < 3 ; ++ k )
	{
		Eigen::VectorXd::SegmentReturnType seg ( res.segment ( 3*k, 3 ) ) ;
		ussbm.splitRowMultiply( k, rhs, seg ) ;
	}
	EXPECT_EQ( 2*expected_2, res ) ;

	ussbm.cacheTranspose() ;

	res.setZero() ;
	for( unsigned k = 0 ; k < 3 ; ++ k )
	{
		Eigen::VectorXd::SegmentReturnType seg ( res.segment ( 3*k, 3 ) ) ;
		ussbm.splitRowMultiply( k, rhs, seg ) ;
	}
	EXPECT_EQ( expected_2, res ) ;


	bogus::SparseBlockMatrix< Eigen::Matrix3d > nsbm ( ssbm );
	//std::cout << nsbm << std::endl ;
	res.setZero() ;
	for( unsigned k = 0 ; k < 3 ; ++ k )
	{
		Eigen::VectorXd::SegmentReturnType seg ( res.segment ( 3*k, 3 ) ) ;
		nsbm.splitRowMultiply( k, rhs, seg ) ;
	}
	EXPECT_EQ( expected_2, res ) ;

	bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::COL_MAJOR | bogus::flags::SYMMETRIC > resbm
			= nsbm.transpose() ;
	res.setZero() ;
	for( unsigned k = 0 ; k < 3 ; ++ k )
	{
		Eigen::VectorXd::SegmentReturnType seg ( res.segment ( 3*k, 3 ) ) ;
		resbm.splitRowMultiply( k, rhs, seg ) ;
	}
	EXPECT_EQ( expected_2, res ) ;
	//std::cout << resbm << std::endl ;
}


TEST( BlockMV, ColMajor )
{
	Eigen::VectorXd expected_1(12) ;
	expected_1 << 13, 14, 15, 0, 0, 0, 0, 0, 0, 49, 54, 59 ;
	Eigen::RowVectorXd expected_2(8) ;
	expected_2 << 6, 4, 2, 12, 21, 17, 13, 51 ;

	typedef Eigen::MatrixXd BlockT ;
	bogus::SparseBlockMatrix< BlockT, bogus::flags::COL_MAJOR > sbm ;
	sbm.setRows( 4, 3 ) ;
	sbm.setCols( 2, 4 ) ;

	BlockT baseBlock( 3, 4 ) ;
	baseBlock << 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1 ;

	sbm.insertBack( 3, 0 ) = 2 * baseBlock ;
	sbm.insertBack( 0, 1 ) = baseBlock ;
	sbm.insertBack( 3, 1 ) = 3 * baseBlock ;
	sbm.finalize() ;

	Eigen::VectorXd rhs( sbm.cols() ) ;
	Eigen::RowVectorXd lhs( sbm.rows() ) ;
	EXPECT_EQ( 8, rhs.rows() ) ;
	EXPECT_EQ( 12, lhs.cols() ) ;
	rhs << 1, 2, 3, 4, 5, 6, 7, 8 ;
	lhs << 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 ;

	EXPECT_EQ( expected_1, sbm * rhs ) ;
	//	std::cout << sbm << std::endl ;
	bogus::SparseBlockMatrix< BlockT, bogus::flags::UNCOMPRESSED > rsbm ( sbm );
	EXPECT_EQ( 4, rsbm.rowsOfBlocks() ) ;
	EXPECT_EQ( 4, rsbm.majorIndex().outerSize() ) ;
	EXPECT_TRUE( rsbm.majorIndex().valid ) ;
	EXPECT_EQ( 2, rsbm.colsOfBlocks() ) ;
	EXPECT_EQ( expected_1, rsbm * rhs ) ;
	bogus::SparseBlockMatrix< BlockT > ursbm ( sbm );
	EXPECT_EQ( 4, ursbm.rowsOfBlocks() ) ;
	EXPECT_EQ( 4, ursbm.majorIndex().outerSize() ) ;
	EXPECT_TRUE( ursbm.majorIndex().valid ) ;
	EXPECT_EQ( 2, ursbm.colsOfBlocks() ) ;
	EXPECT_EQ( expected_1, ursbm * rhs ) ;
	bogus::SparseBlockMatrix< BlockT, bogus::flags::COL_MAJOR > tsbm ( rsbm.transpose() );
	EXPECT_EQ( 8, tsbm.rows() ) ;
	EXPECT_EQ( 12, tsbm.cols() ) ;
	EXPECT_EQ( 2, tsbm.rowsOfBlocks() ) ;
	EXPECT_EQ( 4, tsbm.colsOfBlocks() ) ;
	EXPECT_EQ( expected_1, tsbm.transpose() * rhs ) ;
	bogus::SparseBlockMatrix< BlockT, bogus::flags::COL_MAJOR | bogus::flags::UNCOMPRESSED > ttsbm ( tsbm.transpose() );
	EXPECT_EQ( 12, ttsbm.rows() ) ;
	EXPECT_EQ( 8, ttsbm.cols() ) ;
	EXPECT_EQ( 4, ttsbm.rowsOfBlocks() ) ;
	EXPECT_EQ( 2, ttsbm.colsOfBlocks() ) ;
	EXPECT_EQ( expected_1, ttsbm * rhs ) ;
	bogus::SparseBlockMatrix< BlockT, bogus::flags::COL_MAJOR > uttsbm ( tsbm.transpose() );
	EXPECT_EQ( 12, uttsbm.rows() ) ;
	EXPECT_EQ( 8, uttsbm.cols() ) ;
	EXPECT_EQ( 4, uttsbm.rowsOfBlocks() ) ;
	EXPECT_EQ( 2, uttsbm.colsOfBlocks() ) ;
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


TEST_F( SmallSBM, RowMultiply )
{
	Eigen::VectorXd res ( sbm.rows() ), res2( sbm.rows() ) ;
	res.setZero() ; res2.setZero() ;

	typedef bogus::Segmenter< bogus::internal::DYNAMIC, Eigen::VectorXd, SBMT::Index >
			SegmenterType ;
	SegmenterType resSeg( res, sbm.rowOffsets() ), rhsSeg( rhs, sbm.colOffsets() ) ;

	typedef SegmenterType::ReturnType Seg ;
	for( int i = 0 ; i < sbm.rowsOfBlocks() ; ++i ) {
		Seg seg = resSeg[i] ;
		sbm.rowMultiply< false >( i, rhs, seg ) ;
	}

	EXPECT_EQ( expected_1, res ) ;
	EXPECT_EQ( sbm*rhs, res ) ;

	for( int i = 0 ; i < sbm.colsOfBlocks() ; ++i ) {
		sbm.colMultiply< false >( i, rhsSeg[i], res2 ) ;
	}
	EXPECT_EQ( expected_1, res2 ) ;

	Eigen::VectorXd u = Eigen::VectorXd::Zero( rhs.rows() ) ;

	for( int i = 0 ; i < sbm.rowsOfBlocks() ; ++i ) {
		sbm.colMultiply< true >( i, resSeg[i], u ) ;
	}

	EXPECT_EQ( u, expected_2 ) ;
	EXPECT_EQ( u, sbm.transpose() * res ) ;

	rhs.setZero() ;
	for( int i = 0 ; i < sbm.colsOfBlocks() ; ++i ) {
		Seg seg = rhsSeg[i] ;
		sbm.rowMultiply< true >( i, res, seg ) ;
	}
	EXPECT_EQ( rhs, expected_2 ) ;
}

TEST_F( SmallSBM, Compound )
{
	typedef Eigen::Matrix< double, 3, 3> OBlockT ;
	typedef bogus::SparseBlockMatrix< OBlockT > OSBMT ;

	OSBMT osbm ;
	osbm.setRows( 5 );
	osbm.setCols( 2 );

	osbm.insertBack(0,0).setIdentity() ;
	osbm.insertBack(2,1).setOnes() ;
	osbm.insertBack(3,1).setConstant(2) ;
	osbm.finalize();

	bogus::CompoundBlockMatrix< true, SBMT, OSBMT > comp ( sbm, osbm) ;
	ASSERT_EQ( comp.rows(), sbm.rows() ) ;
	ASSERT_EQ( comp.rows(), osbm.rows() ) ;
	ASSERT_EQ( comp.cols(), sbm.cols() + osbm.cols() ) ;
	ASSERT_EQ( comp.rowsOfBlocks(), sbm.rowsOfBlocks() ) ;
	ASSERT_EQ( comp.rowsOfBlocks(), osbm.rowsOfBlocks() ) ;
	ASSERT_EQ( comp.colsOfBlocks(), sbm.colsOfBlocks() + osbm.colsOfBlocks() ) ;

	rhs.setOnes( comp.cols() ) ;
	Eigen::VectorXd res( sbm.rows() )	 ;

	res = comp * rhs ;
	Eigen::VectorXd exp = sbm*rhs.head( sbm.cols() )  + osbm*rhs.tail( osbm.cols() ) ;
	ASSERT_EQ( res, exp ) ;

	rhs.head( sbm.cols() ) = sbm.transpose() * res ;
	rhs.tail( osbm.cols() ) = osbm.transpose() * res ;
	ASSERT_EQ( rhs, comp.transpose() * res ) ;

	// row, col-multiply

	res.setZero() ;
	Eigen::VectorXd res2 = res ;

	typedef bogus::Segmenter< bogus::internal::DYNAMIC, Eigen::VectorXd, SBMT::Index >
			SegmenterType ;
	SegmenterType resSeg( res, comp.rowOffsets() ), rhsSeg( rhs, comp.colOffsets() ) ;

	typedef SegmenterType::ReturnType Seg ;
	for( int i = 0 ; i < comp.rowsOfBlocks() ; ++i ) {
		Seg seg = resSeg[i] ;
		comp.rowMultiply< false >( i, rhs, seg ) ;
	}

	EXPECT_EQ( comp*rhs, res ) ;

	for( int i = 0 ; i < comp.colsOfBlocks() ; ++i ) {
		comp.colMultiply< false >( i, rhsSeg[i], res2 ) ;
	}
	EXPECT_EQ( res, res2 ) ;

	Eigen::VectorXd u = Eigen::VectorXd::Zero( rhs.rows() ) ;

	for( int i = 0 ; i < comp.rowsOfBlocks() ; ++i ) {
		comp.colMultiply< true >( i, resSeg[i], u ) ;
	}
	EXPECT_EQ( u, comp.transpose() * res ) ;

	rhs.setZero() ;
	for( int i = 0 ; i < comp.colsOfBlocks() ; ++i ) {
		Seg seg = rhsSeg[i] ;
		comp.rowMultiply< true >( i, res, seg ) ;
	}
	EXPECT_EQ( rhs, u ) ;


}
