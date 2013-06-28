#include <iostream>

//#define BOGUS_DONT_PARALLELIZE

#include "Core/Block.impl.hpp"
#include "Core/Block.io.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>

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

	EXPECT_EQ( 2u, sbm.blockPtr( 3, 1 ) ) ;
	EXPECT_EQ( sbm.InvalidBlockPtr, sbm.blockPtr( 0, 0 ) ) ;
	EXPECT_EQ( sbm.InvalidBlockPtr, sbm.blockPtr( 1, 1 ) ) ;

	Eigen::VectorXd rhs ( sbm.cols() ) ;
	Eigen::VectorXd res ( sbm.rows() ) ;

	rhs.setOnes() ;
	res.setZero() ;
	sbm.multiply< false >( rhs, res ) ;

	EXPECT_EQ( expected_1, res ) ;
	EXPECT_EQ( expected_1, sbm*rhs ) ;

	rhs.setZero() ;
	sbm.multiply< true >( res, rhs ) ;

	EXPECT_EQ( expected_2, rhs ) ;
	EXPECT_EQ( expected_2, ( res.transpose() * sbm ).transpose() ) ;
	EXPECT_EQ( expected_2, ( sbm.transpose() * res )  ) ;
	EXPECT_EQ( expected_3, ( rhs.transpose() * sbm.transpose() ).transpose()  ) ;

	bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::flags::COMPRESSED > csbm( sbm ) ;
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

	bogus::SparseBlockMatrix< Eigen::MatrixXd > structofsbm ;
	structofsbm.cloneStructure( sbm ) ;
	EXPECT_TRUE( structofsbm.majorIndex().valid ) ;
	EXPECT_FALSE( structofsbm.transposeCached() ) ;
	EXPECT_TRUE( structofsbm.minorIndex().valid ) ;
	structofsbm = sbm ;

	bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::flags::COMPRESSED > copyofsbm ;
	copyofsbm = sbm ;
	EXPECT_TRUE( copyofsbm.majorIndex().valid ) ;
	EXPECT_TRUE( copyofsbm.transposeCached() ) ;
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
}

TEST( SparseBlock, Symmetric )
{
	Eigen::VectorXd expected_1(9), expected_2(9) ;
	expected_1 << 9, 7, 5, 2, 4, 6, 9, 9, 9 ;
	expected_2 << 6, 4, 2, 2, 4, 6, 0, 0, 0 ;

	Eigen::VectorXd res, rhs ;

	bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::SYMMETRIC | bogus::flags::COMPRESSED > ssbm ;
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

	bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::SYMMETRIC | bogus::flags::COMPRESSED | bogus::flags::COL_MAJOR > ssbm_col_major = ssbm ;
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


TEST( SparseBlock, ColMajor )
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

	Eigen::VectorXd rhs( sbm.cols() ) ;
	Eigen::RowVectorXd lhs( sbm.rows() ) ;
	EXPECT_EQ( 8, rhs.rows() ) ;
	EXPECT_EQ( 12, lhs.cols() ) ;
	rhs << 1, 2, 3, 4, 5, 6, 7, 8 ;
	lhs << 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 ;

	EXPECT_EQ( expected_1, sbm * rhs ) ;
	//	std::cout << sbm << std::endl ;
	bogus::SparseBlockMatrix< BlockT, bogus::flags::COMPRESSED > rsbm ( sbm );
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
	bogus::SparseBlockMatrix< BlockT, bogus::flags::COL_MAJOR | bogus::flags::COMPRESSED > ttsbm ( tsbm.transpose() );
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

TEST( SparseBlock, MMult )
{
	typedef Eigen::Matrix< double, 3, 4 > BlockT ;
	BlockT sample ;
	sample<< 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, 3, 4 ;

	Eigen::VectorXd rhs( 12 ), expected_1( 8 ) ;
	rhs << 4, 2, 1, 0, 7, 6, 9, 8, 2, 4, 6, 8 ;

	expected_1 << 30410, 31550, 32690, 33830, 46674, 46674, 46674, 46674 ;

	{
		bogus::SparseBlockMatrix< BlockT, bogus::flags::COMPRESSED > sbm ;
		sbm.setRows( 4, 3 ) ;
		sbm.setCols( 2, 4 ) ;

		sbm.insertBack( 0, 1 ) = BlockT::Ones() ;
		sbm.insertBack( 1, 0 ) = sample ;
		sbm.insertBack( 2, 0 ) = 3 * BlockT::Ones() ;
		sbm.insertBack( 2, 1 ) = 5 * BlockT::Ones() ;

		sbm.finalize() ;
		ASSERT_FALSE( sbm.minorIndex().valid ) ;

		bogus::SparseBlockMatrix< Eigen::Matrix3d > mm = sbm*sbm.transpose() ;
		EXPECT_EQ( mm.block(0,2), mm.block(2,0).transpose() ) ;
		EXPECT_EQ( mm.block(2,1), mm.block(1,2).transpose() ) ;
		EXPECT_EQ( mm.diagonal(2), Eigen::Matrix3d::Constant( 136 ) ) ;

		mm.setFromProduct< false >( sbm*sbm.transpose() ) ;
		EXPECT_EQ( mm.block(0,2), mm.block(2,0).transpose() ) ;
		EXPECT_EQ( mm.block(2,1), mm.block(1,2).transpose() ) ;
		EXPECT_EQ( mm.diagonal(2), Eigen::Matrix3d::Constant( 136 ) ) ;

		bogus::SparseBlockMatrix< Eigen::Matrix4d, bogus::flags::COMPRESSED | bogus::flags::SYMMETRIC > mm_t = sbm.transpose()*sbm ;
		ASSERT_EQ( 3u, mm_t.nBlocks() ) ;
		EXPECT_EQ( mm_t.diagonal(1), Eigen::Matrix4d::Constant( 78 ) ) ;
		mm_t.setFromProduct< false > ( sbm.transpose()*sbm ) ;
		ASSERT_EQ( 3u, mm_t.nBlocks() ) ;
		EXPECT_EQ( mm_t.diagonal(1), Eigen::Matrix4d::Constant( 78 ) ) ;

		EXPECT_TRUE( mm_t.minorIndex().valid ) ;

		bogus::SparseBlockMatrix< Eigen::MatrixXd > mm_t_mt = mm_t  * sbm.transpose() ;
		EXPECT_EQ( expected_1, mm_t_mt*rhs ) ;
		EXPECT_EQ( expected_1, mm_t * ( sbm.transpose() * rhs )  ) ;
		bogus::SparseBlockMatrix< Eigen::MatrixXd > mm_t_mt_chk = sbm.transpose() * mm;
		EXPECT_EQ( expected_1, mm_t_mt_chk*rhs );
		mm_t_mt.setFromProduct< false > ( mm_t  * sbm.transpose() ) ;
		EXPECT_EQ( expected_1, mm_t_mt*rhs );
		bogus::SparseBlockMatrix< Eigen::MatrixXd > mm_t_mt_t = sbm * mm_t.transpose() ;
		EXPECT_EQ( expected_1, mm_t_mt_t.transpose() *rhs );

		bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::flags::SYMMETRIC > mm_t2( mm_t * mm_t ) ;
		EXPECT_EQ( mm_t2*expected_1, mm_t * ( mm_t * expected_1 ) );


		bogus::SparseBlockMatrix< Eigen::Matrix< double, 4, 3 > > mm_t_mt2 = 2 * mm_t_mt - sbm.transpose() * sbm * sbm.transpose() ;
		EXPECT_EQ( expected_1, mm_t_mt2*rhs ) ;
	}
	{
		bogus::SparseBlockMatrix< BlockT, bogus::flags::COL_MAJOR | bogus::flags::COMPRESSED > sbm ;
		sbm.setRows( 4, 3 ) ;
		sbm.setCols( 2, 4 ) ;

		sbm.insertBack( 1, 0 ) = sample ;
		sbm.insertBack( 2, 0 ) = 3 * BlockT::Ones() ;
		sbm.insertBack( 0, 1 ) = BlockT::Ones() ;
		sbm.insertBack( 2, 1 ) = 5 * BlockT::Ones() ;

		sbm.finalize() ;
		ASSERT_FALSE( sbm.minorIndex().valid ) ;

		bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::COL_MAJOR > mm = sbm*sbm.transpose() ;

		EXPECT_EQ( mm.block(0,2), mm.block(2,0).transpose() ) ;
		EXPECT_EQ( mm.block(2,1), mm.block(1,2).transpose() ) ;
		EXPECT_EQ( mm.diagonal(2), Eigen::Matrix3d::Constant( 136 ) ) ;

		mm.setFromProduct< false >( sbm*sbm.transpose() ) ;
		EXPECT_EQ( mm.block(0,2), mm.block(2,0).transpose() ) ;
		EXPECT_EQ( mm.block(2,1), mm.block(1,2).transpose() ) ;
		EXPECT_EQ( mm.diagonal(2), Eigen::Matrix3d::Constant( 136 ) ) ;

		bogus::SparseBlockMatrix< Eigen::Matrix4d, bogus::flags::COL_MAJOR | bogus::flags::COMPRESSED | bogus::flags::SYMMETRIC > mm_t = sbm.transpose()*sbm ;
		ASSERT_EQ( mm_t.nBlocks(), 3u ) ;
		EXPECT_EQ( mm_t.diagonal(1), Eigen::Matrix4d::Constant( 78 ) ) ;

		mm_t.setFromProduct< false > ( sbm.transpose()*sbm ) ;
		ASSERT_EQ( mm_t.nBlocks(), 3u ) ;
		EXPECT_EQ( mm_t.diagonal(1), Eigen::Matrix4d::Constant( 78 ) ) ;

		EXPECT_TRUE( mm_t.minorIndex().valid ) ;

		bogus::SparseBlockMatrix< Eigen::MatrixXd > mm_t_mt = mm_t  * sbm.transpose() ;
		EXPECT_EQ( expected_1, mm_t_mt*rhs ) ;
		EXPECT_EQ( expected_1, mm_t * ( sbm.transpose() * rhs )  ) ;
		bogus::SparseBlockMatrix< Eigen::MatrixXd > mm_t_mt_chk = sbm.transpose() * mm;
		EXPECT_EQ( expected_1, mm_t_mt_chk*rhs );
		mm_t_mt.setFromProduct< false > ( mm_t  * sbm.transpose() ) ;
		EXPECT_EQ( expected_1, mm_t_mt*rhs );
		bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::flags::COL_MAJOR > mm_t_mt_t = sbm * mm_t.transpose() ;
		EXPECT_EQ( expected_1, mm_t_mt_t.transpose() *rhs );

		bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::flags::SYMMETRIC > mm_t2( mm_t * mm_t ) ;
		EXPECT_EQ( mm_t2*expected_1, mm_t * ( mm_t * expected_1 ) );

		bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::flags::SYMMETRIC > mm_t3( (2 * mm_t * mm_t).transpose() ) ;
		EXPECT_EQ( mm_t3*expected_1, mm_t * ( mm_t * expected_1 * 2 ) );
	}

}

TEST( SparseBlock, Inv )
{
	const Eigen::Vector2d expected_1 ( .5, .5 ) ;
	Eigen::VectorXd expected_2( 6 ) ;
	expected_2 << .25, .25, .25, .125, .125, .125 ;

	bogus::DenseLDLT< double, 2 > ldlt(  2 * Eigen::Matrix2d::Identity() ) ;
	EXPECT_EQ( expected_1, ldlt * Eigen::Vector2d::Ones() )  ;

	Eigen::Matrix2d matRes = ldlt * Eigen::Matrix2d::Ones() ;
	EXPECT_EQ( expected_1, matRes.diagonal() )  ;

	typedef bogus::LU< Eigen::MatrixBase< Eigen::MatrixXd > > BlockT ;
	bogus::SparseBlockMatrix< BlockT > isbm ;
	isbm.setRows( 2, 3 ) ;
	isbm.setCols( 2, 3 ) ;
	isbm.insertBack( 0, 0 ).compute( 4 * Eigen::Matrix3d::Identity() ) ;
	isbm.insertBack( 1, 1 ).compute( 8 * Eigen::Matrix3d::Identity() ) ;
	isbm.finalize() ;

	Eigen::VectorXd rhs( 6 ) ;
	rhs.setOnes() ;

	EXPECT_EQ( expected_2, ( isbm * rhs ) ) ;

}

TEST( SparseBlock, Sparse )
{
	Eigen::VectorXd expected_1( 6 ) ;
	expected_1 << 2, 4, 6, 16, 40, 6 ;
	Eigen::VectorXd expected_2( 6 ) ;
	expected_2 << 0.5, 1, 1.5, 1, 0.625, 6 ;

	typedef Eigen::SparseMatrix< double > BlockT ;
	bogus::SparseBlockMatrix< BlockT > sbm ;
	sbm.reserve( 10 ) ;
	sbm.setRows( 2, 3 ) ;
	sbm.setCols( 2, 3 ) ;

	BlockT& b1 = sbm.insertBackAndResize( 0, 0 ) ;
	b1.insert( 0, 0 ) = 2;
	b1.insert( 1, 1 ) = 2;
	b1.insert( 2, 2 ) = 2;

	BlockT& b2 = sbm.insertBackAndResize( 1, 1 ) ;
	b2.insert( 0, 0 ) = 4;
	b2.insert( 1, 1 ) = 8;
	b2.insert( 2, 2 ) = 1;

	sbm.finalize() ;
	Eigen::VectorXd rhs( 6 ) ;
	rhs << 1, 2, 3, 4, 5, 6;

	EXPECT_EQ( expected_1, ( sbm * rhs ) ) ;
	EXPECT_EQ( expected_1, ( sbm.transpose() * rhs ) ) ;

	bogus::SparseBlockMatrix< BlockT > sbm2 = sbm * sbm.transpose() ;

#if EIGEN_VERSION_AT_LEAST(3,1,0)
	typedef bogus::SparseLDLT< double > InvBlockT ;
	bogus::SparseBlockMatrix< InvBlockT > isbm ;
	isbm.cloneStructure( sbm ) ;

	for( unsigned i = 0 ; i < isbm.nBlocks() ; ++i )
	{
		isbm.block(i).compute( sbm.block(i) );
	}

	EXPECT_EQ( expected_2, ( isbm * rhs ) ) ;
#endif
}

TEST( SparseBlock, Scalar )
{

	bogus::SparseBlockMatrix< Eigen::VectorXd > sbm ;
	sbm.setRows( 2, 5 ) ;
	sbm.setCols( 2, 1 ) ;

	sbm.insertBackAndResize( 0, 0 ) << 1, 2, 3, 4, 5 ;
	sbm.insertBackAndResize( 1, 1 ) << 10, 9, 8, 7, 6 ;

	sbm.finalize() ;
	bogus::SparseBlockMatrix< double > sm ;
	sm = sbm.transpose() * sbm ;

	EXPECT_EQ( sbm.block(0,0).squaredNorm(), sm.block(0,0) ) ;
	EXPECT_EQ( sbm.block(1,1).squaredNorm(), sm.block(1,1) ) ;

	Eigen::Vector2d rhs( 1, 1 ), res( 0, 0 ), expected_1 (55, 330 ) ;
	EXPECT_EQ( expected_1, sm * rhs );
	EXPECT_EQ( expected_1, sm.transpose() * rhs );

	sm.splitRowMultiply( 1, rhs, res ) ;
	EXPECT_EQ( Eigen::Vector2d::Zero(), res ) ;

	bogus::SparseBlockMatrix< float, bogus::flags::SYMMETRIC > ssm ;
	ssm = sbm.transpose() * sbm ;
	EXPECT_EQ( expected_1, ssm * rhs );
	EXPECT_EQ( expected_1, ssm.transpose() * rhs );

	EXPECT_EQ( Eigen::Vector2d::Zero(), res ) ;

	ssm.scale( .5 ) ;
	EXPECT_EQ( expected_1, 2 * ( ssm * rhs ) );

}

TEST( SparseBlock, Add )
{
	typedef Eigen::Matrix< double, 3, 4 > BlockT ;
	BlockT sample ;
	sample<< 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, 3, 4 ;

	Eigen::VectorXd rhs( 8 ) ;
	rhs << 4, 2, 1, 0, 7, 6, 9, 8 ;

	{
		bogus::SparseBlockMatrix< BlockT, bogus::flags::COMPRESSED > sbm1 ;
		sbm1.setRows( 4, 3 ) ;
		sbm1.setCols( 2, 4 ) ;

		sbm1.insertBack( 0, 1 ) = BlockT::Ones() ;
		sbm1.insertBack( 1, 0 ) = sample ;
		sbm1.insertBack( 2, 0 ) = 3 * BlockT::Ones() ;
		sbm1.insertBack( 3, 1 ) = 5 * BlockT::Ones() ;

		sbm1.finalize() ;

		const Eigen::VectorXd res1 = sbm1*rhs ;

		bogus::SparseBlockMatrix< BlockT, bogus::flags::COL_MAJOR > sbm2 ;
		sbm2.setRows( 4, 3 ) ;
		sbm2.setCols( 2, 4 ) ;

		sbm2.insertBackOuterInner( 0, 1 ) = BlockT::Ones() ;
		sbm2.insertBackOuterInner( 1, 0 ) = sample ;
		sbm2.insertBackOuterInner( 1, 1 ) = 3 * BlockT::Ones() ;
		sbm2.insertBackOuterInner( 1, 2 ) = 5 * BlockT::Ones() ;

		sbm2.finalize() ;

		const Eigen::VectorXd res2 = sbm2*rhs ;

		sbm1 += sbm2  ;
		EXPECT_EQ( res1+res2, sbm1*rhs ) ;

		sbm1 -= sbm2 ;
		EXPECT_EQ( res1, sbm1*rhs ) ;

	}

	{
		bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::flags::COMPRESSED > sbm ;
		unsigned dims[3] = {3,3,2} ;
		sbm.setRows( 3, dims ) ;
		sbm.setCols( 3, dims ) ;

		sbm.insertBackAndResize( 0, 0 ).setZero() ;
		sbm.insertBackAndResize( 1, 0 ) = sample.block< 3, 3 >( 0, 0 ) ;
		sbm.insertBackAndResize( 1, 1 ).setZero() ;
		sbm.insertBackAndResize( 2, 0 ) = sample.block< 2, 3 >( 1, 1 ) ;
		sbm.insertBackAndResize( 2, 2 ).setZero() ;

		sbm.finalize() ;

		bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::flags::COMPRESSED | bogus::flags::SYMMETRIC > ssbm (sbm) ;
		const Eigen::VectorXd sym_res =  (ssbm*rhs) ;

		sbm += sbm.transpose() ;
		EXPECT_EQ( sym_res, sbm*rhs ) ;

		sbm.add< true >( ssbm ) ;
		EXPECT_EQ( 2*sym_res, sbm*rhs ) ;

		ssbm -= sbm.transpose() ;
		EXPECT_EQ( -sym_res, ssbm*rhs ) ;

		bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::flags::COL_MAJOR > sbm_bis
				= ssbm + sbm.transpose()  ;
		EXPECT_EQ( sym_res, sbm_bis*rhs ) ;

		bogus::SparseBlockMatrix< Eigen::MatrixXd > sbm_ter
				= ( sbm_bis.transpose() - ssbm ).transpose()  ;
		EXPECT_EQ( 2*sym_res, sbm_ter*rhs ) ;

		sbm_ter = -sbm_bis ;
		EXPECT_EQ( -sym_res, sbm_ter*rhs ) ;
		sbm_ter = (-sbm_bis).transpose() ;
		EXPECT_EQ( -sym_res, sbm_ter*rhs ) ;
		sbm_ter = (-2*sbm_bis).transpose() ;
		EXPECT_EQ( -2*sym_res, sbm_ter*rhs ) ;
		sbm_ter = sbm_bis*2 ;
		EXPECT_EQ( sym_res, (.5*sbm_ter)*rhs ) ;
		sbm_ter = sbm_ter -  ( 2 * sbm_bis * 3 ) * .5 ;
		EXPECT_EQ( sym_res, ( rhs.transpose() * ( - sbm_ter ).transpose() ).transpose()  ) ;

		sbm_ter = 2*sbm_bis + ssbm.transpose() + .5*sbm.transpose() - sbm_bis ;
		EXPECT_EQ( sym_res, sbm_ter*rhs ) ;

		sbm_ter = 2 * sbm_bis * ssbm.transpose() + .5*sbm - 2 * sbm_bis.transpose() * ssbm ;
		EXPECT_EQ( sym_res, sbm_ter*rhs ) ;
	}

}
