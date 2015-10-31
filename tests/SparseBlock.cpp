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


TEST( SparseBlock, MatrixVector )
{
		//Eigen::initParallel()
		Eigen::MatrixXf A = Eigen::MatrixXf::Zero(1,1); A = A*A;

	Eigen::VectorXd expected_1(15), expected_2(8), expected_3(15) ;
	expected_1 << 4, 4, 4, 0, 0, 0, 0, 0, 0, 20, 20, 20, 0, 0, 0 ;
	expected_2 << 120, 120, 120, 120, 192, 192, 192, 192 ;
	expected_3 << 768, 768, 768, 0, 0, 0, 0, 0, 0, 3264, 3264, 3264, 0, 0, 0 ;

	typedef Eigen::MatrixXd BlockT ;
	EXPECT_TRUE( bogus::IsTransposable< BlockT >::Value ) ;
	bogus::SparseBlockMatrix< BlockT, bogus::flags::UNCOMPRESSED > sbm ;
	sbm.setRows( 5, 3 ) ;
	sbm.setCols( 2, 4 ) ;

	sbm.insertBack( 0, 1 ) = BlockT::Ones( 3,4 ) ;
	sbm.insertBack( 3, 0 ) = 2 * BlockT::Ones( 3,4 ) ;
	sbm.insertBack( 3, 1 ) = 3 * BlockT::Ones( 3,4 ) ;

	sbm.finalize();
	//std::cout << sbm << std::endl ;

//	for( int row = 0 ; row < sbm.rowsOfBlocks() ; ++ row ) {
//		for( bogus::SparseBlockMatrix< BlockT, bogus::flags::UNCOMPRESSED >::InnerIterator it ( sbm.innerIterator( row ) ) ; it ; ++it ) {
//			std::cout << "Block at (" << row << ", " << it.inner() << ") is \n "
//					  << sbm.block( it.ptr() ) << "\n" ;
//		}
//	}

	EXPECT_EQ( 2u, sbm.blockPtr( 3, 1 ) ) ;
	EXPECT_EQ( sbm.InvalidBlockPtr, sbm.blockPtr( 0, 0 ) ) ;
	EXPECT_EQ( sbm.InvalidBlockPtr, sbm.blockPtr( 1, 1 ) ) ;

	Eigen::VectorXd rhs ( sbm.cols() ) ;
	Eigen::VectorXd res ( sbm.rows() ) ;

	rhs.setOnes() ;
	res.setZero() ;
	sbm.multiply< false >( rhs, res ) ;

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

TEST( SparseBlock, Symmetric )
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



TEST( SparseBlock, MMult )
{
	typedef Eigen::Matrix< double, 3, 4 > BlockT ;
	BlockT sample ;
	sample<< 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, 3, 4 ;

	Eigen::VectorXd rhs( 12 ), expected_1( 8 ) ;
	rhs << 4, 2, 1, 0, 7, 6, 9, 8, 2, 4, 6, 8 ;

	expected_1 << 30410, 31550, 32690, 33830, 46674, 46674, 46674, 46674 ;

	{
		bogus::SparseBlockMatrix< BlockT > sbm ;
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

		bogus::SparseBlockMatrix< Eigen::Matrix4d, bogus::flags::UNCOMPRESSED | bogus::flags::SYMMETRIC > mm_t = sbm.transpose()*sbm ;
		ASSERT_EQ( 3u, mm_t.nBlocks() ) ;
		EXPECT_TRUE( mm_t.diagonal(1).isApprox( Eigen::Matrix4d::Constant( 78 ) ) ) ;
		mm_t.setFromProduct< false > ( sbm.transpose()*sbm ) ;
		ASSERT_EQ( 3u, mm_t.nBlocks() ) ;
		EXPECT_TRUE( mm_t.diagonal(1).isApprox( Eigen::Matrix4d::Constant( 78 ) ) ) ;

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
		bogus::SparseBlockMatrix< BlockT, bogus::flags::COL_MAJOR | bogus::flags::UNCOMPRESSED > sbm ;
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

		bogus::SparseBlockMatrix< Eigen::Matrix4d, bogus::flags::COL_MAJOR | bogus::flags::UNCOMPRESSED | bogus::flags::SYMMETRIC > mm_t = sbm.transpose()*sbm ;
		ASSERT_EQ( mm_t.nBlocks(), 3u ) ;
		EXPECT_TRUE( mm_t.diagonal(1).isApprox( Eigen::Matrix4d::Constant( 78 ) ) );

		mm_t.setFromProduct< false > ( sbm.transpose()*sbm ) ;
		ASSERT_EQ( mm_t.nBlocks(), 3u ) ;
		EXPECT_TRUE( mm_t.diagonal(1).isApprox( Eigen::Matrix4d::Constant( 78 ) ) );

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

	typedef bogus::LU< Eigen::MatrixBase< Eigen::MatrixXd > > BlockT ;
	EXPECT_FALSE( bogus::IsTransposable< BlockT >::Value ) ;
	bogus::SparseBlockMatrix< BlockT > isbm ;
	isbm.setRows( 2, 3 ) ;
	isbm.setCols( 2, 3 ) ;
	isbm.insertBack( 0, 0 ).compute( 4 * Eigen::Matrix3d::Identity() ) ;
	isbm.insertBack( 1, 1 ).compute( 8 * Eigen::Matrix3d::Identity() ) ;
	isbm.finalize() ;

	Eigen::VectorXd rhs( 6 ) ;
	rhs.setOnes() ;

	EXPECT_EQ( expected_2, ( isbm * rhs ) ) ;
	EXPECT_EQ( expected_2, ( isbm * rhs ) ) ;

}

TEST( SparseBlock, Sparse )
{
	Eigen::VectorXd expected_1( 6 ) ;
	expected_1 << 2, 4, 6, 16, 40, 6 ;
	Eigen::VectorXd expected_2( 6 ) ;
	expected_2 << 0.5, 1, 1.5, 1, 0.625, 6 ;

	typedef Eigen::SparseMatrix< double > BlockT ;
	EXPECT_TRUE( bogus::IsTransposable< BlockT >::Value ) ;
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

#ifdef BOGUS_WITH_EIGEN_SPARSE_LDLT
	{

		typedef bogus::SparseLDLT< double > InvBlockT ;
		EXPECT_FALSE( bogus::IsTransposable< InvBlockT >::Value ) ;
		bogus::SparseBlockMatrix< InvBlockT > isbm ;
		isbm.cloneStructure( sbm ) ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
		for( std::ptrdiff_t i = 0 ; i < (std::ptrdiff_t) isbm.nBlocks() ; ++i )
		{
			isbm.block(i).compute( sbm.block(i) );
		}

		EXPECT_EQ( expected_2, ( isbm * rhs ) ) ;
	}
#endif

#ifdef BOGUS_WITH_EIGEN_SPARSE_LU
	{
		typedef bogus::SparseLU< double > InvBlockT ;
		EXPECT_FALSE( bogus::IsTransposable< InvBlockT >::Value ) ;
		bogus::SparseBlockMatrix< InvBlockT > isbm ;
		isbm.cloneStructure( sbm ) ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
		for( std::ptrdiff_t i = 0 ; i < (std::ptrdiff_t) isbm.nBlocks() ; ++i )
		{
			isbm.block(i).compute( sbm.block(i) );
		}

		EXPECT_EQ( expected_2, ( isbm * rhs ) ) ;
	}
#endif
}

TEST( SparseBlock, Scalar )
{
	ASSERT_TRUE( bogus::IsTransposable< double >::Value ) ;
	ASSERT_TRUE( bogus::IsTransposable< float >::Value ) ;
	ASSERT_TRUE( bogus::IsTransposable< int >::Value ) ;
	ASSERT_TRUE( bogus::IsTransposable< char >::Value ) ;

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

	Eigen::Vector2d rhs( 1, 1 ), expected_1 (55, 330 ) ;
	Eigen::Matrix< double, 1, 1 > res1d ; res1d.setZero() ;

	EXPECT_TRUE( expected_1.isApprox( sm * rhs ) );
	EXPECT_TRUE( expected_1.isApprox( sm.transpose() * rhs ) );


	sm.splitRowMultiply( 1, rhs, res1d ) ;
	EXPECT_TRUE( res1d.isZero() ) ;

	bogus::SparseBlockMatrix< float, bogus::flags::SYMMETRIC > ssm ;
	Eigen::Vector2f frhs( 1, 1 ) ;
	ssm = sbm.transpose() * sbm ;
	EXPECT_TRUE( expected_1.cast< float >().isApprox( ssm * frhs ) );
	EXPECT_TRUE( expected_1.cast< float >().isApprox( ssm.transpose() * frhs ) );

	ssm.scale( .5 ) ;
	EXPECT_TRUE( expected_1.cast< float >().isApprox( 2 * ( ssm * frhs ) ) );

}

TEST( SparseBlock, Add )
{
	typedef Eigen::Matrix< double, 3, 4 > BlockT ;
	BlockT sample ;
	sample<< 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, 3, 4 ;

	Eigen::VectorXd rhs( 8 ) ;
	rhs << 4, 2, 1, 0, 7, 6, 9, 8 ;

	{
		bogus::SparseBlockMatrix< BlockT > sbm1 ;
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

		sbm2.insertByOuterInner< true >( 0, 1 ) = BlockT::Ones() ;
		sbm2.insertByOuterInner< true >( 1, 0 ) = sample ;
		sbm2.insertByOuterInner< true >( 1, 1 ) = 3 * BlockT::Ones() ;
		sbm2.insertByOuterInner< true >( 1, 2 ) = 5 * BlockT::Ones() ;

		sbm2.finalize() ;

		const Eigen::VectorXd res2 = sbm2*rhs ;

		sbm1 += sbm2  ;
		EXPECT_EQ( res1+res2, sbm1*rhs ) ;

		sbm1 -= sbm2 ;
		EXPECT_EQ( res1, sbm1*rhs ) ;

	}

	{
		bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::flags::UNCOMPRESSED > sbm ;
		unsigned dims[3] = {3,3,2} ;
		sbm.setRows( 3, dims ) ;
		sbm.setCols( 3, dims ) ;

		sbm.insertBackAndResize( 0, 0 ).setZero() ;
		sbm.insertBackAndResize( 1, 0 ) = sample.block< 3, 3 >( 0, 0 ) ;
		sbm.insertBackAndResize( 1, 1 ).setZero() ;
		sbm.insertBackAndResize( 2, 0 ) = sample.block< 2, 3 >( 1, 1 ) ;
		sbm.insertBackAndResize( 2, 2 ).setZero() ;

		sbm.finalize() ;

		bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::flags::SYMMETRIC > ssbm (sbm) ;
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

		const Eigen::VectorXd mm_res =  (sbm_bis * ( ssbm.transpose() * rhs )) ;
		sbm_ter.setFromProduct< false >( sbm_bis * ssbm.transpose() ) ;
		EXPECT_EQ( mm_res, sbm_ter * rhs ) ;
		sbm_ter.setFromProduct< true  >( sbm_bis * ssbm.transpose() ) ;
		EXPECT_EQ( mm_res, sbm_ter * rhs ) ;
		sbm_ter.setFromProduct< false >( ssbm.transpose() * sbm_bis ) ;
		EXPECT_EQ( mm_res, sbm_ter * rhs ) ;
		sbm_ter.setFromProduct< true  >( ssbm.transpose() * sbm_bis ) ;
		EXPECT_EQ( mm_res, sbm_ter * rhs ) ;

		sbm_ter = 2 * sbm_bis * ssbm.transpose() + .5*sbm - 2 * sbm_bis.transpose() * ssbm ;
		EXPECT_EQ( sym_res, sbm_ter*rhs ) ;
	}

}

TEST( SparseBlock, YoDawg )
{
	typedef bogus::SparseBlockMatrix< Eigen::Matrix3d > BlockType ;
	EXPECT_TRUE( bogus::IsTransposable< BlockType >::Value ) ;

	bogus::SparseBlockMatrix< BlockType > sbm ;

	sbm.setRows( 2, 6 ) ;
	sbm.setCols( 2, 6 ) ;
	sbm.reserve( 3 ) ;

	BlockType &b00 = sbm.insertBack( 0, 0 ) ;
	BlockType &b10 = sbm.insertBack( 1, 0 ) ;
	BlockType &b11 = sbm.insertBack( 1, 1 ) ;

	sbm.finalize() ;

	b10.setRows( 2, 3 );
	b10.setCols( 2, 3 );
	b10.insertBack( 0, 0 ).setIdentity() ;
	b10.insertBackAndResize( 1, 0 ) << 1, 2, 3, 4, 5, 6, 7, 8, 9 ;
	b10.insertBack( 1, 1 ).setOnes() ;
	b10.finalize() ;

	b00 = b10 * b10.transpose() ;
	b11 = .5 * ( b10 + b10.transpose() ) ;

	Eigen::VectorXd rhs(12), expected_1(12), expected_2(12) ;
	rhs << 1, 2, 3, 4, 5, 0, 6, 5, 4, 3, 2, 1 ;
	expected_1 << 25 , 35 , 45 , 257 , 572 , 887 , 16 , 19 , 22 , 43 , 83.5 , 124 ;
	expected_2 << 49 ,  64 ,  79 , 263 , 578 , 893 ,  15 ,  17 ,  19 ,  20 , 42.5 ,  65 ;

	const Eigen::VectorXd res1 = sbm * rhs  ;
	EXPECT_EQ( res1, expected_1 ) ;
	const Eigen::VectorXd res2 = sbm.transpose() * rhs ;
	EXPECT_EQ( res2, expected_2 ) ;

	EXPECT_EQ( res1 + res2, ( sbm + sbm.transpose() ) * rhs ) ;
	EXPECT_EQ( sbm*res2, ( sbm * sbm.transpose() ) * rhs ) ;
	EXPECT_EQ( 2*res2, ( 2* sbm ).transpose() * rhs ) ;

	bogus::SparseBlockMatrix< BlockType, bogus::SYMMETRIC > sbm2 = sbm + sbm.transpose() ;
	bogus::SparseBlockMatrix< BlockType, bogus::SYMMETRIC > sbm3 = sbm * sbm.transpose() ;
	EXPECT_EQ( res1 + res2, sbm2 * rhs ) ;
	EXPECT_EQ( sbm*res2, sbm3 * rhs ) ;
}

TEST( SparseBlock, Permutation)
{
	bogus::SparseBlockMatrix< Eigen::Matrix2i > sbm ;

	sbm.setCols( 8 );
	sbm.setRows( 8 );
	std::size_t perm[8] = { 3,1,4,2,0,7,6,5 } ;

	for( unsigned i = 0 ; i < 8 ; ++ i )
	{
		for( unsigned j = 0 ; j < 8 ; ++ j )
		{
			sbm.insertBack( i, j ) << i, 0, 0, j ;
		}
	}
	sbm.finalize();

	sbm.applyPermutation( &perm[0] ) ;

	for( unsigned i = 0 ; i < 8 ; ++ i )
	{
		for( unsigned j = 0 ; j < 8 ; ++ j )
		{
			ASSERT_EQ( sbm.block( i, j )( 0, 0 ), (int)perm[i] ) ;
			ASSERT_EQ( sbm.block( i, j )( 1, 1 ), (int)perm[j] ) ;
		}
	}

	sbm += sbm.transpose() ;

	bogus::SparseBlockMatrix< Eigen::Matrix2i, bogus::SYMMETRIC > ssbm ( sbm ) ;
	ssbm.applyPermutation( &perm[0] ) ;

	for( unsigned i = 0 ; i < 8 ; ++ i )
	{
		for( unsigned j = 0 ; j <= i ; ++ j )
		{
			ASSERT_EQ( ssbm.block( i, j )( 0, 0 ), (int)(perm[perm[i]] + perm[perm[j]]) ) ;
			ASSERT_EQ( ssbm.block( i, j )( 1, 1 ), (int)(perm[perm[j]] + perm[perm[i]]) ) ;
		}
	}
}
