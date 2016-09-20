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
		mm_t_mt -= mm_t  * sbm.transpose() ;
		        Eigen::VectorXd one_minus_one = (mm_t_mt*rhs) ;
		EXPECT_TRUE( one_minus_one.isZero() ) ;

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

	{
		const unsigned N = 3 ;

		typedef bogus::SparseBlockMatrix< BlockT > SBM ;
		SBM sbm[3] ;

		for( unsigned i = 0 ; i < N ; ++i ) {
			sbm[i].setRows( 1, 3 ) ;
			sbm[i].setCols( 2, 4 ) ;

			sbm[i].insertBack( 0, 0 ) = BlockT::Ones() ;
			sbm[i].insertBack( 0, 1 ) = i*sample ;

			sbm[i].finalize() ;
		}

		Eigen::VectorXd rhs1( sbm[0].rows() ) ;
		rhs1.setOnes() ;

		Eigen::VectorXd res1( rhs1.rows() ) ;

		typedef bogus::Product< SBM, bogus::Transpose<SBM> > Prod ;
		bogus::NarySum< Prod > sum( sbm[0].rows(), sbm[0].rows() )  ;

		sum.multiply< false >( rhs1, res1, 1, 0) ;
		ASSERT_EQ( 0, res1.norm() ) ;

		for( unsigned i = 0 ; i < N ; ++i ) {
			sum += i * ( sbm[i] * sbm[i].transpose() ) ;
			res1 += i * ( sbm[i] * sbm[i].transpose() ) * rhs1 ;
		}

		ASSERT_TRUE( res1.isApprox( sum * rhs1 ) ) ;
		ASSERT_TRUE( res1.isApprox( sum.transpose() * rhs1 ) ) ;


		bogus::NarySum< Prod > sum2 ( sbm[0] * sbm[0].transpose() )  ;
		sum2 += sum ;
		sum2 -= sbm[0] * sbm[0].transpose() ;
		ASSERT_TRUE( res1.isApprox( sum2 * rhs1 ) ) ;

		sum2 -= sum ;
		bogus::SparseBlockMatrix< Eigen::Matrix3d > sumBM = sum2 + sum.transpose() ;
		ASSERT_TRUE( res1.isApprox( sumBM * rhs1 ) ) ;

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

	typedef bogus::SparseBlockMatrix< BlockType, bogus::SYMMETRIC >	SymSBM ;
	ASSERT_TRUE( (bogus::Addition< SymSBM, bogus::Transpose< SymSBM > >::is_self_transpose) ) ;
	ASSERT_FALSE( (bogus::Product< SymSBM, bogus::Transpose< SymSBM > >::is_self_transpose) ) ;

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

TEST(SparseBlock, FixedSize)
{
	bogus::SparseBlockMatrix< Eigen::Matrix< double, 2, 4 > > sbm(1,1) ;

	ASSERT_EQ( 1, sbm.rowsOfBlocks() ) ;
	ASSERT_EQ( 1, sbm.colsOfBlocks() ) ;
	ASSERT_EQ( 2, sbm.rows() ) ;
	ASSERT_EQ( 4, sbm.cols() ) ;

	sbm.insertBack(0,0) << 1,2,3,4,5,6,7,8 ;
	sbm.finalize() ;


	bogus::SparseBlockMatrix< Eigen::Matrix< double, 2, 4 > > copy = sbm ;
	ASSERT_EQ( sbm.block(0), copy.block(0) ) ;

	bogus::SparseBlockMatrix< Eigen::Matrix< double, 4, 2 > > tsbm = sbm.transpose() ;
	ASSERT_EQ( sbm.block(0), tsbm.block(0).transpose() ) ;

}


