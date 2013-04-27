#include "Core/Block.impl.hpp"
#include "Core/BlockSolvers.impl.hpp"

#include <gtest/gtest.h>

TEST( ConjugateGradient, CG )
{
	const Eigen::Vector3d expected_1( .5, .5, .5 ) ;
	const Eigen::Vector3d expected_2( 2, 1, 3 ) ;

	typedef bogus::SparseBlockMatrix< Eigen::Matrix3d > Mat ;
	Mat sbm ;
	sbm.setRows( 1, 3 ) ;
	sbm.setCols( 1, 3 ) ;
	sbm.insertBack(0,0) = Eigen::Vector3d::Constant( 2 ).asDiagonal() ;
	sbm.finalize() ;

	bogus::ConjugateGradient< Mat > cg( sbm ) ;

	Eigen::Vector3d rhs, res ;
	rhs.setOnes( ) ;

	res.setZero() ;
	cg.solve( rhs, res ) ;
	EXPECT_EQ( expected_1, res ) ;

	res.setZero() ;
	cg.solve_BiCG( rhs, res ) ;
	EXPECT_EQ( expected_1, res ) ;

	res.setZero() ;
	cg.solve_BiCGSTAB( rhs, res ) ;
	EXPECT_EQ( expected_1, res ) ;

	sbm.block(0) << 0, 1, 0, 1, 0, 0, 0, 0, 1 ;
	rhs << 1, 2, 3 ;

	res.setZero() ;
	cg.solve_BiCG( rhs, res ) ;
	EXPECT_TRUE( expected_2.isApprox( res, 1.e-6 ) ) ;

	res.setZero() ;
	cg.setMaxIters( 12 );
    cg.solve_BiCGSTAB( rhs, res ) ;
	EXPECT_TRUE( expected_2.isApprox( res, 1.e-6 ) ) ;
}
