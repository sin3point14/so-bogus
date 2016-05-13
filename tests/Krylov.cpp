
/*
 * Any copyright is dedicated to the Public Domain.
 * http://creativecommons.org/publicdomain/zero/1.0/
*/

#include <iostream>

#include "Core/Block.impl.hpp"
#include "Core/BlockSolvers/Krylov.impl.hpp"

#include "ResidualInfo.hpp"

#include <Eigen/Eigenvalues>

#include <gtest/gtest.h>


TEST( Krylov, CG )
{
	ResidualInfo ri ;

	const Eigen::Vector3d expected_1( .5, .5, .5 ) ;
	const Eigen::Vector3d expected_2( 2, 1, 3 ) ;

	typedef bogus::SparseBlockMatrix< Eigen::Matrix3d > Mat ;
	Mat sbm ;
	sbm.setRows( 1, 3 ) ;
	sbm.setCols( 1, 3 ) ;
	sbm.insertBack(0,0) = Eigen::Vector3d::Constant( 2 ).asDiagonal() ;
	sbm.finalize() ;

	bogus::Krylov< Mat > cg( sbm ) ;
	ri.bindTo( cg.callback() );

	Eigen::Vector3d rhs, res ;
	rhs.setOnes( ) ;

	res.setZero() ;
	ri.setMethodName("CG");
	cg.solve( rhs, res ) ;
	EXPECT_EQ( expected_1, res ) ;

	res.setZero() ;
	ri.setMethodName("BiCG");
	cg.solve_BiCG( rhs, res ) ;
	EXPECT_EQ( expected_1, res ) ;

	res.setZero() ;
	ri.setMethodName("CGS");
	cg.solve_CGS( rhs, res ) ;
	EXPECT_EQ( expected_1, res ) ;

	res.setZero() ;
	ri.setMethodName("BiCGSTAB");
	cg.solve_BiCGSTAB( rhs, res ) ;
	EXPECT_EQ( expected_1, res ) ;

	sbm.block(0) << 0, 1, 0, 1, 0, 0, 0, 0, 1 ;
	rhs << 1, 2, 3 ;

	res.setZero() ;
	ri.setMethodName("BiCG");
	cg.solve( rhs, res, bogus::krylov::BiCG ) ;
	EXPECT_TRUE( expected_2.isApprox( res, 1.e-6 ) ) ;

	res.setZero() ;
	ri.setMethodName("CGS");
	cg.solve_CGS( rhs, res ) ;
	EXPECT_TRUE( expected_2.isApprox( res, 1.e-6 ) ) ;

	res.setOnes() ;
	ri.setMethodName("GMRES");
	cg.solve_GMRES( rhs, res ) ;
	EXPECT_TRUE( expected_2.isApprox( res, 1.e-6 ) ) ;

	res.setOnes() ;
	ri.setMethodName("TFQMR");
	cg.solve_TFQMR( rhs, res ) ;
	EXPECT_TRUE( expected_2.isApprox( res, 1.e-6 ) ) ;

	res.setZero() ;
	cg.setMaxIters( 12 );
	ri.setMethodName("BiCGSTAB");
	cg.solve_BiCGSTAB( rhs, res ) ;
	EXPECT_TRUE( expected_2.isApprox( res, 1.e-6 ) ) ;
}

TEST( Krylov, Preconditioner )
{
	ResidualInfo ri ;

	typedef Eigen::Matrix< double, 5, 5 > Block ;
	Block M_L ;
	M_L <<  1.72194,           0 ,            0,            0,            0,
		   0.027804,      1.78422,            0,            0,            0,
			0.26607,     0.097189,            0,            0,            0,
			0.96157,      0.71437,      0.98738,      1.66828,            0,
		   0.024571,     0.046486,      0.94515,      0.38009,     1.087634 ;

	typedef bogus::SparseBlockMatrix< Block > Mat ;
	Mat sbm ;
	sbm.setRows( 1, 5 ) ;
	sbm.setCols( 1, 5 ) ;
	sbm.insertBack(0,0) = .5 * ( M_L + M_L.transpose() ) ;
	sbm.finalize() ;

	Eigen::Matrix< double, 5, 1 > rhs, res ;
	rhs.setOnes() ;
	double err ;

	res.setZero() ;
	bogus::Krylov< Mat > cg( sbm ) ;
	ri.setMethodName("CG");
	err = cg.solve_CG( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;

	res.setZero() ;
	bogus::Krylov< Mat, bogus::DiagonalPreconditioner > pcg( sbm ) ;
	ri.bindTo( pcg.callback() );

	ri.setMethodName("CG");
	err = pcg.solve_CG( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;

	res.setZero() ;
	ri.setMethodName("BiCG");
	err = pcg.solve_BiCG( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;

	res.setZero() ;
	ri.setMethodName("CGS");
	err = pcg.solve_CGS( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;

	res.setZero() ;
	ri.setMethodName("GMRES");
	err = pcg.solve_GMRES( rhs, res ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;
	EXPECT_GT( 1.e-16, err ) ;

	res.setZero() ;
	pcg.setMaxIters( 40 ) ;
	ri.setMethodName("GMRES(3)");
	err = pcg.asGMRES().setRestart(3).solve( rhs, res ) ;
	EXPECT_GT( 1.e-15, err ) ;
	EXPECT_GT( 1.e-14, ( sbm*res - rhs ).squaredNorm() ) ;

	res.setZero() ;
	ri.setMethodName("TFQMR");
	err = pcg.solve_TFQMR( rhs, res ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;
	EXPECT_GT( 1.e-16, err ) ;
	res.setZero() ;

	res.setZero() ;
	bogus::Krylov< Mat, bogus::DiagonalLUPreconditioner > lucg( sbm ) ;
	lucg.setMaxIters( 1 );
	ri.setMethodName("CG");
	err = lucg.solve_CG( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;

	// Making the block positive definite
	sbm.block(0) *= sbm.block(0) ;

	res.setZero() ;
	bogus::Krylov< Mat, bogus::DiagonalLDLTPreconditioner > ldltcg( sbm ) ;
	ldltcg.setMaxIters( 1 );
	ri.setMethodName("CG");
	err = ldltcg.solve_CG( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;

	res.setZero() ;
	ldltcg.setMaxIters( 10 );
	ri.setMethodName("GMRES");
	err = ldltcg.solve_GMRES( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;

#ifndef BOGUS_BLOCK_WITHOUT_EIGEN_SPARSE
	typedef Eigen::SparseMatrix< double > SparseBlock ;
	typedef bogus::SparseBlockMatrix< SparseBlock > SparseMat ;
	SparseMat ssbm ;
	ssbm.cloneStructure( sbm ) ;
	ssbm.block(0) =  sbm.block(0).sparseView() ;

	res.setZero() ;
	bogus::Krylov< SparseMat > scg( ssbm ) ;
	ri.setMethodName("CG");
	err = scg.solve_CG( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;

#ifdef BOGUS_WITH_EIGEN_STABLE_SPARSE_API

	res.setZero() ;
	bogus::Krylov< SparseMat, bogus::DiagonalPreconditioner > spcg( ssbm ) ;
	ri.setMethodName("CG");
	err = spcg.solve_CG( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;

#ifdef BOGUS_WITH_EIGEN_SPARSE_LDLT
	res.setZero() ;
	bogus::Krylov< SparseMat, bogus::DiagonalLDLTPreconditioner > sldltcg( ssbm ) ;
	sldltcg.setMaxIters( 1 );
	ri.bindTo( sldltcg.callback() );
	ri.setMethodName("CG");
	err = sldltcg.solve_CG( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;
	res.setZero() ;
	ri.setMethodName("BiCGSTAB");
	err = sldltcg.solve_BiCGSTAB( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;
#endif

#ifdef BOGUS_WITH_EIGEN_SPARSE_LU
	res.setZero() ;
	bogus::Krylov< SparseMat, bogus::DiagonalLUPreconditioner > slucg( ssbm ) ;
	slucg.setMaxIters( 1 );
	ri.bindTo( slucg.callback() );
	ri.setMethodName("CG");
	err = slucg.solve_CG( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;
	res.setZero() ;
	ri.setMethodName("BiCGSTAB");
	err = slucg.solve_BiCGSTAB( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;
#endif

#endif

#endif
}

TEST( Krylov, MultipleRhs )
{
	ResidualInfo ri ;

	typedef Eigen::Matrix< double, 10, 10 > Block ;
	Block A ;
	A <<
	1044,  -167, -1016,  -735,   329, -1278,  -454,  1273,  -963, 	-651	,
	  631, -1409,  1633,  1989,  -596,  -493,  -969,    92,  1949, 	939		,
	-1842,  1161, -1341,  1407, -1720, -1099,   516, -1437, -1645, -492		,
	 1741,  1069,  1011, -1593, -1096,  -748,  1584,    28,  1801,	1417	,
	 1951,   -25,  1194,  -608, -1555,   304, -1320,  1022,  1881,	1795	,
	 1535,    22,  1037,  -302,  -339,  -551,   478, -1784,  1826,	-1543	,
	   84, -1653,  1752,  1786,   945,   221, -1767,   347,  1237,	1260	,
	 1259,   774, -1990,  1899,   607, -1085,  -566, -1454,   771,	-1668	,
	 -869,   648,  -864,  -366, -1188,  1193,   -98,  1766,  1381,	-73		,
	 1426, -1462,  -426,  -367,   906,  1747, -1946,  -651,  1336,	-263	;

	const Block rhs = Block::Identity() ;
	Block res ;
	double err ;

	//bogus::Signal< unsigned, double > callback ;
	//callback.connect( &ackCurrentResidual );

	res.setZero() ;
	ri.setMethodName("BiCG");
	err = bogus::krylov::solvers::BiCG< Block >( A, 10 ).solve( rhs, res ) ;
	EXPECT_GT( 1.e-15, err ) ;
	EXPECT_GT( 1.e-15, ( A*res - rhs ).squaredNorm() ) ;

	res.setZero() ;
	ri.setMethodName("CGS");
	err = bogus::krylov::solvers::CGS< Block >( A, 10 ).solve( rhs, res ) ;
	EXPECT_GT( 1.e-15, err ) ;
	EXPECT_GT( 1.e-15, ( A*res - rhs ).squaredNorm() ) ;

	res.setZero() ;
	ri.setMethodName("GMRES");
	err = bogus::krylov::solvers::GMRES< Block >( A, 10 )
			.parallelizeRhs(true).solve( rhs, res ) ;
	EXPECT_GT( 1.e-15, err ) ;
	EXPECT_GT( 1.e-15, ( A*res - rhs ).squaredNorm() ) ;

	res.setZero() ;
	ri.setMethodName("TFQMR");
	err = bogus::krylov::solvers::TFQMR< Block >( A, 20 ).solve( rhs, res ) ;
	EXPECT_GT( 1.e-15, err ) ;
	EXPECT_GT( 1.e-15, ( A*res - rhs ).squaredNorm() ) ;

	typedef bogus::krylov::solvers::GMRES< Block > GMRES ;
	typedef bogus::SparseBlockMatrix< GMRES > Mat ;

	Mat sbm ;
	sbm.setRows( 1, 10 ) ;
	sbm.setCols( 1, 10 ) ;
	sbm.insertBack(0,0) = GMRES( A, 50 ) ;
	sbm.finalize() ;

	bogus::SparseBlockMatrix< Block > mrhs ;
	mrhs.cloneStructure( sbm ) ;
	mrhs.block( 0 ) = rhs ;

	bogus::SparseBlockMatrix< Block > mres =  sbm*mrhs ;
	EXPECT_TRUE ( ( A * mres.block( 0 ) ).isApprox( rhs ) ) ;

}

TEST( Krylov, ProductLhs )
{
	ResidualInfo ri ;

	typedef Eigen::Matrix< double, 5, 5 > Block ;
	Block M_L ;
	M_L <<  1.72194,           0 ,            0,            0,            0,
		   0.027804,      1.78422,            0,            0,            0,
			0.26607,     0.097189,            0,            0,            0,
			0.96157,      0.71437,      0.98738,      1.66828,            0,
		   0.024571,     0.046486,      0.94515,      0.38009,     1.087634 ;


	typedef bogus::SparseBlockMatrix< Block > Mat ;
	Mat sbm ;
	sbm.setRows( 1, 5 ) ;
	sbm.setCols( 1, 5 ) ;
	sbm.insertBack(0,0) = .5 * ( M_L + 2*M_L.transpose() ) ;
	sbm.finalize() ;

	typedef bogus::Product< Mat, bogus::Transpose< Mat > > Prod ;
	Prod prod = sbm * sbm.transpose() ;

	Eigen::Matrix< double, 5, 1 > rhs, res ;
	rhs.setOnes() ;
	double err ;

	res.setZero() ;
	bogus::Krylov< Prod > cg( prod ) ;
	ri.bindTo( cg.callback() );
	ri.setMethodName("CG");
	err = cg.solve_CG( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-8, (prod * res - rhs).lpNorm<Eigen::Infinity>() ) ;
}

