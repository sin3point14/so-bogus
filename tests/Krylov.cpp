
#include <iostream>

#include "Core/Block.impl.hpp"
#include "Core/BlockSolvers/Krylov.impl.hpp"

#include <Eigen/Eigenvalues>

#include <gtest/gtest.h>

const char* g_meth  ;

static void ackCurrentResidual( unsigned GSIter, double err )
{
	std::cout << g_meth << " " << GSIter << " ==> " << err << std::endl ;
}

TEST( Krylov, CG )
{

	const Eigen::Vector3d expected_1( .5, .5, .5 ) ;
	const Eigen::Vector3d expected_2( 2, 1, 3 ) ;

	typedef bogus::SparseBlockMatrix< Eigen::Matrix3d > Mat ;
	Mat sbm ;
	sbm.setRows( 1, 3 ) ;
	sbm.setCols( 1, 3 ) ;
	sbm.insertBack(0,0) = Eigen::Vector3d::Constant( 2 ).asDiagonal() ;
	sbm.finalize() ;

	bogus::Krylov< Mat > cg( sbm ) ;
	cg.callback().connect( &ackCurrentResidual );

	Eigen::Vector3d rhs, res ;
	rhs.setOnes( ) ;

	res.setZero() ;
	g_meth = "CG" ;
	cg.solve( rhs, res ) ;
	EXPECT_EQ( expected_1, res ) ;

	res.setZero() ;
	g_meth = "BiCG" ;
	cg.solve_BiCG( rhs, res ) ;
	EXPECT_EQ( expected_1, res ) ;

	res.setZero() ;
	g_meth = "CGS" ;
	cg.solve_CGS( rhs, res ) ;
	EXPECT_EQ( expected_1, res ) ;

	res.setZero() ;
	g_meth = "BiCGSTAB" ;
	cg.solve_BiCGSTAB( rhs, res ) ;
	EXPECT_EQ( expected_1, res ) ;

	sbm.block(0) << 0, 1, 0, 1, 0, 0, 0, 0, 1 ;
	rhs << 1, 2, 3 ;

	res.setZero() ;
	g_meth = "BiCG" ;
	cg.solve( rhs, res, bogus::krylov::BiCG ) ;
	EXPECT_TRUE( expected_2.isApprox( res, 1.e-6 ) ) ;

	res.setZero() ;
	g_meth = "CGS" ;
	cg.solve_CGS( rhs, res ) ;
	EXPECT_TRUE( expected_2.isApprox( res, 1.e-6 ) ) ;

	res.setOnes() ;
	g_meth = "GMRES" ;
	cg.solve_GMRES( rhs, res ) ;
	EXPECT_TRUE( expected_2.isApprox( res, 1.e-6 ) ) ;

	res.setOnes() ;
	g_meth = "TFQMR" ;
	cg.solve_TFQMR( rhs, res ) ;
	EXPECT_TRUE( expected_2.isApprox( res, 1.e-6 ) ) ;

	res.setZero() ;
	cg.setMaxIters( 12 );
	g_meth = "BiCGSTAB" ;
	cg.solve_BiCGSTAB( rhs, res ) ;
	EXPECT_TRUE( expected_2.isApprox( res, 1.e-6 ) ) ;
}

TEST( Krylov, Preconditioner )
{

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
	g_meth = "CG" ;
	err = cg.solve_CG( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;

	res.setZero() ;
	bogus::Krylov< Mat, bogus::DiagonalPreconditioner > pcg( sbm ) ;
	pcg.callback().connect( &ackCurrentResidual );

	g_meth = "CG" ;
	err = pcg.solve_CG( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;

	res.setZero() ;
	g_meth = "BiCG" ;
	err = pcg.solve_BiCG( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;

	res.setZero() ;
	g_meth = "CGS" ;
	err = pcg.solve_CGS( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;

	res.setZero() ;
	g_meth = "GMRES" ;
	err = pcg.solve_GMRES( rhs, res ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;
	EXPECT_GT( 1.e-16, err ) ;

	res.setZero() ;
	pcg.setMaxIters( 30 ) ;
	g_meth = "GMRES(3)" ;
	err = pcg.solve_GMRES( rhs, res, 3 ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-12, ( sbm*res - rhs ).squaredNorm() ) ;

	res.setZero() ;
	g_meth = "TFQMR" ;
	err = pcg.solve_TFQMR( rhs, res ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;
	EXPECT_GT( 1.e-16, err ) ;
	res.setZero() ;

	res.setZero() ;
	bogus::Krylov< Mat, bogus::DiagonalLUPreconditioner > lucg( sbm ) ;
	lucg.setMaxIters( 1 );
	g_meth = "CG" ;
	err = lucg.solve_CG( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;

	// Making the block positive definite
	sbm.block(0) *= sbm.block(0) ;

	res.setZero() ;
	bogus::Krylov< Mat, bogus::DiagonalLDLTPreconditioner > ldltcg( sbm ) ;
	ldltcg.setMaxIters( 1 );
	g_meth = "CG" ;
	err = ldltcg.solve_CG( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;

	res.setZero() ;
	ldltcg.setMaxIters( 10 );
	g_meth = "GMRES" ;
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
	g_meth = "CG" ;
	err = scg.solve_CG( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;

#ifdef BOGUS_WITH_EIGEN_STABLE_SPARSE_API

	res.setZero() ;
	bogus::Krylov< SparseMat, bogus::DiagonalPreconditioner > spcg( ssbm ) ;
	g_meth = "CG" ;
	err = spcg.solve_CG( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;

#ifdef BOGUS_WITH_EIGEN_SPARSE_LDLT
	res.setZero() ;
	bogus::Krylov< SparseMat, bogus::DiagonalLDLTPreconditioner > sldltcg( ssbm ) ;
	sldltcg.setMaxIters( 1 );
	sldltcg.callback().connect( &ackCurrentResidual );
	g_meth = "CG" ;
	err = sldltcg.solve_CG( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;
	res.setZero() ;
	g_meth = "BiCGSTAB" ;
	err = sldltcg.solve_BiCGSTAB( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;
#endif

#ifdef BOGUS_WITH_EIGEN_SPARSE_LU
	res.setZero() ;
	bogus::Krylov< SparseMat, bogus::DiagonalLUPreconditioner > slucg( ssbm ) ;
	slucg.setMaxIters( 1 );
	slucg.callback().connect( &ackCurrentResidual );
	g_meth = "CG" ;
	err = slucg.solve_CG( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;
	res.setZero() ;
	g_meth = "BiCGSTAB" ;
	err = slucg.solve_BiCGSTAB( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	EXPECT_GT( 1.e-16, ( sbm*res - rhs ).squaredNorm() ) ;
#endif

#endif

#endif
}

