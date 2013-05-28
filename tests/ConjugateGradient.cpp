#include "Core/Block.impl.hpp"
#include "Core/BlockSolvers.impl.hpp"

#include <Eigen/Eigenvalues>

#include <gtest/gtest.h>

static void ackCurrentResidual( unsigned GSIter, double err )
{
	std::cout << "CG: " << GSIter << " ==> " << err << std::endl ;
}

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
	//cg.callback().connect( &ackCurrentResidual );

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

TEST( ConjugateGradient, Preconditioner )
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
	bogus::ConjugateGradient< Mat > cg( sbm ) ;
	err = cg.solve( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;

	res.setZero() ;
	bogus::ConjugateGradient< Mat, bogus::DiagonalPreconditioner > pcg( sbm ) ;
	err = pcg.solve( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;

	res.setZero() ;
	err = pcg.solve_BiCG( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;

	res.setZero() ;
	bogus::ConjugateGradient< Mat, bogus::DiagonalLUPreconditioner > lucg( sbm ) ;
	lucg.setMaxIters( 1 );
	err = lucg.solve( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;

	// Making the block positive definite
	sbm.block(0) *= sbm.block(0) ;

	res.setZero() ;
	bogus::ConjugateGradient< Mat, bogus::DiagonalLDLTPreconditioner > ldltcg( sbm ) ;
	ldltcg.setMaxIters( 1 );
	err = ldltcg.solve( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;

#if EIGEN_VERSION_AT_LEAST(3,1,0)
	typedef Eigen::SparseMatrix< double > SparseBlock ;
	typedef bogus::SparseBlockMatrix< SparseBlock > SparseMat ;
	SparseMat ssbm ;
	ssbm.setRows( 1, 5 ) ;
	ssbm.setCols( 1, 5 ) ;
	ssbm.insertBack(0,0) =  sbm.block(0).sparseView() ;
	ssbm.finalize() ;

	res.setZero() ;
	bogus::ConjugateGradient< SparseMat > scg( ssbm ) ;
	err = scg.solve( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;

	res.setZero() ;
	bogus::ConjugateGradient< SparseMat, bogus::DiagonalPreconditioner > spcg( ssbm ) ;
	err = spcg.solve( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;

	res.setZero() ;
	bogus::ConjugateGradient< SparseMat, bogus::DiagonalLDLTPreconditioner > sldltcg( ssbm ) ;
	sldltcg.setMaxIters( 1 );
	err = sldltcg.solve( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
	res.setZero() ;
	err = sldltcg.solve_BiCGSTAB( rhs, res ) ;
	EXPECT_GT( 1.e-16, err ) ;
#endif
}

