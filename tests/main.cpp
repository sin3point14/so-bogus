#include <iostream>

#include "Block.hpp"

#include <Eigen/Core>

using namespace std;

int main()
{
	cout << "Hello World!" << endl;
	typedef Eigen::MatrixXd BlockT ;
	bogus::SparseBlockMatrix< BlockT > sbm ;
	sbm.setRows( 5, 3 ) ;
	sbm.setCols( 2, 4 ) ;

	sbm.insertBack( 0, 1 ) = BlockT::Ones( 3,4 ) ;
	sbm.insertBack( 3, 0 ) = 2 * BlockT::Ones( 3,4 ) ;
	sbm.insertBack( 3, 1 ) = 3 * BlockT::Ones( 3,4 ) ;

	sbm.finalize();
	cout << sbm << std::endl ;

	Eigen::VectorXd rhs ( sbm.cols() ) ;
	Eigen::VectorXd res ( sbm.rows() ) ;

	rhs.setOnes() ;
	res.setZero() ;
	sbm.multiply( rhs, res ) ;

	std::cout << res.transpose() << std::endl ;
	std::cout << ( sbm*rhs ).transpose() << std::endl ;

	rhs.setZero() ;
	sbm.multiply( res, rhs, true ) ;

	std::cout << rhs.transpose() << std::endl ;
	std::cout << ( res.transpose() * sbm ) << std::endl ;

	sbm.cacheTranspose();
	std::cout << ( res.transpose() * sbm ) << std::endl ;

	bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::BlockMatrixFlags::SYMMETRIC | bogus::BlockMatrixFlags::COMPRESSED > ssbm ;
	ssbm.setRows( 3 ) ;

	ssbm.insertBack( 0, 0 ) = Eigen::Matrix3d::Ones() ;
	ssbm.insertBack( 1, 0 ) = 2*Eigen::Matrix3d::Ones() ;
	ssbm.insertBack( 2, 2 ) = 3*Eigen::Matrix3d::Ones() ;

	ssbm.finalize() ;

	rhs.resize( ssbm.rows() ) ;
	rhs.setOnes() ;

	std::cout << ( ssbm * rhs ).transpose() << std::endl ;

	ssbm.cacheTranspose() ;

	std::cout << ( ssbm * rhs ).transpose() << std::endl ;
	res.resize( 3 ) ;

	for( unsigned k = 0 ; k < 3 ; ++ k )
	{
		res.setZero() ;
		ssbm.splitRowMultiply( k, rhs, res ) ;
		std::cout << res.transpose() << std::endl ;
	}


	bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::BlockMatrixFlags::SYMMETRIC | bogus::BlockMatrixFlags::COMPRESSED > ssbm_copy = ssbm ;
	std::cout << ssbm_copy << std::endl ;

	bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::BlockMatrixFlags::SYMMETRIC > ussbm ;
	ussbm.setRows( 3 ) ;

	ussbm.insertBack( 0, 0 ) = 3*Eigen::Matrix3d::Ones() ;
	ussbm.insertBack( 2, 0 ) = 2*Eigen::Matrix3d::Ones() ;
	ussbm.insertBack( 2, 2 ) = 1*Eigen::Matrix3d::Ones() ;

	ussbm.finalize() ;

	rhs.resize( ussbm.rows() ) ;
	rhs.setOnes() ;

	assert( ussbm.computeMinorIndex() ) ;

	for( unsigned k = 0 ; k < 3 ; ++ k )
	{
		res.setZero() ;
		ussbm.splitRowMultiply( k, rhs, res ) ;
		std::cout << res.transpose() << std::endl ;
	}

	ussbm.cacheTranspose() ;

	for( unsigned k = 0 ; k < 3 ; ++ k )
	{
		res.setZero() ;
		ussbm.splitRowMultiply( k, rhs, res ) ;
		std::cout << res.transpose() << std::endl ;
	}

	bogus::SparseBlockMatrix< Eigen::RowVector4d > copyofsbm ;
	copyofsbm.cloneStructure( sbm ) ;
	std::cout << copyofsbm << std::endl;

	return 0;
}

