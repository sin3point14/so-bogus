#include <iostream>

#include "Block.hpp"

#include <Eigen/Core>

using namespace std;

int main()
{
	cout << "Hello World!" << endl;
	typedef Eigen::Matrix< double, 3, 4 > BlockT ;
	bogus::SparseBlockMatrix< BlockT, false > sbm ;
	sbm.setRows( 5, 3 ) ;
	sbm.setCols( 2, 4 ) ;

	sbm.insertBack( 0, 1 ) = BlockT::Ones() ;
	sbm.insertBack( 3, 0 ) = 2 * BlockT::Ones() ;
	sbm.insertBack( 3, 1 ) = 3 * BlockT::Ones() ;

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

	return 0;
}

