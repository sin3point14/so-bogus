#include <iostream>

#include "Block.hpp"

#include <Eigen/Core>

using namespace std;

int main()
{
	cout << "Hello World!" << endl;

	bogus::SparseBlockMatrix< Eigen::Matrix3d, false > sbm ;
	sbm.setRows( 5, 3 ) ;
	sbm.setCols( 2, 4 ) ;

	sbm.insertBack( 0, 1 ) = Eigen::Matrix3d::Ones() ;
	sbm.insertBack( 3, 0 ) = 2 * Eigen::Matrix3d::Ones() ;
	sbm.insertBack( 3, 1 ) = 3 * Eigen::Matrix3d::Ones() ;

	sbm.finalize();
	cout << sbm << std::endl ;

	return 0;
}

