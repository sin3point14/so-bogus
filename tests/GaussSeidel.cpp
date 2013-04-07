
#include "Block.hpp"
#include "GaussSeidel.hpp"
#include "SecondOrder.hpp"

#include <gtest/gtest.h>

TEST( GaussSeidel, Small )
{
    typedef bogus::Coulomb3D Law ;
	typedef bogus::SparseBlockMatrix< Eigen::Matrix3d > Mat ;
	Law law ;
	Mat mat ;

	bogus::GaussSeidel< Mat > gs( mat ) ;
	Eigen::VectorXd x, b ;

	double res = gs.solve( law, b, x ) ;
	std::cout << res << std::endl ;
}
