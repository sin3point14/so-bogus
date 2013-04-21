
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/array.hpp>

#include <Core/Eigen/EigenSerialization.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include <gtest/gtest.h>

TEST( Serialization, Eigen )
{

	Eigen::MatrixXd m (3,3) ;
	m << 1, 2, 3, 4, 5, 6, 7, 8, 9 ;
	Eigen::Matrix3d m33 (3,3) ;
	m33 << 2, 4, 6, 8, 5, 6, 7, 8, 9 ;
	Eigen::VectorXd v (3) ;
	v << 1, 2, 3 ;
	Eigen::RowVectorXd rv (3) ;
	rv << 1, 2, 3 ;

	{
		std::ofstream ofs("/tmp/bogus_serialization_test");
		boost::archive::text_oarchive oa(ofs);
		oa << m << m33 << v << rv;
	}

	Eigen::MatrixXd m_ ;
	Eigen::Matrix3d m33_ ;
	Eigen::VectorXd v_;
	Eigen::RowVectorXd rv_;
	{
		std::ifstream ifs("/tmp/bogus_serialization_test");
		boost::archive::text_iarchive ia(ifs);
		ia >> m_ >> m33_ >> v_ >> rv_ ;
	}

	ASSERT_EQ( m, m_ ) ;
	ASSERT_EQ( m33, m33_ ) ;
	ASSERT_EQ( v, v_ ) ;
	ASSERT_EQ( rv, rv_ ) ;
}

#if EIGEN_VERSION_AT_LEAST(3,1,0)
TEST( Serialization, EigenSparse )
{
	Eigen::SparseMatrix< double > sm ( 3, 3 ) ;
	sm.reserve ( 3 );
	sm.insert(0,0) = 1 ;
	sm.insert(1,1) = 2 ;
	sm.insert(2,2) = 3 ;
	sm.makeCompressed() ;

	{
		std::ofstream ofs("/tmp/bogus_serialization_test");
		boost::archive::text_oarchive oa(ofs);
		oa << sm;
	}

	Eigen::SparseMatrix< double > sm_;
	{
		std::ifstream ifs("/tmp/bogus_serialization_test");
		boost::archive::text_iarchive ia(ifs);
		ia >> sm_ ;
	}

	ASSERT_EQ( Eigen::Vector3d( 1., 2., 3.), sm_ * Eigen::Vector3d::Ones() ) ;
}
#endif
