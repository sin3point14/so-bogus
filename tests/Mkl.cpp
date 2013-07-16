
#include <Core/Block.impl.hpp>
#include <gtest/gtest.h>

#ifdef BOGUS_WITH_MKL

TEST( Mkl, bsr_mv )
{
  typedef Eigen::Matrix< double, 3, 3, Eigen::RowMajor > BlockType ;

    Eigen::VectorXd expected_1(9), expected_2(9), expected_3(9) ;
    expected_1 << 9, 7, 5, 2, 4, 6, 9, 9, 9 ;
    expected_2 << 3, 3, 3, 2, 4, 6, 9, 9, 9 ;
    expected_3 << 9, 7, 5, 0, 0, 0, 9, 9, 9 ;

    Eigen::VectorXd rhs ;

    bogus::SparseBlockMatrix< BlockType, bogus::COMPRESSED > sbm ;
    sbm.setRows( 3 ) ;
    sbm.setCols( 3 ) ;

    sbm.insertBack( 0, 0 ) = BlockType::Ones() ;
    BlockType secBlock ;
    secBlock << 2, 0, 0, 2, 2, 0, 2, 2, 2 ;
    sbm.insertBack( 1, 0 ) = secBlock ;
    sbm.insertBack( 2, 2 ) = 3*BlockType::Ones() ;

    sbm.finalize() ;

    rhs.resize( sbm.rows() ) ;
    rhs.setOnes() ;
    EXPECT_EQ( rhs.rows(), 9 ) ;

    EXPECT_EQ( expected_2, sbm * rhs ) ;
    EXPECT_EQ( expected_3, sbm.transpose() * rhs ) ;

    bogus::SparseBlockMatrix< BlockType, bogus::SYMMETRIC | bogus::COMPRESSED >
            ssbm = sbm ;

    EXPECT_EQ( expected_1, ssbm * rhs ) ;
    EXPECT_EQ( expected_1, ssbm.transpose() * rhs ) ;

    {
        Eigen::MatrixXd mrhs( rhs.rows(), 2 ) ;
        mrhs.col(0) = rhs ;
        mrhs.col(1) = 2 * rhs ;
        Eigen::MatrixXd mres = ssbm * mrhs ;
        EXPECT_EQ( expected_1, mres.col(0) ) ;
        EXPECT_EQ( 2*expected_1, mres.col(1) ) ;
    }
}

#endif

