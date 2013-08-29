
#include <Core/Block.impl.hpp>
#include <Core/Block.io.hpp>
#include <gtest/gtest.h>

namespace bogus {

template < typename Matrix, typename RhsT, typename ResT >
void multiply( const Matrix& matrix, const RhsT& rhs, ResT& res )
{
    mv_set< false >( matrix, rhs, res ) ;
}

}


TEST( Map, SparseBlock )
{

  typedef Eigen::MatrixXd BlockT ;

  bogus::MappedSparseBlockMatrix< BlockT > sbm ;


  // Start by creating a normal SparseBlockMatrix
  Eigen::VectorXd expected_1(15), expected_2(8), expected_3(15) ;
  expected_1 << 4, 4, 4, 0, 0, 0, 0, 0, 0, 20, 20, 20, 0, 0, 0 ;
  expected_2 << 120, 120, 120, 120, 192, 192, 192, 192 ;
  expected_3 << 768, 768, 768, 0, 0, 0, 0, 0, 0, 3264, 3264, 3264, 0, 0, 0 ;

  EXPECT_TRUE( bogus::IsTransposable< BlockT >::Value ) ;
  bogus::SparseBlockMatrix< BlockT > origsbm ;
  origsbm.setRows( 5, 3 ) ;
  origsbm.setCols( 2, 4 ) ;

  origsbm.insertBack( 0, 1 ) = BlockT::Ones( 3,4 ) ;
  origsbm.insertBack( 3, 0 ) = 2 * BlockT::Ones( 3,4 ) ;
  origsbm.insertBack( 3, 1 ) = 3 * BlockT::Ones( 3,4 ) ;

  origsbm.finalize();

  // Map it

  sbm.cloneDimensions( origsbm ) ;
  sbm.mapTo(
              origsbm.nBlocks(),
              origsbm.data(),
              origsbm.majorIndex().rowIndex(),
              origsbm.majorIndex().columns()
            );

  EXPECT_EQ( 2u, sbm.blockPtr( 3, 1 ) ) ;
  EXPECT_EQ( sbm.InvalidBlockPtr, sbm.blockPtr( 0, 0 ) ) ;
  EXPECT_EQ( sbm.InvalidBlockPtr, sbm.blockPtr( 1, 1 ) ) ;

  Eigen::VectorXd rhs ( sbm.cols() ) ;
  Eigen::VectorXd res ( sbm.rows() ) ;

  rhs.setOnes() ;
  res.setZero() ;
  sbm.multiply< false >( rhs, res ) ;


  Eigen::MatrixXd nah ;
  nah.setZero( sbm.rows(), sbm.cols() ) ;
  bogus::multiply( nah, rhs, res ) ;
  bogus::multiply( sbm, rhs, res ) ;

  EXPECT_EQ( expected_1, res ) ;
  EXPECT_EQ( expected_1, sbm*rhs ) ;

  rhs.setZero() ;
  sbm.multiply< true >( res, rhs ) ;

  EXPECT_EQ( expected_2, rhs ) ;
  EXPECT_EQ( expected_2, ( res.transpose() * sbm ).transpose() ) ;
  EXPECT_EQ( expected_2, ( sbm.transpose() * res )  ) ;
  EXPECT_EQ( expected_3, ( rhs.transpose() * sbm.transpose() ).transpose()  ) ;

  EXPECT_FALSE( sbm.minorIndex().valid );
  EXPECT_FALSE( sbm.transposeIndex().valid );

  sbm.computeMinorIndex() ;
  EXPECT_TRUE( sbm.minorIndex().valid );
  EXPECT_EQ( expected_3, ( rhs.transpose() * sbm.transpose() ).transpose()  ) ;

  bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::SYMMETRIC > prod( sbm * sbm.transpose() );
  EXPECT_EQ( sbm * expected_2, prod * res ) ;
  EXPECT_EQ( sbm.transpose() * expected_3, ( sbm.transpose() * sbm ) * rhs ) ;
  bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::COL_MAJOR > sum( sbm.transpose() + sbm.transpose() );
  EXPECT_EQ( 2 * expected_2, sum * res ) ;
  EXPECT_EQ( 2 * expected_3, 2 * sbm * rhs ) ;

  // Symmetric
  bogus::MappedSparseBlockMatrix< Eigen::MatrixXd, bogus::SYMMETRIC > mprod ;
  mprod.mapTo( prod ) ;
  EXPECT_EQ( sbm * expected_2, mprod * res ) ;
  EXPECT_EQ( sbm * expected_2, mprod.transpose() * res ) ;
  EXPECT_EQ( prod * ( sbm * expected_2 ), ( mprod * mprod ) * res ) ;

  // Col-major
  bogus::MappedSparseBlockMatrix< Eigen::MatrixXd, bogus::COL_MAJOR > msum( sum ) ;
  EXPECT_EQ( 2 * expected_2, msum * res ) ;

}


