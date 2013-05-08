#ifndef BOGUS_FRICTION_PROBLEM_HPP
#define BOGUS_FRICTION_PROBLEM_HPP

#include "../Core/Block.hpp"
#include "../Core/BlockSolvers.fwd.hpp"

#include "../Core/Utils/LinearSolverBase.hpp"

#include <Eigen/Core>

namespace bogus
{

template< unsigned Dimension >
struct PrimalFrictionProblem
{
    // Primal Data
    //! M^-1
    SparseBlockMatrix< Eigen::MatrixXd, flags::COMPRESSED  > M ;
    //! E
    bogus::SparseBlockMatrix< Eigen::Matrix< double, Dimension, Dimension >, bogus::flags::COMPRESSED > E ;
    //! H
    typedef Eigen::Matrix< double, Dimension, Eigen::Dynamic > HBlock ;
    SparseBlockMatrix< HBlock > H;

    const double *f ;
    const double *w ;
    const double *mu ;

    // Cached data

    //! M^-1
    SparseBlockMatrix< LU< Eigen::MatrixBase< Eigen::MatrixXd > >, flags::COMPRESSED > MInv ;

    //! M^-1 * H'
    typedef Eigen::Matrix< double, Eigen::Dynamic, Dimension > HtBlock ;
    SparseBlockMatrix< HtBlock, bogus::flags::COL_MAJOR > MInvHt ;

    Eigen::VectorXd MInvf ;
} ;


template< unsigned Dimension >
struct DualFrictionProblem
{
    typedef SparseBlockMatrix< Eigen::Matrix< double, Dimension, Dimension >,
							   flags::SYMMETRIC | flags::COMPRESSED > WType ;
   	typedef GaussSeidel< WType > GaussSeidelType ;

	//! W
    WType W ;

    Eigen::VectorXd b ;
    const double *mu ;

    void computeFrom( PrimalFrictionProblem< Dimension >& primal ) ;

    double solveWith( GaussSeidelType &gs, double * r, const bool staticProblem = false ) const ;

} ;

} //namespace bogus

#endif
