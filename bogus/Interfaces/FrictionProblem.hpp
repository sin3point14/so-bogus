#ifndef BOGUS_FRICTION_PROBLEM_HPP
#define BOGUS_FRICTION_PROBLEM_HPP

#include "../Core/Block.hpp"
#include "../Core/BlockSolvers.fwd.hpp"

#include "../Core/Utils/LinearSolverBase.hpp"

#include <Eigen/Core>

namespace bogus
{

struct PrimalFrictionProblem
{
    // Primal Data
    //! M^-1
    SparseBlockMatrix< Eigen::MatrixXd, flags::COMPRESSED  > M ;
    //! E
    bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::COMPRESSED > E ;
    //! H
    typedef Eigen::Matrix< double, 3, Eigen::Dynamic > HBlock ;
    SparseBlockMatrix< HBlock > H;

    const double *f ;
    const double *w ;
    const double *mu ;

    // Cached data

    //! M^-1
    SparseBlockMatrix< LU< Eigen::MatrixBase< Eigen::MatrixXd > >, flags::COMPRESSED > MInv ;

    //! M^-1 * H'
    typedef Eigen::Matrix< double, Eigen::Dynamic, 3 > HtBlock ;
    SparseBlockMatrix< HtBlock, bogus::flags::COL_MAJOR > MInvHt ;

    Eigen::VectorXd MInvf ;
} ;


struct DualFrictionProblem
{
    //! W
    typedef bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::SYMMETRIC | bogus::flags::COMPRESSED > WType ;
    WType W ;

    Eigen::VectorXd b ;

    const double *mu ;

    void computeFrom( PrimalFrictionProblem& primal ) ;

    double solveWith(const GaussSeidel< WType > &gs, const bool staticProblem, double * r ) const ;

} ;

} //namespace bogus

#endif
