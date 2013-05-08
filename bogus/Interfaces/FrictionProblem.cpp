#include "FrictionProblem.hpp"

#include "../Core/Block.impl.hpp"
#include "../Core/BlockSolvers.impl.hpp"
#include "../Core/SecondOrder.impl.hpp"

namespace bogus {

template< unsigned Dimension >
void DualFrictionProblem< Dimension >::computeFrom(PrimalFrictionProblem<Dimension> &primal )
{

    // M^-1
    primal.MInv.cloneStructure( primal.M ) ;
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
    for( int i = 0 ; i < (int) primal.M.nBlocks()  ; ++ i )
    {
        primal.MInv.block(i).compute( primal.M.block(i) ) ;
    }

    // M^-1 * H'
    primal.MInvHt = primal.MInv * primal.H.transpose() ;

    //W
    W = primal.H * primal.MInvHt ;

    W.cacheTranspose() ;

    // M^-1 f, b
    primal.MInvf = primal.MInv * Eigen::VectorXd::Map( primal.f, primal.H.cols() ) ;
    b = ( primal.E.transpose() * Eigen::VectorXd::Map( primal.w, Dimension*primal.H.rowsOfBlocks()) )
            - primal.H * ( primal.MInvf );

    mu = primal.mu ;
}


template< unsigned Dimension >
double DualFrictionProblem< Dimension >::solveWith( GaussSeidel<WType> &gs, double *r,
                                       const bool staticProblem ) const
{
	gs.setMatrix( W );
    Eigen::Map< Eigen::VectorXd > r_map ( r, W.rows() ) ;

    double res = staticProblem
            ? gs.solve( bogus::SOCLaw< Dimension, double, false > ( W.rowsOfBlocks(), mu ), b, r_map )
            : gs.solve( bogus::SOCLaw< Dimension, double, true  > ( W.rowsOfBlocks(), mu ), b, r_map ) ;

    return res ;
}

template struct DualFrictionProblem< 2u > ;
template struct DualFrictionProblem< 3u > ;
template struct PrimalFrictionProblem< 2u > ;
template struct PrimalFrictionProblem< 3u > ;

} //namespace bogus
