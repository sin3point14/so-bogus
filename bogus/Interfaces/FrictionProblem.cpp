#include "FrictionProblem.hpp"

#include "../Core/Block.impl.hpp"
#include "../Core/BlockSolvers/GaussSeidel.impl.hpp"
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

	//W
	W = primal.H * ( primal.MInv * primal.H.transpose() ) ;

	W.cacheTranspose() ;

	// M^-1 f, b
	b = primal.E.transpose() * Eigen::VectorXd::Map( primal.w, primal.H.rows())
			- primal.H * ( primal.MInv * Eigen::VectorXd::Map( primal.f, primal.H.cols() ) );

	mu = primal.mu ;
}


template< unsigned Dimension >
double DualFrictionProblem< Dimension >::solveWith( GaussSeidelType &gs, double *r,
									   const bool staticProblem ) const
{
	typedef bogus::SOCLaw< Dimension, double, true  > CoulombLawType	;
	typedef bogus::SOCLaw< Dimension, double, false > SOCLawType	;

	gs.setMatrix( W );
	Eigen::Map< Eigen::VectorXd > r_map ( r, W.rows() ) ;

	double res = staticProblem
			? gs.solve( SOCLawType     ( W.rowsOfBlocks(), mu ), b, r_map )
			: gs.solve( CoulombLawType ( W.rowsOfBlocks(), mu ), b, r_map ) ;

	return res ;
}

#ifdef BOGUS_INSTANTIATE_2D_SOC
template struct DualFrictionProblem< 2u > ;
template struct PrimalFrictionProblem< 2u > ;
#endif

#ifdef BOGUS_INSTANTIATE_3D_SOC
template struct DualFrictionProblem< 3u > ;
template struct PrimalFrictionProblem< 3u > ;
#endif

#ifdef BOGUS_INSTANTIATE_DYNAMIC_SOC
template struct DualFrictionProblem< Eigen::Dynamic > ;
template struct PrimalFrictionProblem< Eigen::Dynamic > ;
#endif

} //namespace bogus
