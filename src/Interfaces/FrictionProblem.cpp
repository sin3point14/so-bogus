/*
 * This file is part of So-bogus, a C++ sparse block matrix library and
 * Second Order Cone solver.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * So-bogus is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.

 * So-bogus is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with So-bogus.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "FrictionProblem.hpp"

#include "../Core/Block.impl.hpp"
#include "../Core/BlockSolvers/GaussSeidel.impl.hpp"
#include "../Extra/SecondOrder.impl.hpp"

namespace bogus {

template< unsigned Dimension >
void DualFrictionProblem< Dimension >::computeFrom(PrimalFrictionProblem<Dimension> &primal )
{

	// M^-1
	primal.MInv.cloneStructure( primal.M ) ;
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( std::ptrdiff_t i = 0 ; i < (std::ptrdiff_t) primal.M.nBlocks()  ; ++ i )
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

template< unsigned Dimension >
double DualFrictionProblem< Dimension >::evalWith( const GaussSeidelType &gs,
                                               const double *r, const double *u,
									   const bool staticProblem ) const
{
	typedef bogus::SOCLaw< Dimension, double, true  > CoulombLawType	;
	typedef bogus::SOCLaw< Dimension, double, false > SOCLawType	;

	Eigen::Map< const Eigen::VectorXd > r_map ( r, W.rows() ) ;
	Eigen::Map< const Eigen::VectorXd > u_map ( u, W.rows() ) ;

	double res = staticProblem
			? gs.eval( SOCLawType     ( W.rowsOfBlocks(), mu ), u_map, r_map )
			: gs.eval( CoulombLawType ( W.rowsOfBlocks(), mu ), u_map, r_map ) ;

	return res ;
}

template< unsigned Dimension >
double DualFrictionProblem< Dimension >::solveCadoux(GaussSeidelType &gs, double *r, const unsigned cadouxIterations,
        const Signal<unsigned, double> *callback ) const
{
	const std::ptrdiff_t n = W.rowsOfBlocks() ;

	bogus::SOCLaw< Dimension, double, true  > CoulombLaw( n, mu ) ;
	bogus::SOCLaw< Dimension, double, false > SOCLaw	( n, mu ) ;

	gs.setMatrix( W );
	Eigen::Map< Eigen::VectorXd > r_map ( r, W.rows() ) ;

	Eigen::VectorXd s( W.rows() ) ;

	double res = -1 ;
	const double tol = gs.tol() ;
	gs.setTol( 1.e-1 * tol ) ;	//We might experience slow convergence is GS not precise enough

	for( unsigned cdxIter = 0 ; cdxIter < cadouxIterations ; ++cdxIter )
	{
		s = W * r_map + b ;

		res = gs.eval( CoulombLaw, s, r_map ) ;

		if( callback ) callback->trigger( cdxIter, res ) ;
		if( res < tol ) break ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
		for( std::ptrdiff_t i = 0 ; i < n ; ++i )
		{
			s[ Dimension*i ] = s.segment< Dimension >( Dimension*i ).norm() * mu[i] ;
			s.segment< Dimension -1  >( Dimension*i+1 ).setZero() ;
		}

		s += b ;

		gs.solve( SOCLaw, s, r_map, false ) ;

	}

	gs.setTol( tol ) ;

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
