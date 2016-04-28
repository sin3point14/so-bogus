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

#include "FrictionProblem.impl.hpp"

#include "../Core/BlockSolvers/GaussSeidel.impl.hpp"
#include "../Core/BlockSolvers/ProjectedGradient.impl.hpp"
#include "../Core/BlockSolvers/ADMM.impl.hpp"


namespace bogus {

// PrimalFrictionProblem

template< unsigned Dimension >
void PrimalFrictionProblem< Dimension >::computeMInv( )
{
	// M^-1
	MInv.cloneStructure( M ) ;
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( std::ptrdiff_t i = 0 ; i < (std::ptrdiff_t) M.nBlocks()  ; ++ i )
	{
		MInv.block(i).compute( M.block(i) ) ;
	}
}


template< unsigned Dimension >
double PrimalFrictionProblem< Dimension >::solveWith( ADMMType &admm, double lambda, double* v, double * r ) const
{
	const Eigen::VectorXd f = Eigen::VectorXd::Map( this->f, M.rows() ) ;
	const Eigen::VectorXd w = E.transpose() * Eigen::VectorXd::Map( this->w, H.rows() ) ;

	Eigen::VectorXd::MapType r_map( r, H.rows() ) ;
	Eigen::VectorXd::MapType v_map( v, H.cols() ) ;

	bogus::QuadraticProxOp< MInvType > prox( MInv, lambda, f ) ;

	Eigen::ArrayXd inv_mu = 1./Eigen::ArrayXd::Map( mu, H.rowsOfBlocks()  );
	bogus::SOCLaw< Dimension, double, false > law( inv_mu.rows(), inv_mu.data() ) ;

	admm.setMatrix( H ) ;
	return admm.solve( law, prox, w, v_map, r_map ) ;
}


template< unsigned Dimension >
double PrimalFrictionProblem< Dimension >::solveWith( DualAMAType &dama, double* v, double * r, const bool staticProblem ) const
{
	const Eigen::VectorXd f = Eigen::VectorXd::Map( this->f, M.rows() ) ;
	const Eigen::VectorXd w = E.transpose() * Eigen::VectorXd::Map( this->w, H.rows() ) ;

	Eigen::VectorXd::MapType r_map( r, H.rows() ) ;
	Eigen::VectorXd::MapType v_map( v, H.cols() ) ;

	Eigen::ArrayXd inv_mu = 1./Eigen::ArrayXd::Map( mu, H.rowsOfBlocks()  );
	bogus::SOCLaw< Dimension, double, false > law( inv_mu.rows(), inv_mu.data() ) ;

	dama.setMatrix( H ) ;

	if( staticProblem ) {
		bogus::SOCLaw< Dimension, double, false > law( H.rowsOfBlocks(), mu ) ;
		return dama.solve( law, M, f, w, v_map, r_map ) ;
	} else {
		bogus::SOCLaw< Dimension, double, true  > law( H.rowsOfBlocks(), mu ) ;
		return dama.solve( law, M, f, w, v_map, r_map ) ;
	}

}


// DualFrictionProblem

template< unsigned Dimension >
void DualFrictionProblem< Dimension >::computeFrom(const PrimalFrictionProblem<Dimension> &primal )
{

	//W
	W = primal.H * ( primal.MInv * primal.H.transpose() ) ;

	// M^-1 f, b
	b = primal.E.transpose() * Eigen::VectorXd::Map( primal.w, primal.H.rows())
			- primal.H * ( primal.MInv * Eigen::VectorXd::Map( primal.f, primal.H.cols() ) );

	mu = Eigen::VectorXd::Map( primal.mu, W.rowsOfBlocks() ) ;
}

template< unsigned Dimension >
double DualFrictionProblem< Dimension >::solveWith( GaussSeidelType &gs, double *r,
										 const bool staticProblem ) const
{
	gs.setMatrix( W );

	return friction_problem::solve( *this, gs, r, staticProblem ) ;
}

template< unsigned Dimension >
double DualFrictionProblem< Dimension >::solveWith( ProjectedGradientType &pg,
													double *r ) const
{
	pg.setMatrix( W );

	return friction_problem::solve( *this, pg, r, true ) ;
}

template< unsigned Dimension >
double DualFrictionProblem< Dimension >::evalWith( const GaussSeidelType &gs,
													 const double *r,
													 const bool staticProblem ) const
{
	return friction_problem::eval( *this, gs, r, staticProblem ) ;
}

template< unsigned Dimension >
double DualFrictionProblem< Dimension >::evalWith( const ProjectedGradientType &gs,
													 const double *r ) const
{
	return friction_problem::eval( *this, gs, r, true) ;
}


template< unsigned Dimension >
double DualFrictionProblem< Dimension >::solveCadoux(GaussSeidelType &gs, double *r, const unsigned cadouxIterations,
		const Signal<unsigned, double> *callback ) const
{
	gs.setMatrix( W );

	return friction_problem::solveCadoux( *this, gs, r, cadouxIterations, callback ) ;
}

template< unsigned Dimension >
double DualFrictionProblem< Dimension >::solveCadoux(ProjectedGradientType &pg, double *r, const unsigned cadouxIterations,
		const Signal<unsigned, double> *callback ) const
{
	pg.setMatrix( W );

	return friction_problem::solveCadoux( *this, pg, r, cadouxIterations, callback ) ;
}

template< unsigned Dimension >
double DualFrictionProblem< Dimension >::solveCadouxVel(GaussSeidelType &gs, double *u, const unsigned cadouxIterations,
		const Signal<unsigned, double> *callback ) const
{
	gs.setMatrix( W );

	return friction_problem::solveCadouxVel( *this, gs, u, cadouxIterations, callback ) ;
}

template< unsigned Dimension >
double DualFrictionProblem< Dimension >::solveCadouxVel(ProjectedGradientType &pg, double *u, const unsigned cadouxIterations,
		const Signal<unsigned, double> *callback ) const
{
	pg.setMatrix( W );

	return friction_problem::solveCadouxVel( *this, pg, u, cadouxIterations, callback ) ;
}

template< unsigned Dimension >
void DualFrictionProblem< Dimension >::applyPermutation(
		const std::vector<std::size_t> &permutation)
{
	assert( !permuted() ) ;

	m_permutation = permutation ;

	m_invPermutation.resize( m_permutation.size() );
	for( std::size_t i = 0 ; i < m_permutation.size() ; ++i )
		m_invPermutation[ m_permutation[i] ] = i ;

	W.applyPermutation( data_pointer(m_permutation) ) ;
	friction_problem::applyPermutation< Dimension >( m_permutation, b, W.colOffsets() ) ;
	bogus::applyPermutation( m_permutation.size(), data_pointer(m_permutation), mu ) ;
}

template< unsigned Dimension >
void DualFrictionProblem< Dimension >::undoPermutation()
{
	if( !permuted() )
		return ;

	W.applyPermutation( data_pointer(m_invPermutation) ) ;
	friction_problem::applyPermutation< Dimension >( m_invPermutation, b, W.colOffsets() ) ;
	bogus::applyPermutation( m_invPermutation.size(), data_pointer(m_invPermutation), mu ) ;

	m_permutation.clear() ;
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
