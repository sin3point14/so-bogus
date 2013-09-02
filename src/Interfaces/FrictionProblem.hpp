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

#ifndef BOGUS_FRICTION_PROBLEM_HPP
#define BOGUS_FRICTION_PROBLEM_HPP

#include "../Core/Block.hpp"
#include "../Core/BlockSolvers.fwd.hpp"

namespace bogus
{

template< unsigned Dimension >
struct PrimalFrictionProblem
{
	// Primal Data
	//! M -- mass matrix
	SparseBlockMatrix< Eigen::MatrixXd > M ;
	//! E -- local rotation matrix ( world <-> contact basis )
	bogus::SparseBlockMatrix< Eigen::Matrix< double, Dimension, Dimension > > E ;

	typedef Eigen::Matrix< double, Dimension, Eigen::Dynamic > HBlock ;
	//! H -- deformation gradient ( generalized coordinates <-> 3D world )
	SparseBlockMatrix< HBlock, UNCOMPRESSED > H;

	//! External forces
	const double *f ;
	//! Free velocity ( such that u = Hv + w )
	const double *w ;
	//! Coulomb friction coefficients
	const double *mu ;

	// Cached data

	//! M^-1
	SparseBlockMatrix< LU< Eigen::MatrixBase< Eigen::MatrixXd > > > MInv ;

} ;


template< unsigned Dimension >
struct DualFrictionProblem
{
	typedef SparseBlockMatrix< Eigen::Matrix< double, Dimension, Dimension, Eigen::RowMajor >,
							   SYMMETRIC > WType ;
	typedef GaussSeidel< WType > GaussSeidelType ;

	//! W -- Delassus operator
	WType W ;

	//! Rhs ( such that u = Wr + b )
	Eigen::VectorXd b ;

	//! Coulomb friction coefficients
	const double *mu ;

	//! Computes this DualFrictionProblem from the given \p primal
	void computeFrom( PrimalFrictionProblem< Dimension >& primal ) ;

	//! Solves this problem
	/*!
	  \param gs The GaussSeidel< WType > solver to use
	  \param r  Both the initial guess and the result
	  \param staticProblem If true, solve this problem as a \b SOCQP instead of a Coulomb Friction problem
	  \returns the error as returned by the GaussSeidel::solve() function
	  */
	double solveWith( GaussSeidelType &gs, double * r, const bool staticProblem = false ) const ;

} ;

} //namespace bogus

#endif
