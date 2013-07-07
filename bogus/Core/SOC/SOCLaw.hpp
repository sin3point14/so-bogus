/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_SOCLAW_HPP
#define BOGUS_SOCLAW_HPP

#include "../BlockSolvers.fwd.hpp"
#include "FischerBurmeister.hpp"
#include "LocalSOCSolver.hpp"

#include <vector>

namespace bogus
{

//! Non-smooth laws based on Second Order Cone complementarity. To be used within as the first argument to GaussSeidel::solve().
/*!
	\tparam Dimension the dimension of the local problem. Specializations exist form dimension 2 and 3.
	\tparam Scalar the scalar type
	\tparam DeSaxceCOV Whether to perform the \cite DSF98 change of variable when solving the local problem.
	Should be \c true for modeling Coulomb friction, or \c false for standard SOC complementarity. \sa solveLocal()
	\tparam Strat local_soc_solver::Strategy for solving the local problems. Unavailable for dimensions other than 2 and 3.
  */
template < unsigned Dimension, typename Scalar, bool DeSaxceCOV,
		   local_soc_solver::Strategy Strat >
class SOCLaw
{
public:
	typedef LocalProblemTraits< Dimension, Scalar > Traits ;
	enum{ dimension = Dimension } ;

	//! Constructor
	/*!
	  \param n the size of the global problem ( number of contacts )
	  \param mu array containing the apertures of each second order cone ( friction coefficients )
	  */
	SOCLaw( const unsigned n, const double * mu ) ;

	//! \return \f$ \vert fb( mu, x, y ) \vert^2_2 \f$, where fb is the SOC Fischer-Burmeister function
	Scalar eval( const unsigned problemIndex,
				 const typename Traits::Vector &x,
				 const typename Traits::Vector &y ) const
	{
		typedef FischerBurmeister< Traits::dimension, typename Traits::Scalar, DeSaxceCOV > FBFunction ;

		typename Traits::Vector fb( x.rows() ) ;
		FBFunction::compute( m_mu[problemIndex], x, y, fb ) ;

		return fb.squaredNorm() ;
	}

	//! Solves the local problem
	/*!
	  \f[
		\left\{
		  \begin{array}{rcl}
			y &=& \mathrm{ <DS> } \left( A x + b \right ) \\
			K_{ \frac 1 \mu } \ni y & \perp & x \in K_{ \mu }
		  \end{array}
		\right.
	  \f]
	  where \f$ \mu \f$ is \c m_mu[\p problemIndex] and \c <DS> is the optional De Saxce change of variable.

	  That is, if \p DeSaxceCOV is false then \c <DS> is the identity function, otherwise
	  \f[ \mathrm{ <DS> }( x ) := x + \mu \vert x_T \vert \left( 1, 0, ... \right)^{\top} \f]
	  \param scaling Used as a scaling factor for \p x when calculating the error function
	*/
	bool solveLocal(
			const unsigned problemIndex,
			const typename Traits::Matrix &A,
			const typename Traits::Vector &b,
			typename Traits::Vector &x,
			const Scalar scaling
			) const ;

private:

	const double * m_mu ;
	const unsigned m_n ;
	Scalar m_localTol ;

} ;

//! Predefined non-smooth law for 2D Coulomb friction
typedef SOCLaw< 2u, double,  true > Coulomb2D ;
//! Predefined non-smooth law for 3D Coulomb friction
typedef SOCLaw< 3u, double,  true > Coulomb3D ;
//! Predefined non-smooth law for 2D SOC complementarity
typedef SOCLaw< 2u, double, false > SOC2D ;
//! Predefined non-smooth law for 3D SOC complementarity
typedef SOCLaw< 3u, double, false > SOC3D ;

}

#endif // SOCLAW_HPP
