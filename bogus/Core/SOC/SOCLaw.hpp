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

	//! Constructor
	/*!
	  \param n the size of the global problem ( number of contacts )
	  \param mu array containing the apertures of each second order cone ( friction coefficients )
	  */
	SOCLaw( const unsigned n, const double * mu ) ;

	//! \return \f$ \frac 1 {1 + n} \vert fb( mu, x, y ) \vert^2_2 \f$, where fb is the SOC Fischer-Burmeister function
	template< typename VectorT, typename OtherVectorT >
	Scalar eval( const VectorT &x, const OtherVectorT &y ) const
	{
		typedef FischerBurmeister< Traits::dimension, typename Traits::Scalar, DeSaxceCOV > FBFunction ;

		assert( (unsigned) x.rows() == m_n * Traits::dimension ) ;
		assert( (unsigned) y.rows() == m_n * Traits::dimension ) ;

		Scalar sum = 0. ;
		typename Traits::Vector lx, ly, fb ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for private( lx, ly, fb ) reduction ( + : sum )
#endif
		for( int i = 0 ; i < (int) m_n ; ++ i )
		{
			lx = Traits::segment( i, x ) ;
			ly = Traits::segment( i, y ) ;
			FBFunction::compute( m_mu[i], lx, ly, fb ) ;
			sum += fb.squaredNorm() ;
		}

		return sum / ( 1 + m_n );
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
