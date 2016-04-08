/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_PYRAMIDLAW_HPP
#define BOGUS_PYRAMIDLAW_HPP

#include <cmath>

namespace bogus
{

//! Pyramid local solver that can be used within GaussSeidel and ProjectedGradient solvers
/*!
  For demonstration purposes. Since blocks are 1x1, other libraries are probably more suited.
	\tparam Scalar the scalar type
	\tparam Dimension the dimension of the blocks of the global matrix
  */
template < DenseIndexType Dimension, typename Scalar >
class PyramidLaw
{
public:
	typedef LocalProblemTraits< Dimension, Scalar > Traits ;
	enum{ dimension = Dimension } ;

	//! Constructor
	/*!
	  \param n the size of the global problem ( number of contacts )
	  \param mu array containing the apertures of each second order cone ( friction coefficients )
	  */
	PyramidLaw( const unsigned n, const double * mu ) ;

	//! \return \f$ \vert fb( mu, x, y ) \vert^2_2 \f$, where fb is the SOC Fischer-Burmeister function
	Scalar eval( const unsigned problemIndex,
				 const typename Traits::Vector &x,
				 const typename Traits::Vector &y ) const ;

	bool solveLocal(
			const unsigned problemIndex,
			const typename Traits::Matrix &A,
			const typename Traits::Vector &b,
			typename Traits::Vector &x,
			const Scalar scaling
			) const ;

	//! Projects x on \f$ R^+ \f$
	void projectOnConstraint( const unsigned problemIndex, typename Traits::Vector &x ) const ;

private:

	const double * m_mu ;
	const unsigned m_n ;

} ;


}

#endif
