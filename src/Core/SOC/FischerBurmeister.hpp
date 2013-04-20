/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver          
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>                       
 *
 * This Source Code Form is subject to the terms of the Mozilla Public 
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef BOGUS_FISCHER_BURMEISTER_HPP
#define BOGUS_FISCHER_BURMEISTER_HPP

#include "../Utils/NumTraits.hpp"
#include "../Utils/EigenMatrixTraits.hpp"

namespace bogus {

template< unsigned Dimension, typename Scalar >
struct FBBaseFunction
{
	typedef MatrixTraits< Dimension, Scalar > Traits ;
	typedef typename Traits::Vector Vector ;
	typedef typename Traits::Matrix Matrix ;

	static void compute( const Scalar mu, const Vector& x, const Vector& y, Vector& fb ) ;

	static void computeJacobian(
				const Scalar mu, const Vector& x, const Vector& y,
			Vector& fb, Matrix& dFb_dx, Matrix& dFb_dy ) ;

private:
	template <bool JacobianAsWell >
	static void compute(
			const Vector& x, const Vector& y, Vector& fb,
			Matrix& dFb_dx, Matrix& dFb_dy ) ;

} ;

template< unsigned Dimension, typename Scalar, bool DeSaxceCOV >
class FischerBurmeister
{

public:
  typedef MatrixTraits< Dimension, Scalar > Traits ;
  typedef FBBaseFunction< Dimension, Scalar > BaseFunction ;
  typedef typename Traits::Vector Vector ;
  typedef typename Traits::Matrix Matrix ;

  FischerBurmeister(
	const Scalar mu,
	const Matrix& A,
	const Vector& b,
	const Scalar scaling )
	  : m_mu( mu ), m_scaling( scaling ), m_A( A ), m_b( b )
  {}

  void compute( const Vector& x, Vector& fb ) const ;
  void computeJacobian( const Vector& x, Vector& fb, Matrix& dFb_dx ) const ;

  static void compute( const Scalar mu, const Vector& x, const Vector& y, Vector& fb ) ;

private:
  Scalar m_mu ;
  Scalar m_scaling ;
  const Matrix& m_A ;
  const Vector& m_b ;

} ;



}

#endif
