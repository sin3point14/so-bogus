/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_PRODUCT_GAUSS_SEIDEL_HPP
#define BOGUS_PRODUCT_GAUSS_SEIDEL_HPP

#include "GaussSeidelBase.hpp"
#include "Coloring.hpp"

#include <vector>

namespace bogus
{

//! Matrix-free version of the GaussSeidel iterative solver.
/*!
  Assumes that the system matrix is defines as the product (M M')
  \warning Parallelization is supported, but dangerous. If in doubt, use setMaxThreads(1)
  \sa GaussSeidel
  */
template < typename BlockMatrixType >
class ProductGaussSeidel : public GaussSeidelBase< ProductGaussSeidel, BlockMatrixType >
{
public:
	typedef GaussSeidelBase< bogus::ProductGaussSeidel, BlockMatrixType > Base ;

	typedef typename Base::GlobalProblemTraits GlobalProblemTraits ;
	typedef typename GlobalProblemTraits::Scalar Scalar ;

	//! Default constructor -- you will have to call setMatrix() before using the solve() function
	ProductGaussSeidel( ) : Base() { }
	//! Constructor with the system matrix
	explicit ProductGaussSeidel( const BlockObjectBase< BlockMatrixType > & matrix ) : Base()
	{  setMatrix( matrix ) ; }

	//! Sets the system matrix and initializes internal structures
	ProductGaussSeidel& setMatrix( const BlockObjectBase< BlockMatrixType > & matrix ) ;

	//! Finds an approximate solution for a constrained linear problem
	template < typename NSLaw, typename RhsT, typename ResT >
	Scalar solve( const NSLaw &law, const RhsT &b, ResT &x, bool tryZeroAsWell = true ) const ;

	/*!
	 *  Solves
	 *   y = M x +         p + b
	 *   0 = H'x + C/alpha p + c
	 *   s.t. law(x,y)
	 *  \warning Requires: m_evalEvery multiple of solveEvery ;
	 */
	template < typename NSLaw, typename RhsT, typename ResT, typename LSDerived, typename HDerived >
	Scalar solveWithLinearConstraints( const NSLaw &law,
									   const BlockObjectBase< LSDerived >& Cinv,
									   const BlockObjectBase<  HDerived >& H,
									   const Scalar alpha,
									   const RhsT &b, const RhsT &c,
									   ResT &x,
									   bool tryZeroAsWell = true, unsigned solveEvery = 1) const ;

	/*!
	 *  Solves
	 *   y = M x + W x + b
	 *   s.t. law(x,y)
	 *
	 *  with W arbitrary linear operator ( matrix or expression )
	 *  \warning Requires: m_evalEvery multiple of solveEvery ;
	 */
	template < typename NSLaw, typename RhsT, typename ResT, typename WDerived >
	Scalar solveWithLinearConstraints( const NSLaw &law,
									   const BlockObjectBase < WDerived >& W,
									   const RhsT &b, ResT &x,
									   bool tryZeroAsWell = true, unsigned solveEvery = 1 ) const ;

protected:

	void updateLocalMatrices();

	template < typename NSLaw,  typename VecT, typename ResT >
	void innerLoop (
		bool parallelize, const NSLaw &law, const VecT& b,
		std::vector< unsigned char > &skip, Scalar &ndxRef,
		VecT& Mx, ResT &x	) const ;

	typedef typename Base::Index Index ;

	using Base::m_matrix ;
	using Base::m_maxIters ;
	using Base::m_tol ;
	using Base::m_scaling ;
	using Base::m_maxThreads ;
	using Base::m_evalEvery ;
	using Base::m_skipTol ;
	using Base::m_skipIters ;
	using Base::m_localMatrices ;
	using Base::m_regularization ;

} ;

} //namespace bogus


#endif
