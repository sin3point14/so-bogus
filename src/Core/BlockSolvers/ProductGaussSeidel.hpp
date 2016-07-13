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

namespace block_solvers_impl {

//! Wrapper for the block-diagonal matrix of the ProductGaussSeidel solver
/*! Specializations allows to:
 *   - Use the Identity matrix by default without requiring the user to
 *     call any function or specify template parameters
 *   - Specify a pointer to an arbitrary matrix when the third template
 *     parameter of ProductGaussSeidel is explicitely set
 */
template < typename MainMatrixType, typename Type>
struct DiagonalMatrixWrapper {

	const Type* ptr ;

	DiagonalMatrixWrapper()
	    : ptr( BOGUS_NULL_PTR(const Type) ) {}
	DiagonalMatrixWrapper( const Type& diag)
	    : ptr(&diag) { }

	bool valid() const { return ptr; }
	const Type& get() const {
		return *ptr ;
	}
} ;
template< typename MainMatrixType >
struct DiagonalMatrixWrapper < MainMatrixType, typename MainMatrixType::Scalar >
{
	typedef typename MainMatrixType::Scalar Scalar;
	ConstantArray<Scalar> ptr ;
	DiagonalMatrixWrapper( const Scalar& s = 1)
	    : ptr(s) { }

	bool valid() const { return true; }
	const ConstantArray<Scalar>& get() const {
		return ptr ;
	}
} ;
} //block_solvers_impl

//! Matrix-free version of the GaussSeidel iterative solver.
/*!
  Assumes that the system matrix is defined as the product (M D M'),
  with D a block-diagonal matrix whose block sizes coincide with those of
  the columns of M.
  \warning Parallelization is supported, but dangerous. If in doubt, use setMaxThreads(1)
  \note Works best when D and the columns of M are quite sparse
  \sa GaussSeidel
  */
template < typename BlockMatrixType, typename DiagonalType = typename BlockMatrixType::Scalar >
class ProductGaussSeidel
        : public GaussSeidelBase< ProductGaussSeidel<BlockMatrixType, DiagonalType>, BlockMatrixType >
{
public:
	typedef GaussSeidelBase< ProductGaussSeidel, BlockMatrixType > Base ;

	typedef typename Base::GlobalProblemTraits GlobalProblemTraits ;
	typedef typename GlobalProblemTraits::Scalar Scalar ;

	typedef block_solvers_impl::DiagonalMatrixWrapper<BlockMatrixType, DiagonalType> DiagWrapper ;

	//! Default constructor -- you will have to call setMatrix() before using the solve() function
	ProductGaussSeidel()
	    : Base()
	{ }
	//! Constructor with the system matrix
	explicit ProductGaussSeidel( const BlockObjectBase< BlockMatrixType > & matrix )
	    : Base()
	{  setMatrix( matrix ) ; }
	//! Constructor with both
	ProductGaussSeidel( const BlockObjectBase< BlockMatrixType > & matrix, const DiagonalType& diagonal )
	    : Base(), m_diagonal( DiagWrapper( diagonal ) )
	{  setMatrix( matrix ) ; }

	//! Sets the system matrix (M) and initializes internal structures
	ProductGaussSeidel& setMatrix( const BlockObjectBase< BlockMatrixType > & matrix ) ;

	//! Sets the system diagonal (D) and initializes internal structures
	ProductGaussSeidel& setDiagonal( const DiagonalType &diagonal )
	{ m_diagonal = DiagWrapper( diagonal ) ; }

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

	DiagWrapper m_diagonal ;

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
