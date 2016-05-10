/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_BLOCK_GAUSS_SEIDEL_BASE_IMPL_HPP
#define BOGUS_BLOCK_GAUSS_SEIDEL_BASE_IMPL_HPP

#include "GaussSeidelBase.hpp"
#include "ConstrainedSolverBase.impl.hpp"

#include "../Utils/Threads.hpp"

namespace bogus
{

template < template <typename> class GaussSeidelImpl, typename BlockMatrixType >
void GaussSeidelBase< GaussSeidelImpl, BlockMatrixType >::updateLocalMatrices( )
{

	Base::updateScalings() ;

	if( !m_matrix )
		return ;

	const Index n = m_matrix->rowsOfBlocks() ;
	m_localMatrices.resize( n ) ;
	m_regularization.resize( n ) ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( Index i = 0 ; i <  n ; ++i )
	{
		const typename BlockMatrixType::BlockPtr ptr = explicitMatrix().diagonalBlockPtr( i ) ;

		if( ptr == BlockMatrixType::InvalidBlockPtr ) {
			resize( m_localMatrices[i], m_matrix->blockRows(i), m_matrix->blockCols(i) ) ;
			set_zero( m_localMatrices[i] ) ;
		} else {
			m_localMatrices[i] = GlobalProblemTraits::asConstMatrix( explicitMatrix().block( ptr ) ) ;
		}

		if( m_autoRegularization > 0. )
		{
			m_regularization(i) = std::max( 0., m_autoRegularization - m_localMatrices[i].eigenvalues().real().minCoeff() ) ;
			m_localMatrices[i].diagonal().array() += m_regularization(i) ;
		} else m_regularization(i) = 0. ;

	}
}

template < template <typename> class GaussSeidelImpl, typename BlockMatrixType >
template < typename NSLaw, typename RhsT, typename ResT >
typename GaussSeidelBase< GaussSeidelImpl, BlockMatrixType >::Scalar
GaussSeidelBase< GaussSeidelImpl, BlockMatrixType >::evalAndKeepBest(
		const NSLaw &law, const RhsT &b, const ResT &x,
		typename GlobalProblemTraits::DynVector& y,
		typename GlobalProblemTraits::DynVector& x_best, Scalar &err_best ) const
{

	y = b ;
	m_matrix->template multiply< false >( x, y, 1, 1 ) ;
	const Scalar err = Base::eval( law, y, x ) ;

	if( err < err_best )
	{
		x_best = x ;
		err_best = err ;
	}

	return err ;
}

template < template <typename> class GaussSeidelImpl, typename BlockMatrixType >
template < typename NSLaw, typename RhsT, typename ResT >
void GaussSeidelBase< GaussSeidelImpl, BlockMatrixType >::tryZero(
		const NSLaw &law, const RhsT &b, ResT &x,
		typename GlobalProblemTraits::DynVector& x_best, Scalar &err_best ) const
{
	x.setZero() ;
	const Scalar err_zero = Base::eval( law, b, x ) ;

	if( err_zero < err_best )
	{
		err_best = err_zero ;
		x_best.setZero() ;
	} else {
		x = x_best ;
	}
}


} //namespace bogus


#endif
