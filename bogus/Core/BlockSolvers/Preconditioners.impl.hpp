/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver          
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>                       
 *
 * This Source Code Form is subject to the terms of the Mozilla Public 
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "Preconditioners.hpp"
#include "../Utils/NumTraits.hpp"

#ifndef BOGUS_PRECONDITIONERS_IMPL_HPP
#define BOGUS_PRECONDITIONERS_IMPL_HPP

namespace bogus {

template < typename BlockMatrixType >
DiagonalPreconditioner< BlockMatrixBase< BlockMatrixType > >::DiagonalPreconditioner( const BlockMatrixBase< BlockMatrixType > &matrix )
{
	typedef BlockMatrixTraits< BlockMatrixType > Traits ;
	typedef typename Traits::BlockType BlockType ;
	typedef typename Traits::Index Index ;
	typedef typename MatrixTraits< BlockType >::Scalar Scalar ;

	m_diagonal.resize( matrix.rows() ) ;
	for( Index i = 0, cur_row = 0 ; i < matrix.rowsOfBlocks() ; ++i )
	{
		const BlockType &b = matrix.diagonal( i ) ;
		m_diagonal.segment( cur_row, b.rows() ) = b.diagonal() ;
		for( Index k = 0 ; k < b.rows() ; ++k, ++cur_row )
		{
			if( NumTraits< Scalar >::isSquareZero( m_diagonal[ cur_row ] ) )
			{
				m_diagonal[ cur_row ] = 1 ;
			}
		}
	}
}

template < typename BlockMatrixType >
template < bool transpose, typename ResT, typename RhsT >
void DiagonalPreconditioner< BlockMatrixBase< BlockMatrixType > >::apply( const RhsT& rhs, ResT &res ) const
{
	res = rhs.cwiseQuotient( m_diagonal ) ;
}


} //naemspace bogus

#endif

