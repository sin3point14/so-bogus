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
class DiagonalPreconditioner< BlockMatrixBase< BlockMatrixType > >
{
public:
	typedef typename BlockMatrixTraits< BlockMatrixType >::BlockType LocalMatrixType ;
	typedef ProblemTraits< LocalMatrixType > GlobalProblemTraits ;
	typedef typename GlobalProblemTraits::DynVector Vector ;

	explicit DiagonalPreconditioner( const BlockMatrixBase< BlockMatrixType > &matrix )
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
			for( Index k = 0 ; k < (Index) b.rows() ; ++k, ++cur_row )
			{
				if( NumTraits< Scalar >::isSquareZero( m_diagonal[ cur_row ] ) )
				{
					m_diagonal[ cur_row ] = 1 ;
				}
			}
		}
	}

	template < bool transpose, typename ResT, typename RhsT >
	void apply( const RhsT& rhs, ResT &res ) const
	{
		res = rhs.cwiseQuotient( m_diagonal ) ;
	}

private:
	Vector m_diagonal ;
} ;

template < typename BlockMatrixType, typename FactorizationType >
class DiagonalFactorizationPreconditioner
{
public:
	typedef typename BlockMatrixTraits< BlockMatrixType >::BlockType BlockType ;
	typedef typename BlockMatrixTraits< BlockMatrixType >::Index Index ;

	template < bool transpose, typename ResT, typename RhsT >
	void apply( const RhsT& rhs, ResT &res ) const
	{
		res.setZero() ;
		m_fact.template multiply< transpose >( rhs, res )  ;
	}

protected:
	DiagonalFactorizationPreconditioner( const BlockMatrixBase< BlockMatrixType > &matrix )
	{
		m_fact.cloneDimensions( matrix ) ;

		for( Index i = 0 ; i < matrix.rowsOfBlocks() ; ++i )
		{
			m_fact.insertBack( i, i ).compute( matrix.diagonal( i ) ) ;
		}
		m_fact.finalize() ;
	}

private:

	 SparseBlockMatrix< FactorizationType > m_fact ;
} ;

template < typename BlockMatrixType >
class DiagonalLUPreconditioner< BlockMatrixBase< BlockMatrixType > >
		: public DiagonalFactorizationPreconditioner< BlockMatrixType ,
		typename MatrixTraits< typename BlockMatrixTraits< BlockMatrixType >::BlockType >::LUType >
{
public:
	typedef DiagonalFactorizationPreconditioner< BlockMatrixType,
		typename MatrixTraits< typename BlockMatrixTraits< BlockMatrixType >::BlockType >::LUType >
	Base ;

	explicit DiagonalLUPreconditioner( const BlockMatrixBase< BlockMatrixType > &matrix )
		: Base( matrix )
	{}

} ;

template < typename BlockMatrixType >
class DiagonalLDLTPreconditioner< BlockMatrixBase< BlockMatrixType > >
		: public DiagonalFactorizationPreconditioner< BlockMatrixType ,
		typename MatrixTraits< typename BlockMatrixTraits< BlockMatrixType >::BlockType >::LDLTType >
{
public:
	typedef DiagonalFactorizationPreconditioner< BlockMatrixType,
		typename MatrixTraits< typename BlockMatrixTraits< BlockMatrixType >::BlockType >::LDLTType >
	Base ;

	explicit DiagonalLDLTPreconditioner( const BlockMatrixBase< BlockMatrixType > &matrix )
		: Base( matrix )
	{}

} ;

} //naemspace bogus

#endif
