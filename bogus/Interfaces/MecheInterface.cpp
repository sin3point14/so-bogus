/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "MecheInterface.hpp"

#include "../Core/Block.impl.hpp"
#include "../Core/Block.io.hpp"
#include "../Core/BlockSolvers.impl.hpp"
#include "../Core/SecondOrder.impl.hpp"

#include <algorithm>

#ifdef BOGUS_WITH_BOOST_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <fstream>
#endif

namespace bogus
{

struct MecheFrictionProblem::Data
{
	// Primal Data
	//! M^-1
	SparseBlockMatrix< Eigen::MatrixXd, flags::COMPRESSED  > M ;
	//! E
	bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::COMPRESSED > E ;
	//! H
	typedef Eigen::Matrix< double, 3, Eigen::Dynamic > HBlock ;
	SparseBlockMatrix< HBlock > H;

	const double *f ;
	const double *w ;
	const double *mu ;

	//Dual

	//! W
	typedef bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::SYMMETRIC | bogus::flags::COMPRESSED > WType ;
	WType W ;

	Eigen::VectorXd b ;

	// Cached data

	//! M^-1
	SparseBlockMatrix< LU< Eigen::MatrixBase< Eigen::MatrixXd > >, flags::COMPRESSED > MInv ;

	//! M^-1 * H'
	typedef Eigen::Matrix< double, Eigen::Dynamic, 3 > HtBlock ;
	SparseBlockMatrix< HtBlock, bogus::flags::COL_MAJOR > MInvHt ;

	Eigen::VectorXd MInvf ;
} ;

MecheFrictionProblem::MecheFrictionProblem()
	: m_data( 0 ), m_f( 0 ), m_w( 0 ), m_mu( 0 ),
      m_out( &std::cout )
{
}


MecheFrictionProblem::~MecheFrictionProblem()
{
	destroy() ;
}

void MecheFrictionProblem::destroy()
{
	delete[] m_f ;
	m_f = 0 ;
	delete[] m_w ;
	m_w = 0 ;
	delete[] m_mu ;
	m_mu = 0 ;
	delete m_data ;
	m_data = 0 ;
}

void MecheFrictionProblem::ackCurrentResidual( unsigned GSIter, double err )
{
	if( m_out )
	{
		*m_out << "Finished iteration " << GSIter
		       << " with residual " << err
		       << std::endl ;
	}
}

void MecheFrictionProblem::fromPrimal (
		unsigned int NObj, //!< number of subsystems
		const unsigned int * ndof, //!< array of size \a NObj, the number of degree of freedom of each subsystem
		const double *const * MassMat, //!< array of pointers to the mass matrix of each subsystem
		const double * f_in, //!< the constant term in \f$ M v + f= {}^t \! H r \f$
		unsigned int n_in, //!< number of contact points
		const double * mu_in, //!< array of size \a n giving the friction coeffs
		const double * E_in, //!< array of size \f$ n \times d \times d \f$ giving the \a n normals followed by the \a n tangent vectors (and by again \a n tangent vectors if \a d is 3). Said otherwise, \a E is a \f$ (nd) \times d \f$ matrix, stored column-major, formed by \a n blocks of size \f$ d \times d \f$ with each block being an orthogonal matrix (the transition matrix from the world space coordinates \f$ (x_1, x_2, x_3) \f$ to the local coordinates \f$ (x_N, x_{T1}, x_{T2}) \f$
		const double * w_in, //!< array of size \a nd, the constant term in \f$ u = H v + w \f$
		const int * const ObjA, //!< array of size \a n, the first object involved in the \a i-th contact (must be an internal object) (counted from 0)
		const int * const ObjB, //!< array of size \a n, the second object involved in the \a i-th contact (-1 for an external object) (counted from 0)
		const double *const HA[], //!< array of size \a n, containing pointers to a dense, colum-major matrix of size <c> d*ndof[ObjA[i]] </c> corresponding to the H-matrix of <c> ObjA[i] </c>
		const double *const HB[] //!< array of size \a n, containing pointers to a dense, colum-major matrix of size <c> d*ndof[ObjA[i]] </c> corresponding to the H-matrix of <c> ObjB[i] </c> (\c NULL for an external object)
		)
{
	delete m_data ;
	m_data = new Data() ;

	// Build M^-1

	m_data->M.reserve( NObj ) ;
	m_data->M.setRows( NObj, ndof ) ;
	m_data->M.setCols( NObj, ndof ) ;

	for( unsigned i = 0 ; i < NObj ; ++i )
	{
		m_data->M.insertBack( i, i ) = Eigen::MatrixXd::Map( MassMat[i], ndof[i], ndof[i] ) ;
	}
	m_data->M.finalize() ;

	// E
	Eigen::Map< const Eigen::Matrix< double, Eigen::Dynamic, 3 > > E_flat( E_in, 3*n_in, 3 ) ;

	m_data->E.reserve( n_in ) ;
	m_data->E.setRows( n_in, 3 ) ;
	m_data->E.setCols( n_in, 3 ) ;
	for( unsigned i = 0 ; i < n_in ; ++i )
	{
		m_data->E.insertBack( i, i ) = E_flat.block< 3,3 > ( 3*i, 0 ) ;
	}
	m_data->E.finalize() ;
	m_data->E.cacheTranspose() ;


	// Build H
	m_data->H.reserve( 2*n_in ) ;
	m_data->H.setRows( n_in, 3 ) ;
	m_data->H.setCols( NObj, ndof ) ;
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( int i = 0 ; i < (int) n_in ; ++i )
	{
		const Eigen::Matrix3d Et = m_data->E.diagonal(i).transpose() ;
		if( ObjB[i] == -1 )
		{
			m_data->H.insertBack( i, ObjA[i] ) =  Et *
					Eigen::MatrixXd::Map( HA[i], 3, ndof[ ObjA[i] ] ) ;
		} else if( ObjB[i] == ObjA[i] )
		{
			m_data->H.insertBack( i, ObjA[i] ) =  Et *
					( Eigen::MatrixXd::Map( HA[i], 3, ndof[ ObjA[i] ] ) -
					Eigen::MatrixXd::Map( HB[i], 3, ndof[ ObjA[i] ] ) ) ;
		} else {
			m_data->H.insertBack( i, ObjA[i] ) =  Et *
					Eigen::MatrixXd::Map( HA[i], 3, ndof[ ObjA[i] ] ) ;
			m_data->H.insertBack( i, ObjB[i] ) =  - Et *
					Eigen::MatrixXd::Map( HB[i], 3, ndof[ ObjB[i] ] ) ;
		}
	}
	m_data->H.finalize() ;

	m_data->f = f_in ;
	m_data->w = w_in ;
	m_data->mu = mu_in ;
}

unsigned MecheFrictionProblem::nDegreesOfFreedom() const
{
	return m_data ? m_data->M.rows() : 0u ;
}

unsigned MecheFrictionProblem::nContacts() const
{
	return m_data ? m_data->H.rowsOfBlocks() : 0u ;
}



double MecheFrictionProblem::solve(double *r,
		double *v,
		bool deterministic, //!< Whether the Gauss-Seidel should be eterministic
		double tol,                  //!< Gauss-Seidel tolerance. 0. means GS's default
		unsigned maxIters, //!< Max number of iterations. 0 means GS's default
		bool staticProblem)
{
	assert( m_data ) ;
	const unsigned m = m_data->H.cols() ;
	const unsigned n = m_data->H.rowsOfBlocks() ;

	// If dual has not been computed yet
	if( m_data->W.rowsOfBlocks() != n )
	{
		// M^-1
		m_data->MInv.cloneStructure( m_data->M ) ;
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
		for( int i = 0 ; i < (int) m_data->M.nBlocks()  ; ++ i )
		{
			m_data->MInv.block(i).compute( m_data->M.block(i) ) ;
		}

		// M^-1 * H'
		m_data->MInvHt = m_data->MInv * m_data->H.transpose() ;

		//W
		m_data->W = m_data->H * m_data->MInvHt ;
		m_data->W.cacheTranspose() ;

		// M^-1 f, b
		m_data->MInvf = m_data->MInv * Eigen::VectorXd::Map( m_data->f, m ) ;
		m_data->b = ( m_data->E.transpose() * Eigen::VectorXd::Map( m_data->w, 3*n) )
				- m_data->H * ( m_data->MInvf );
	}


	// r to local coords
	Eigen::VectorXd r_loc = m_data->E.transpose() * Eigen::VectorXd::Map( r, 3*n ) ; ;

	bogus::GaussSeidel< Data::WType > gs( m_data->W ) ;
	if( tol != 0. ) gs.setTol( tol );
	if( maxIters != 0 ) gs.setMaxIters( maxIters );
	gs.setDeterministic( deterministic );

	gs.callback().connect( *this, &MecheFrictionProblem::ackCurrentResidual );

	double res = staticProblem
			? gs.solve( bogus::SOC3D    ( n, m_data->mu ), m_data->b, r_loc )
			: gs.solve( bogus::Coulomb3D( n, m_data->mu ), m_data->b, r_loc ) ;


	// compute v
	if( v )
	{
		Eigen::VectorXd::Map( v, m ) = m_data->MInvHt * r_loc -  m_data->MInvf ;
	}

	// r to world coords
	Eigen::VectorXd::Map( r, 3*n ) = m_data->E * r_loc ;


	return res ;
}

#ifdef BOGUS_WITH_BOOST_SERIALIZATION
bool MecheFrictionProblem::dumpToFile( const char* fileName, const double * r0 ) const
{
	if( !m_data ) return false ;

	std::ofstream ofs( fileName );
	boost::archive::binary_oarchive oa(ofs);
	oa << m_data->M << m_data->H << m_data->E ;
	oa << boost::serialization::make_array( m_data->f , nDegreesOfFreedom() ) ;
	oa << boost::serialization::make_array( m_data->w , 3 * nContacts() ) ;
	oa << boost::serialization::make_array( m_data->mu, nContacts() ) ;
	bool has_r0 = r0 != 0 ;
	oa << has_r0 ; ;
	if( r0 )
	{
		oa << boost::serialization::make_array( r0, 3*nContacts() ) ;
	}

	return true ;
}

bool MecheFrictionProblem::fromFile( const char* fileName, double *& r0 )
{
	std::ifstream ifs( fileName );
	if( !ifs.is_open() ) return false ;

	destroy() ;
	m_data = new Data() ;

	boost::archive::binary_iarchive ia(ifs);
	ia >> m_data->M >> m_data->H >> m_data->E ;

	m_f  = new double[ nDegreesOfFreedom() ] ;
	m_w  = new double[ 3 * nContacts() ] ;
	m_mu = new double[ nContacts() ] ;

	std::cout << fileName << ": " << nDegreesOfFreedom() << " dofs, " << nContacts() << " contacts" << std::endl ;

	ia >> boost::serialization::make_array( m_f , nDegreesOfFreedom() ) ;
	ia >> boost::serialization::make_array( m_w , 3 * nContacts() ) ;
	ia >> boost::serialization::make_array( m_mu, nContacts() ) ;

	m_data->f  = m_f ;
	m_data->w  = m_w ;
	m_data->mu = m_mu ;

	r0 = new double[ 3 * nContacts() ] ;

	bool has_r0 ;
	ia >> has_r0 ;
	if ( has_r0 ) {
		ia >> boost::serialization::make_array( r0, 3*nContacts() ) ;
	} else {
		Eigen::VectorXd::Map( r0, 3*nContacts() ).setZero() ;
	}

	return true ;
}

#else
bool MecheFrictionProblem::dumpToFile( const char*, const double* ) const
{
	std::cerr << "MecheInterface::dumpToFile: Error, bogus compiled without serialization capabilities" ;
	return false ;
}
bool MecheFrictionProblem::fromFile(const char*, double *& ) {
	std::cerr << "MecheInterface::fromFile: Error, bogus compiled without serialization capabilities" ;
	return false ;
}
#endif

}


