/*
 * This file is part of So-bogus, a C++ sparse block matrix library and
 * Second Order Cone solver.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * So-bogus is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.

 * So-bogus is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with So-bogus.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "MecheInterface.hpp"

#include "FrictionProblem.hpp"

#include "../Core/Block.impl.hpp"
#include "../Core/Block.io.hpp"
#include "../Core/BlockSolvers/GaussSeidel.impl.hpp"

#include <algorithm>

#ifdef BOGUS_WITH_BOOST_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <fstream>
#endif

namespace bogus
{

MecheFrictionProblem::MecheFrictionProblem()
	: m_primal( 0 ), m_dual( 0 ),
	  m_f( 0 ), m_w( 0 ), m_mu( 0 ),
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
	delete m_primal ;
	m_primal = 0 ;
	delete m_dual ;
	m_dual = 0 ;
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

void MecheFrictionProblem::reset ()
{
	destroy() ;

	m_primal = new PrimalFrictionProblem<3u>() ;
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
	reset() ;

	// Copy M
	// We don't actually need it after having computed a factorization of M, but we keep it around
	// in case we want to use dumpToFile()

	m_primal->M.reserve( NObj ) ;
	m_primal->M.setRows( NObj, ndof ) ;
	m_primal->M.setCols( NObj, ndof ) ;

	for( unsigned i = 0 ; i < NObj ; ++i )
	{
		m_primal->M.insertBack( i, i ) = Eigen::MatrixXd::Map( MassMat[i], ndof[i], ndof[i] ) ;
	}
	m_primal->M.finalize() ;

	// E
	Eigen::Map< const Eigen::Matrix< double, Eigen::Dynamic, 3 > > E_flat( E_in, 3*n_in, 3 ) ;

	m_primal->E.reserve( n_in ) ;
	m_primal->E.setRows( n_in, 3 ) ;
	m_primal->E.setCols( n_in, 3 ) ;
	for( unsigned i = 0 ; i < n_in ; ++i )
	{
		m_primal->E.insertBack( i, i ) = E_flat.block< 3,3 > ( 3*i, 0 ) ;
	}
	m_primal->E.finalize() ;
	m_primal->E.cacheTranspose() ;


	// Build H
	m_primal->H.reserve( 2*n_in ) ;
	m_primal->H.setRows( n_in, 3 ) ;
	m_primal->H.setCols( NObj, ndof ) ;
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( std::ptrdiff_t i = 0 ; i < (std::ptrdiff_t) n_in ; ++i )
	{
		const Eigen::Matrix3d Et = m_primal->E.diagonal(i).transpose() ;
		if( ObjB[i] == -1 )
		{
			m_primal->H.insertBack( i, ObjA[i] ) =  Et *
					Eigen::MatrixXd::Map( HA[i], 3, ndof[ ObjA[i] ] ) ;
		} else if( ObjB[i] == ObjA[i] )
		{
			m_primal->H.insertBack( i, ObjA[i] ) =  Et *
					( Eigen::MatrixXd::Map( HA[i], 3, ndof[ ObjA[i] ] ) -
					Eigen::MatrixXd::Map( HB[i], 3, ndof[ ObjA[i] ] ) ) ;
		} else {
			m_primal->H.insertBack( i, ObjA[i] ) =  Et *
					Eigen::MatrixXd::Map( HA[i], 3, ndof[ ObjA[i] ] ) ;
			m_primal->H.insertBack( i, ObjB[i] ) =  - Et *
					Eigen::MatrixXd::Map( HB[i], 3, ndof[ ObjB[i] ] ) ;
		}
	}
	m_primal->H.finalize() ;

	m_primal->f = f_in ;
	m_primal->w = w_in ;
	m_primal->mu = mu_in ;
}

unsigned MecheFrictionProblem::nDegreesOfFreedom() const
{
	return m_primal ? m_primal->M.rows() : 0u ;
}

unsigned MecheFrictionProblem::nContacts() const
{
	return m_primal ? m_primal->H.rowsOfBlocks() : 0u ;
}

void MecheFrictionProblem::computeDual( double regularization )
{
	delete m_dual ;
	m_dual = new DualFrictionProblem<3u>() ;
	m_dual->computeFrom( *m_primal );

	if( regularization > 0. )
	{
		for( int i = 0 ; i < m_dual->W.rowsOfBlocks() ; ++ i )
		{
			m_dual->W.diagonal( i ).diagonal() += Eigen::Vector3d::Constant( regularization ) ;
		}
	}
}


double MecheFrictionProblem::solve(double *r,
		double *v,
		int maxThreads, //!< Whether the Gauss-Seidel should be eterministic
		double tol,                  //!< Gauss-Seidel tolerance. 0. means GS's default
		unsigned maxIters, //!< Max number of iterations. 0 means GS's default
		bool staticProblem,
		double regularization,
		bool useInfinityNorm
								   )
{
	assert( m_primal ) ;
	const unsigned m = m_primal->H.cols() ;
	const unsigned n = m_primal->H.rowsOfBlocks() ;

	// If dual has not been computed yet
	if( !m_dual )
	{
		computeDual( staticProblem ? regularization : 0. );
	}


	// r to local coords
	Eigen::VectorXd r_loc = m_primal->E.transpose() * Eigen::VectorXd::Map( r, 3*n ) ;

	bogus::DualFrictionProblem<3u>::GaussSeidelType gs ;
	if( tol != 0. ) gs.setTol( tol );
	if( maxIters != 0 ) gs.setMaxIters( maxIters );
	gs.enableColoring( maxThreads > 1 );
	gs.setMaxThreads( maxThreads );

	gs.callback().connect( *this, &MecheFrictionProblem::ackCurrentResidual );
	gs.setAutoRegularization( regularization ) ;
	gs.useInfinityNorm( useInfinityNorm ) ;

	const double res = m_dual->solveWith( gs, r_loc.data(), staticProblem ) ;

	// compute v
	if( v )
	{
		Eigen::VectorXd::Map( v, m ) = m_primal->MInv * (
					m_primal->H.transpose() * r_loc -
					Eigen::VectorXd::Map( m_primal->f, m_primal->H.cols() ) ) ;
	}

	if( m_out && n != 0 )
	{
		*m_out << "Max coeff : " << r_loc.lpNorm< Eigen::Infinity >() << std::endl ;
	}

	// r to world coords
	Eigen::VectorXd::Map( r, 3*n ) = m_primal->E * r_loc ;


	return res ;
}

void MecheFrictionProblem::setOutStream( std::ostream *out )
{
	m_out = out ;
}

#ifdef BOGUS_WITH_BOOST_SERIALIZATION
bool MecheFrictionProblem::dumpToFile( const char* fileName, const double * r0 ) const
{
	if( !m_primal ) return false ;

	std::ofstream ofs( fileName );
	boost::archive::binary_oarchive oa(ofs);
	oa << m_primal->M << m_primal->H << m_primal->E ;
	oa << boost::serialization::make_array( m_primal->f , nDegreesOfFreedom() ) ;
	oa << boost::serialization::make_array( m_primal->w , 3 * nContacts() ) ;
	oa << boost::serialization::make_array( m_primal->mu, nContacts() ) ;
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

	reset() ;

	boost::archive::binary_iarchive ia(ifs);
	ia >> m_primal->M >> m_primal->H >> m_primal->E ;

	m_f  = new double[ nDegreesOfFreedom() ] ;
	m_w  = new double[ 3 * nContacts() ] ;
	m_mu = new double[ nContacts() ] ;

	std::cout << fileName << ": " << nDegreesOfFreedom() << " dofs, " << nContacts() << " contacts" << std::endl ;

	ia >> boost::serialization::make_array( m_f , nDegreesOfFreedom() ) ;
	ia >> boost::serialization::make_array( m_w , 3 * nContacts() ) ;
	ia >> boost::serialization::make_array( m_mu, nContacts() ) ;

	m_primal->f  = m_f ;
	m_primal->w  = m_w ;
	m_primal->mu = m_mu ;

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


