#include "MecheInterface.hpp"

#include "../Core/Block.impl.hpp"
#include "../Core/GaussSeidel.impl.hpp"
#include "../Core/SecondOrder.impl.hpp"

#include <algorithm>

namespace bogus
{

struct MecheFrictionProblem::Data
{
	// Primal Data
	//! M^-1
	SparseBlockMatrix< LU< Eigen::MatrixBase< Eigen::MatrixXd > >, flags::COMPRESSED > MInv ;
	//! E
	bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::COMPRESSED > E ;
	//! H
	typedef Eigen::Matrix< double, 3, Eigen::Dynamic > HBlock ;
	SparseBlockMatrix< HBlock, flags::COL_MAJOR > H;

	const double *f ;
	const double *w ;
	const double *mu ;

	//Dual

	//! W
	typedef bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::SYMMETRIC | bogus::flags::COMPRESSED > WType ;
	WType W ;

	Eigen::VectorXd b ;

	// Cached data

	//! M^-1 * H'
	typedef Eigen::Matrix< double, Eigen::Dynamic, 3 > HtBlock ;
	SparseBlockMatrix< HtBlock > MInvHt ;

	Eigen::VectorXd MInvf ;

} ;

MecheFrictionProblem::~MecheFrictionProblem()
{
	delete m_data ;
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

	std::vector< unsigned > dofs( NObj ) ;
	std::copy( ndof, ndof + NObj, dofs.begin() ) ;

	// Build M^-1

	m_data->MInv.reserve( NObj ) ;
	m_data->MInv.setRows( dofs ) ;
	m_data->MInv.setCols( dofs ) ;

	for( unsigned i = 0 ; i < NObj ; ++i )
	{
		m_data->MInv.insertBack( i, i ) ;
	}
	m_data->MInv.finalize() ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( unsigned i = 0 ; i < NObj ; ++ i )
	{
		m_data->MInv.block(i).compute( Eigen::MatrixXd::Map( MassMat[i], dofs[i], dofs[i] ) ) ;
	}

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
	m_data->H.setCols( dofs ) ;
	for( unsigned i = 0 ; i < n_in ; ++i )
	{
		const Eigen::Matrix3d Et = m_data->E.diagonal(i).transpose() ;
		if( ObjB[i] == -1 )
		{
			m_data->H.insertBack( i, ObjA[i] ) =  Et *
					Eigen::MatrixXd::Map( HA[i], 3, dofs[ ObjA[i] ] ) ;
		} else if( ObjB[i] == ObjA[i] )
		{
			m_data->H.insertBack( i, ObjA[i] ) =  Et *
					( Eigen::MatrixXd::Map( HA[i], 3, dofs[ ObjA[i] ] ) -
					Eigen::MatrixXd::Map( HB[i], 3, dofs[ ObjA[i] ] ) ) ;
		} else {
			m_data->H.insertBack( i, ObjA[i] ) =  Et *
					Eigen::MatrixXd::Map( HA[i], 3, dofs[ ObjA[i] ] ) ;
			m_data->H.insertBack( i, ObjB[i] ) =  - Et *
					Eigen::MatrixXd::Map( HB[i], 3, dofs[ ObjB[i] ] ) ;
		}
	}
	m_data->H.finalize() ;

	m_data->f = f_in ;
	m_data->w = w_in ;
	m_data->mu = mu_in ;
}


double MecheFrictionProblem::solve(
		double *r,
		double *v,
		bool deterministic, //!< Whether the Gauss-Seidel should be eterministic
		double tol,                  //!< Gauss-Seidel tolerance. 0. means GS's default
		unsigned maxIters //!< Max number of iterations. 0 means GS's default
		)
{
	assert( m_data ) ;
	const unsigned m = m_data->H.cols() ;
	const unsigned n = m_data->H.rowsOfBlocks() ;

	// If dual has not been computed yet
	if( m_data->W.rowsOfBlocks() != n )
	{
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

	double res = gs.solve( bogus::Coulomb3D( n, m_data->mu ), m_data->b, r_loc ) ;

	// compute v
	if( v )
	{
		Eigen::VectorXd::Map( v, m ) = m_data->MInvHt * r_loc -  m_data->MInvf ;
	}

	// r to world coords
	Eigen::VectorXd::Map( r, 3*n ) = m_data->E * r_loc ;


	return res ;
}


}


