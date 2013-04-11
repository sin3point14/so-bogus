#include "MecheInterface.hpp"

#include "../Core/Block.impl.hpp"
#include "../Core/GaussSeidel.impl.hpp"
#include "../Core/SecondOrder.impl.hpp"

#include <algorithm>

namespace bogus
{

double MecheInterface::solveFrictionProblem (
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
		const double *const HB[], //!< array of size \a n, containing pointers to a dense, colum-major matrix of size <c> d*ndof[ObjA[i]] </c> corresponding to the H-matrix of <c> ObjB[i] </c> (\c NULL for an external object)
		double r0[], //!< length \a nd : initialization for \a r (in world space coordinates) + used to return computed r
		double v0[] //!< length \a m: initialization for v + to return computed v
		  )
{
	std::vector< unsigned > dofs( NObj ) ;
	std::copy( ndof, ndof + NObj, dofs.begin() ) ;

	// Build M^-1

	SparseBlockMatrix< LU< Eigen::MatrixBase< Eigen::MatrixXd > >, flags::COMPRESSED > MInv ;
	MInv.reserve( NObj ) ;
	MInv.setRows( dofs ) ;
	MInv.setCols( dofs ) ;

	for( unsigned i = 0 ; i < NObj ; ++i )
	{
	   MInv.insertBack( i, i ) ;
	}
	MInv.finalize() ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
	for( unsigned i = 0 ; i < NObj ; ++ i )
	{
		MInv.block(i).compute( Eigen::MatrixXd::Map( MassMat[i], dofs[i], dofs[i] ) ) ;
	}

	// E
	Eigen::Map< const Eigen::Matrix< double, Eigen::Dynamic, 3 > > E_flat( E_in, 3*n_in, 3 ) ;

	bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::COMPRESSED > E ;
	E.reserve( n_in ) ;
	E.setRows( n_in, 3 ) ;
	E.setCols( n_in, 3 ) ;
	for( unsigned i = 0 ; i < n_in ; ++i )
	{
		E.insertBack( i, i ) = E_flat.block< 3,3 > ( 3*i, 0 ) ;
	}
	E.finalize() ;
	E.cacheTranspose() ;


	// Build H
	typedef Eigen::Matrix< double, 3, Eigen::Dynamic > HBlock ;
	SparseBlockMatrix< HBlock, flags::COL_MAJOR > H;
	H.reserve( 2*n_in ) ;
	H.setRows( n_in, 3 ) ;
	H.setCols( dofs ) ;
	for( unsigned i = 0 ; i < n_in ; ++i )
	{
		const Eigen::Matrix3d Et = E.diagonal(i).transpose() ;
		if( ObjB[i] == -1 )
		{
			H.insertBack( i, ObjA[i] ) =  Et *
					Eigen::MatrixXd::Map( HA[i], 3, dofs[ ObjA[i] ] ) ;
		} else if( ObjB[i] == ObjA[i] )
		{
			H.insertBack( i, ObjA[i] ) =  Et *
					( Eigen::MatrixXd::Map( HA[i], 3, dofs[ ObjA[i] ] ) -
					  Eigen::MatrixXd::Map( HB[i], 3, dofs[ ObjA[i] ] ) ) ;
		} else {
			H.insertBack( i, ObjA[i] ) =  Et *
					Eigen::MatrixXd::Map( HA[i], 3, dofs[ ObjA[i] ] ) ;
			H.insertBack( i, ObjB[i] ) =  - Et *
					Eigen::MatrixXd::Map( HB[i], 3, dofs[ ObjB[i] ] ) ;
		}
	}
	H.finalize() ;

//	std::cout << H << std::endl ;

	// M^-1 * H'
	typedef Eigen::Matrix< double, Eigen::Dynamic, 3 > HtBlock ;
	const SparseBlockMatrix< HtBlock > MInvHt ( MInv * H.transpose() );

//	std::cout << MInvHt << std::endl ;

	//W
	typedef bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::SYMMETRIC | bogus::flags::COMPRESSED > WType ;
	WType W = H * MInvHt ;

	W.cacheTranspose() ; // More efficient Gauss-seidel
	bogus::GaussSeidel< WType > gs( W ) ;

	// b
	const Eigen::VectorXd MInvf = MInv * Eigen::VectorXd::Map( f_in, MInv.rows() ) ;
	const Eigen::VectorXd uf =Eigen::VectorXd::Map( w_in, 3*n_in) ;
	const Eigen::VectorXd b = ( E.transpose() * uf ) - H * ( MInvf );

//	std::cout << W << std::endl ;
//	std::cout << b.transpose() << std::endl ;
//	std::cout << ( Eigen::VectorXd::Map( w_in, 3*n_in) ).transpose() << std::endl ;
//	std::cout << ( E.transpose() * uf ).transpose() << std::endl ;
//	std::cout << ( E.transpose() * Eigen::VectorXd::Map( w_in, 3*n_in) ).transpose() << std::endl ;

	// mu
	std::vector< double > mu( n_in ) ;
	std::copy( mu_in, mu_in + n_in, mu.begin() ) ;

	// r to local coords
	Eigen::Map< Eigen::VectorXd > r( r0, 3*n_in ) ;
	Eigen::VectorXd r_loc = E.transpose() * r ;

	double res = gs.solve( bogus::Coulomb3D( mu ), b, r_loc ) ;

	// compute v
	Eigen::VectorXd::Map( v0, MInv.rows() ) = MInvHt * r_loc -  MInvf ;

	// r to world coords
	r = E * r_loc ;


	return res ;
}


}


