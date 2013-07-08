#ifndef BOGUS_FRICTION_PROBLEM_HPP
#define BOGUS_FRICTION_PROBLEM_HPP

#include "../Core/Block.hpp"
#include "../Core/BlockSolvers.fwd.hpp"

namespace bogus
{

template< unsigned Dimension >
struct PrimalFrictionProblem
{
	// Primal Data
	//! M -- mass matrix
	SparseBlockMatrix< Eigen::MatrixXd, COMPRESSED  > M ;
	//! E -- local rotation matrix ( world <-> contact basis )
	bogus::SparseBlockMatrix< Eigen::Matrix< double, Dimension, Dimension >, COMPRESSED > E ;

	typedef Eigen::Matrix< double, Dimension, Eigen::Dynamic > HBlock ;
	//! H -- deformation gradient ( generalized coordinates <-> 3D world )
	SparseBlockMatrix< HBlock > H;

	//! External forces
	const double *f ;
	//! Free velocity ( such that u = Hv + w )
	const double *w ;
	//! Coulomb friction coefficients
	const double *mu ;

	// Cached data

	//! M^-1
	SparseBlockMatrix< LU< Eigen::MatrixBase< Eigen::MatrixXd > >, COMPRESSED > MInv ;

} ;


template< unsigned Dimension >
struct DualFrictionProblem
{
	typedef SparseBlockMatrix< Eigen::Matrix< double, Dimension, Dimension, Eigen::RowMajor >,
							   SYMMETRIC | COMPRESSED > WType ;
	typedef GaussSeidel< WType > GaussSeidelType ;

	//! W -- Delassus operator
	WType W ;

	//! Rhs ( such that u = Wr + b )
	Eigen::VectorXd b ;

	//! Coulomb friction coefficients
	const double *mu ;

	void computeFrom( PrimalFrictionProblem< Dimension >& primal ) ;

	double solveWith( GaussSeidelType &gs, double * r, const bool staticProblem = false ) const ;

} ;

} //namespace bogus

#endif
