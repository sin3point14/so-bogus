#ifndef BOGUS_CADOUX_HPP
#define BOGUS_CADOUX_HPP

#include "../Core/BlockSolvers/ConstrainedSolverBase.impl.hpp"

#include "../Extra/SecondOrder.impl.hpp"

namespace bogus {

//! Solves a Coulomb friction problem using the Cadoux algorithm ( with fixed-point iteration )
/*!
	  See \cite ACLM11
	  \param W  The Delassus operator
	  \param b  The affine term of the linear system u = Wr + b
	  \param mu The vector of friction coefficients
	  \param minimizer The minimizer to use for each inner SOCQP problem
	  \param r  Both the initial guess and the result
	  \param cadouxIterations Number of fixed-point iterations
	  \param callback 0, or a pointer to a user-defined function that takes ( unsigned iteration, double residual ) as arguments
	  \returns the error as returned by the minimizer eval() function
	  */
template< unsigned Dimension, typename WType, template <typename> class Method >
static typename WType::Scalar solveCadoux(
		const WType& W,
		const typename WType::Scalar* b, const typename WType::Scalar* mu,
		ConstrainedSolverBase< Method, WType > &minimizer,
		typename WType::Scalar *r, const unsigned cadouxIterations,
		const Signal<unsigned, typename WType::Scalar> *callback )
{
	typedef typename WType::Scalar Scalar ;
	const std::ptrdiff_t n = W.rowsOfBlocks() ;

	SOCLaw< Dimension, Scalar, true  > coulombLaw ( n, mu )	;
	SOCLaw< Dimension, Scalar, false > socLaw	  ( n, mu ) ;

	minimizer.setMatrix( W );
	Eigen::Map< Eigen::VectorXd > r_map ( r, W.rows() ) ;
	Eigen::Map< const Eigen::VectorXd > b_map ( b, W.rows() ) ;

	Eigen::VectorXd s( W.rows() ) ;

	Scalar res = -1 ;
	const Scalar tol = minimizer.tol() ;
	minimizer.setTol( 1.e-1 * tol ) ;	//We might experience slow convergence is GS not precise enough

	for( unsigned cdxIter = 0 ; cdxIter < cadouxIterations ; ++cdxIter )
	{
		s = W * r_map + b_map ;

		res = minimizer.eval( coulombLaw, s, r_map ) ;

		if( callback ) callback->trigger( cdxIter, res ) ;
		if( cdxIter > 0 && res < tol ) break ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
		for( std::ptrdiff_t i = 0 ; i < n ; ++i )
		{
			s[ Dimension*i ] = s.segment< Dimension-1 >( Dimension*i+1 ).norm() * mu[i] ;
			s.segment< Dimension-1  >( Dimension*i+1 ).setZero() ;
		}

		s += b_map ;

		minimizer.solve( socLaw, s, r_map ) ;

	}

	minimizer.setTol( tol ) ;

	return res ;
}


} //bogus

#endif
