#ifndef BOGUS_SOCLAW_HPP
#define BOGUS_SOCLAW_HPP

#include <Eigen/Core>

#include "../GaussSeidel/LocalProblem.hpp"
#include "FischerBurmeister.hpp"

#include <iostream>

namespace bogus
{

template < typename LocalMatrixType, bool DeSaxceCOV = false >
struct SOCLaw
{
	typedef LocalProblemTraits< LocalMatrixType > ProblemTraits ;
	typedef typename ProblemTraits::Scalar Scalar ;

	typename ProblemTraits::DynVector mu ;
	Scalar localTol ;

	SOCLaw() ;

	template< typename Derived, typename OtherDerived >
	Scalar eval(
			const Eigen::MatrixBase< Derived > &x,
			const Eigen::MatrixBase< OtherDerived > &y ) const
	{
		typedef FischerBurmeister< ProblemTraits::dimension, typename ProblemTraits::Scalar, DeSaxceCOV > FBFunction ;

		Scalar sum = 0. ;
		typename ProblemTraits::Vector lx, ly, fb ;

		assert( x.rows() == mu.size() * ProblemTraits::dimension ) ;
		assert( y.rows() == mu.size() * ProblemTraits::dimension ) ;

		for( unsigned i = 0 ; i < mu.size() ; ++ i )
		{
			lx = ProblemTraits::segment( i, x ) ;
			ly = ProblemTraits::segment( i, y ) ;
			FBFunction::compute( mu[i], lx, ly, fb ) ;
			sum += fb.squaredNorm() ;
		}

		return sum / ( 1 + mu.size() );
	}

	bool solveLocal(
			const unsigned problemIndex,
			const typename ProblemTraits::Matrix &A,
			const typename ProblemTraits::Vector &b,
			typename ProblemTraits::Vector &x
			) const ;

} ;

typedef SOCLaw< Eigen::Matrix2d,  true > Coulomb2D ;
typedef SOCLaw< Eigen::Matrix3d,  true > Coulomb3D ;
typedef SOCLaw< Eigen::Matrix2d, false > SOC2D ;
typedef SOCLaw< Eigen::Matrix3d, false > SOC3D ;

}

#endif // SOCLAW_HPP
