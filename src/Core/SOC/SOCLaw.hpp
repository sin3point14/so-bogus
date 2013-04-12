#ifndef BOGUS_SOCLAW_HPP
#define BOGUS_SOCLAW_HPP

#include <Eigen/Core>

#include "../GaussSeidel.fwd.hpp"
#include "FischerBurmeister.hpp"
#include "LocalSOCSolver.hpp"

#include <vector>

namespace bogus
{

template < typename LocalMatrixType, bool DeSaxceCOV,
		   local_soc_solver::Strategy Strat >
class SOCLaw
{
public:
	typedef LocalProblemTraits< LocalMatrixType > ProblemTraits ;
	typedef typename ProblemTraits::Scalar Scalar ;

	SOCLaw( const unsigned n, const double * mu ) ;

	template< typename Derived, typename OtherDerived >
	Scalar eval(
			const Eigen::MatrixBase< Derived > &x,
			const Eigen::MatrixBase< OtherDerived > &y ) const
	{
		typedef FischerBurmeister< ProblemTraits::dimension, typename ProblemTraits::Scalar, DeSaxceCOV > FBFunction ;

		assert( (unsigned) x.rows() == m_n * ProblemTraits::dimension ) ;
		assert( (unsigned) y.rows() == m_n * ProblemTraits::dimension ) ;

		Scalar sum = 0. ;
		typename ProblemTraits::Vector lx, ly, fb ;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for private( lx, ly, fb ) reduction ( + : sum )
#endif
		for( unsigned i = 0 ; i < m_n ; ++ i )
		{
			lx = ProblemTraits::segment( i, x ) ;
			ly = ProblemTraits::segment( i, y ) ;
			FBFunction::compute( m_mu[i], lx, ly, fb ) ;
			sum += fb.squaredNorm() ;
		}

		return sum / ( 1 + m_n );
	}

	bool solveLocal(
			const unsigned problemIndex,
			const typename ProblemTraits::Matrix &A,
			const typename ProblemTraits::Vector &b,
			typename ProblemTraits::Vector &xm,
			const Scalar scaling
			) const ;

private:

	const double * m_mu ;
	const unsigned m_n ;
	Scalar m_localTol ;

} ;

typedef SOCLaw< Eigen::Matrix2d,  true > Coulomb2D ;
typedef SOCLaw< Eigen::Matrix3d,  true > Coulomb3D ;
typedef SOCLaw< Eigen::Matrix2d, false > SOC2D ;
typedef SOCLaw< Eigen::Matrix3d, false > SOC3D ;

}

#endif // SOCLAW_HPP
