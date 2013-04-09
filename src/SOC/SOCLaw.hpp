#ifndef BOGUS_SOCLAW_HPP
#define BOGUS_SOCLAW_HPP

#include <Eigen/Core>

#include "../GaussSeidel/LocalProblem.hpp"
#include "FischerBurmeister.hpp"
#include "LocalSOCSolver.hpp"

#include <vector>

namespace bogus
{

template < typename LocalMatrixType, bool DeSaxceCOV,
		   local_soc_solver::Strategy Strat = local_soc_solver::Hybrid  >
class SOCLaw
{
public:
	typedef LocalProblemTraits< LocalMatrixType > ProblemTraits ;
	typedef typename ProblemTraits::Scalar Scalar ;

	SOCLaw( const std::vector< double >& mu ) ;

	template< typename Derived, typename OtherDerived >
	Scalar eval(
			const Eigen::MatrixBase< Derived > &x,
			const Eigen::MatrixBase< OtherDerived > &y ) const
	{
		typedef FischerBurmeister< ProblemTraits::dimension, typename ProblemTraits::Scalar, DeSaxceCOV > FBFunction ;

		Scalar sum = 0. ;
		typename ProblemTraits::Vector lx, ly, fb ;

		assert( (unsigned) x.rows() == m_mu.size() * ProblemTraits::dimension ) ;
		assert( (unsigned) y.rows() == m_mu.size() * ProblemTraits::dimension ) ;

		for( unsigned i = 0 ; i < m_mu.size() ; ++ i )
		{
			lx = ProblemTraits::segment( i, x ) ;
			ly = ProblemTraits::segment( i, y ) ;
			FBFunction::compute( m_mu[i], lx, ly, fb ) ;
			sum += fb.squaredNorm() ;
		}

		return sum / ( 1 + m_mu.size() );
	}

	bool solveLocal(
			const unsigned problemIndex,
			const typename ProblemTraits::Matrix &A,
			const typename ProblemTraits::Vector &b,
			typename ProblemTraits::Vector &xm,
			const Scalar scaling
			) const ;

private:

	const std::vector< double > &m_mu ;
	Scalar m_localTol ;

} ;

typedef SOCLaw< Eigen::Matrix2d,  true > Coulomb2D ;
typedef SOCLaw< Eigen::Matrix3d,  true > Coulomb3D ;
typedef SOCLaw< Eigen::Matrix2d, false > SOC2D ;
typedef SOCLaw< Eigen::Matrix3d, false > SOC3D ;

}

#endif // SOCLAW_HPP
