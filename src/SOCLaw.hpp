#ifndef BOGUS_SOCLAW_HPP
#define BOGUS_SOCLAW_HPP

#include <Eigen/Core>

namespace bogus
{

template < bool DeSaxceCOV = false >
struct SOCLaw
{
	Eigen::VectorXd mu ;

	template< typename Derived, typename OtherDerived >
	typename Derived::Scalar eval( const Eigen::MatrixBase< Derived > &x,
						  const Eigen::MatrixBase< OtherDerived > &y ) const
	{ (void) x ; (void) y ; return 0.; }

	template< typename MatDerived, typename VecDerived >
	bool solveLocal(
			const unsigned problemIndex,
			const Eigen::MatrixBase< MatDerived > &A,
			const Eigen::MatrixBase< VecDerived > &b,
			Eigen::MatrixBase< VecDerived > &x
			) const
	{ (void) problemIndex ; (void) b ; (void) A ; (void) x ; return true; }
} ;

typedef SOCLaw< true > CoulombLaw ;

}

#endif // SOCLAW_HPP
