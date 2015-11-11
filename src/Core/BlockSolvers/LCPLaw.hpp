#ifndef BOGUS_LCPLAW_HPP
#define BOGUS_LCPLAW_HPP

#include <cmath>

namespace bogus
{

//! LCP local solver that can be used within GaussSeidel and ProjectedGradient solvers
/*!
  For demonstration purposes. Since blocks are 1x1, other libraries are probably more suited.
	\tparam Scalar the scalar type
	\tparam Dimension the dimension of the blocks of the global matrix
  */
template < typename Scalar >
class LCPLaw
{
public:
	enum{ dimension = 1 } ;

	typedef LocalProblemTraits< dimension, Scalar > Traits ;

	//! Constructor
	LCPLaw( ) {}

	//! \return \f$ \vert fb( x, y ) \vert^2_2 \f$, where fb is the scalar Fischer-Burmeister function
	Scalar eval( const unsigned problemIndex,
				 const typename Traits::Vector &x,
				 const typename Traits::Vector &y ) const ;

	//! Solves the local problem
	/*!
		0 \leq y \perp a x + b \geq 0
	*/
	bool solveLocal(
			const unsigned problemIndex,
			const typename Traits::Matrix &A,
			const typename Traits::Vector &b,
			typename Traits::Vector &x,
			const Scalar scaling
			) const ;

	//! Projects x on \f$ R^+ \f$
	void projectOnConstraint( const unsigned problemIndex, typename Traits::Vector &x ) const ;

} ;


}

#endif
