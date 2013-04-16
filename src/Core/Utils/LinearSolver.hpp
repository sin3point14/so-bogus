#ifndef BOGUS_LINEAR_SOLVER_HPP
#define BOGUS_LINEAR_SOLVER_HPP

#include <cassert>

namespace bogus {

template < typename LSDerived >
struct LinearSolverTraits {} ;


template < typename Derived >
struct LinearSolverBase
{
	template < typename RhsT >
	typename LinearSolverTraits< Derived >::template Result< RhsT >::Type
	solve( const RhsT& rhs ) const
	{
	   return static_cast< const Derived& >( *this ).solve( rhs ) ;
	}

} ;


template < typename MatrixType >
struct LU : public LinearSolverBase< LU< MatrixType > >
{ } ;

template < typename MatrixType >
struct LDLT : public LinearSolverBase< LDLT< MatrixType > >
{ } ;

}

#endif
