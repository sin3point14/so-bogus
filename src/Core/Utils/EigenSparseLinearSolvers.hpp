#ifndef BOGUS_EIGEN_SPARSE_LINEAR_SOLVERS
#define BOGUS_EIGEN_SPARSE_LINEAR_SOLVERS

#if EIGEN_VERSION_AT_LEAST(3,1,0)

#include "LinearSolver.hpp"

#include <Eigen/Sparse>

#include <tr1/memory>

namespace bogus {

template < typename Derived >
struct LinearSolverTraits< LDLT< Eigen::SparseMatrixBase< Derived > > >
{
  typedef typename Derived::PlainObject MatrixType ;
  typedef Eigen::SimplicialLDLT< MatrixType > FactType ;

  template < typename RhsT > struct Result {
	  typedef Eigen::internal::solve_retval< Eigen::SimplicialCholeskyBase< FactType >, RhsT > Type ;
  } ;
  template < typename RhsT >
  struct Result< Eigen::MatrixBase< RhsT > > {
	  typedef typename Result< RhsT >::Type Type ;
  } ;
} ;


template < typename Derived >
struct LDLT< Eigen::SparseMatrixBase< Derived > >
		: public LinearSolverBase< LDLT< Eigen::SparseMatrixBase< Derived > > >
{
	typedef Eigen::SparseMatrixBase< Derived > MatrixType ;
	typedef LinearSolverTraits< LDLT< MatrixType > > Traits ;

	LDLT() {}
	template< typename OtherDerived >
	explicit LDLT ( const Eigen::SparseMatrixBase< OtherDerived >& mat )
		: m_fact( new typename Traits::FactType( mat ) )
	{}

	template< typename OtherDerived >
	LDLT& compute ( const Eigen::SparseMatrixBase< OtherDerived >& mat )
	{
		m_fact.reset( new typename Traits::FactType( mat ) ) ;
		return *this ;
	}

	template < typename RhsT >
	typename Traits::template Result< Eigen::MatrixBase< RhsT > >::Type
	solve( const Eigen::MatrixBase< RhsT >& rhs ) const
	{
		assert( m_fact ) ;
		return m_fact->solve( rhs ) ;
	}

  private:
	std::tr1::shared_ptr< typename Traits::FactType > m_fact ;
} ;

template < typename Scalar, int _Options = 0, typename _Index = int >
struct SparseLDLT : public LDLT< Eigen::SparseMatrixBase< Eigen::SparseMatrix< Scalar, _Options, _Index > > >
{
	SparseLDLT() {}
	template< typename OtherDerived >
	explicit SparseLDLT ( const Eigen::SparseMatrixBase< OtherDerived >& mat )
		: LDLT< Eigen::SparseMatrixBase< Eigen::SparseMatrix< Scalar, _Options, _Index > > >( mat )
	{}
} ;

} //namespace bogus

#endif

#endif
