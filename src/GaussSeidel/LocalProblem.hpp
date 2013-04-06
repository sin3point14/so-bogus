#ifndef BOGUS_LOCAL_PROBLEM_HPP
#define BOGUS_LOCAL_PROBLEM_HPP

#include <Eigen/Core>

namespace bogus
{

template < typename MatrixType >
struct LocalProblemTraits
{
	typedef MatrixType Matrix ;
	enum { dimension = Matrix::RowsAtCompileTime  } ;

	typedef typename Eigen::internal::traits<Matrix>::Scalar Scalar ;
	typedef Eigen::Matrix< Scalar, dimension, 1 > Vector ;
	typedef Eigen::Matrix< Scalar, Eigen::Dynamic, 1 > DynVector ;

	template< typename VectorType >
	static typename VectorType::template FixedSegmentReturnType< dimension >::Type segment( const unsigned i, VectorType& v )
	{ return v.template segment< dimension > ( i * dimension ) ; }
	template< typename VectorType >
	static typename VectorType::template ConstFixedSegmentReturnType< dimension >::Type segment( const unsigned i, const VectorType& v )
	{ return v.template segment< dimension > ( i * dimension ) ; }
} ;

} //namespace bogus

#endif
