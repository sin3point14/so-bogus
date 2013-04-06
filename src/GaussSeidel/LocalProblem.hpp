#ifndef BOGUS_LOCAL_PROBLEM_HPP
#define BOGUS_LOCAL_PROBLEM_HPP

#include <Eigen/Core>

namespace bogus
{

template < typename MatrixType >
class LocalProblemTraits
{
	typedef MatrixType Matrix ;
	enum { dimension = Matrix::RowsAtCompileTime  } ;

	typename Eigen::internal::traits<Matrix>::Scalar Scalar ;
	typedef Eigen::Matrix< Scalar, dimension, 1 > Vector ;
} ;

template < typename MatrixType >
class LocalProblem
{
public:
	typedef LocalProblemTraits< MatrixType > Traits ;

	// y = Ax + b
	Traits::Matrix  A ;
	Traits::Scalar scaling ;
} ;

}

#endif
