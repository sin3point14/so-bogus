#ifndef BOGUS_KRYLOV_METHODS_HPP
#define BOGUS_KRYLOV_METHODS_HPP

#include "../Utils/LinearSolverBase.hpp"
#include "../Utils/Signal.hpp"
#include "../Utils/NumTraits.hpp"
#include "../Block/Access.hpp"
#include "../BlockSolvers.fwd.hpp"

#define BOGUS_KRYLOV_METHODS \
	BOGUS_PROCESS_KRYLOV_METHOD(CG      )\
	BOGUS_PROCESS_KRYLOV_METHOD(BiCG    )\
	BOGUS_PROCESS_KRYLOV_METHOD(BiCGSTAB)\
	BOGUS_PROCESS_KRYLOV_METHOD(CGS     )\
	BOGUS_PROCESS_KRYLOV_METHOD(GMRES   )\
	BOGUS_PROCESS_KRYLOV_METHOD(TFQMR   )


namespace bogus {

template< template < typename, typename, typename > class Method,
		  typename Matrix, typename Preconditioner, typename Traits >
struct KrylovSolverBase
		: public LinearSolverBase< Method< Matrix, Preconditioner, Traits > >
{
	typedef Method< Matrix, Preconditioner, Traits > Derived ;
	typedef LinearSolverBase< Derived > Base ;
	typedef typename Traits::Scalar Scalar ;

	const Matrix &m_A ;
	const Preconditioner &m_P;
	const Signal< unsigned, Scalar >  &m_callback ;

	Scalar m_tol ;
	unsigned m_maxIters;

	KrylovSolverBase( const Matrix &A,
					  const Preconditioner &P,
					  const Signal< unsigned, Scalar > &callback,
					  Scalar tol, unsigned maxIters )
		: m_A( A ), m_P( P ), m_callback( callback ),
		  m_tol( tol ), m_maxIters( maxIters ),
		  m_scale( 1. / ( 1 + A.size() ) )
	{}

	//! Returns the solution \b x of the linear system \b M \c * \b x \c = \c rhs
	template < typename RhsT >
	typename LinearSolverTraits< Derived >::template Result< RhsT >::Type
	solve( const RhsT& rhs ) const
	{
		typename LinearSolverTraits< Derived >::template Result< RhsT >::Type
				x( m_A.rows(), rhs.cols() ) ;
		Base::solve( rhs, x ) ;
		return x ;
	}


protected:
	Scalar m_scale ;
} ;

// Forward decl + create traits

#define BOGUS_PROCESS_KRYLOV_METHOD( MethodName ) \
	namespace krylov {\
		template< typename Matrix, typename Preconditioner, class Traits > \
		struct MethodName ;\
	} \
	template< typename Matrix, typename Preconditioner, class Traits > \
	struct LinearSolverTraits< krylov::MethodName< Matrix, Preconditioner, Traits > > \
	{ \
	  typedef Matrix MatrixType ; \
	  template < typename RhsT > struct Result { \
		  typedef typename Traits::template MutableClone< RhsT >::Type Type ; \
	  } ; \
	} ;

BOGUS_KRYLOV_METHODS
#undef BOGUS_PROCESS_KRYLOV_METHOD


#define BOGUS_MAKE_KRYLOV_SOLVER_TYPEDEFS( MethodName ) \
	typedef KrylovSolverBase< krylov::MethodName, Matrix, Preconditioner, Traits > Base ; \
	typedef typename Traits::Scalar Scalar ; \
\
	using Base::m_A ;			\
	using Base::m_P ;			\
	using Base::m_maxIters ;		\
	using Base::m_tol ;			\
	using Base::m_scale ;

#define BOGUS_MAKE_KRYLOV_SOLVER_HEADER( MethodName ) \
	BOGUS_MAKE_KRYLOV_SOLVER_TYPEDEFS( MethodName )\
	MethodName( const Matrix &A, \
	const Preconditioner &P, \
	const Signal< unsigned, Scalar > &callback, \
	Scalar tol, unsigned maxIters ) \
	: Base( A, P, callback, tol, maxIters ) \
{}  \


namespace krylov {

  //! Check init guess, reset it to zero if that would give a lower residual
template < typename Matrix, typename Vector, typename RhsT, typename ResT >
typename MatrixTraits< Vector >::Scalar init( const Matrix& A, const RhsT &b, ResT &x, Vector &r0  )
{
	typedef typename MatrixTraits< Vector >::Scalar Scalar ;

	r0 = b ;
	mv_add< false >( A, x, r0, -1 ) ;

	Scalar res = r0.squaredNorm() ;
	const Scalar resAt0 = b.squaredNorm()  ;

	if( res > resAt0 ) {
		r0 = b;
		x.setZero() ;
		res = resAt0 ;
	}
	return res ;
}


// Conjugate Gradient

//! Solves ( m_A * \p x = \p b ) using the Conjugate Gradient algorithm
/*! For symmetric matrices only. Converges for positive definite linear systems.

		<b>Matrix-vector mults/iter: </b> 1
		<b>Preconditionner calls/iter: </b> 1
		<b>Storage requirements: </b> 4n
	*/
template < typename Matrix, typename Preconditioner, typename Traits >
struct CG : public KrylovSolverBase< CG, Matrix, Preconditioner, Traits >
{


	BOGUS_MAKE_KRYLOV_SOLVER_HEADER( CG )

	template < typename RhsT, typename ResT >
	Scalar solve( const RhsT &b, ResT &x ) const
	{
		typedef typename Traits::template MutableClone< RhsT >::Type Vector ;
		Vector r ;

		Scalar res = init( m_A, b, x, r ) * m_scale ;
		if( res < m_tol ) return res ;

		Vector z( r.rows() ) ;
		m_P.template apply< false >( r, z ) ;
		Vector p = z;

		Scalar zr0 = r.dot( z ) ;
		Scalar zr1 ;

		Vector Mp( m_A.rows() ) ;

		for( unsigned k = 0 ; k < m_maxIters ; ++k )
		{
			mv_set< false >(m_A, p, Mp ) ;
			const Scalar alpha = zr0 / ( p.dot( Mp ) ) ;
			x += alpha * p ;
			r -= alpha * Mp ;

			res = r.squaredNorm() * m_scale ;
			this->m_callback.trigger( k, res ) ;
			if( res < m_tol ) break ;

			m_P.template apply< false >( r, z ) ;
			zr1 = z.dot( r ) ;

			p = z + ( zr1 / zr0 ) * p ;

			zr0 = zr1 ;
		}

		return res ;
	}
} ;

// BiConjugate Gradient
//! Solves ( m_A * \p x = \p b ) using the BiConjugate Gradient algorithm
/*! Works for non-symmetric linear systems. Convergence not guaranteed.
		Requires ability to perform transpose multiplication and preconditioning

		<b>Matrix-vector mults/iter: </b> 2 ( inc. 1 transpose )
		<b>Preconditionner calls/iter: </b> 2 ( inc. 1 transpose )
		<b>Storage requirements: </b> 8n
	*/
template < typename Matrix, typename Preconditioner, typename Traits>
struct BiCG : public KrylovSolverBase< BiCG, Matrix, Preconditioner, Traits>
{

	BOGUS_MAKE_KRYLOV_SOLVER_HEADER( BiCG )

	template < typename RhsT, typename ResT >
	Scalar solve( const RhsT &b, ResT &x ) const
	{
		typedef typename Traits::template MutableClone< RhsT >::Type Vector ;
		Vector r ;

		Scalar res = init( m_A, b, x, r ) * m_scale ;
		if( res < m_tol ) return res ;

		Vector x_ = x ;

		Vector r_ = b ;
		mv_add< true >( m_A, x_, r_, 1 );

		Vector p( r.rows() ), p_ ( r_.rows() ) ;
		m_P.template apply< false >( r,  p ) ;
		m_P.template apply< true >( r_, p_ ) ;

		Vector z = p, z_ = p_ ;
		Scalar zr0 = r_.dot( z ) ;
		Scalar zr1 ;

		Vector Mp( m_A.rows() ) ;

		for( unsigned k = 0 ; k < m_maxIters ; ++k )
		{
			mv_set< false >( m_A, p, Mp ) ;
			const Scalar alpha = zr0 / ( p_.dot( Mp ) ) ;
			x  += alpha * p  ;
			x_ += alpha * p_ ;

			r  -= alpha * Mp ;
			mv_add< true >( m_A, p_, r_, -alpha ) ;

			res = r.squaredNorm() * m_scale ;
			this->m_callback.trigger( k, res ) ;
			if( res < m_tol ) break ;

			m_P.template apply< false >( r, z ) ;
			zr1 = z.dot( r_ ) ;

			const Scalar beta = ( zr1 / zr0 ) ;

			p = z + beta * p ;
			m_P.template apply< true >( r_, z_ ) ;
			p_ = z_ + beta * p_ ;

			zr0 = zr1 ;
		}


		return res ;
	}
} ;

// BiCG STAB
//! Solves ( m_A * \p x = \p b ) using the BiConjugate Gradient stabilized algorithm
/*! Works for non-symmetric linear systems. Convergence not guaranteed.
		Supposedly less erratic convergence than BiCG method.

		<b>Matrix-vector mults/iter: </b> 2
		<b>Preconditionner calls/iter: </b> 2
		<b>Storage requirements: </b> 8n
	*/
template < typename Matrix, typename Preconditioner, typename Traits>
struct BiCGSTAB : public KrylovSolverBase< BiCGSTAB, Matrix, Preconditioner, Traits>
{

	BOGUS_MAKE_KRYLOV_SOLVER_HEADER( BiCGSTAB )

	template < typename RhsT, typename ResT >
	Scalar solve( const RhsT &b, ResT &x ) const
	{
		typedef typename Traits::template MutableClone< RhsT >::Type Vector ;
		Vector r ;

		Scalar res = init( m_A, b, x, r ) * m_scale ;
		if( res < m_tol ) return res ;

		const Vector r0h = r ;

		Scalar rho0 = 1, rho1 ;
		Scalar alpha = 1, w = 1 ;

		Vector nu = Vector::Zero( r.rows() );
		Vector p = nu ;
		Vector s, t ( m_A.rows() ) ;
		Vector y( r.rows() ), z( t.rows() ) ;

		for( unsigned k = 0 ; k < m_maxIters ; ++k )
		{
			rho1 = r0h.dot( r ) ;

			const Scalar beta = ( rho1 / rho0 ) * ( alpha / w ) ;
			p = r + beta * ( p - w * nu ) ;
			m_P.template apply< false >( p, y ) ;
			mv_set< false >( m_A, y, nu ) ;

			alpha = rho1 / r0h.dot( nu ) ;
			s = r - alpha * nu ;
			m_P.template apply< false >( s, z ) ;
			mv_set< false >( m_A, z, t ) ;

			const Scalar nt2 = t.squaredNorm() ;
			if ( nt2 < bogus::NumTraits< Scalar >::epsilon( ) )
			{
				w = 1 ;
			} else {
				w = t.dot( s ) / nt2 ;
			}

			x += alpha*y + w*z ;
			r = s - w*t ;

			res = r.squaredNorm() * m_scale;
			this->m_callback.trigger( k, res ) ;
			if( res < m_tol ) break ;

		}

		return res ;
	}
} ;

// CGS
//! Solves ( m_A * \p x = \p b ) using the Conjugate Gradient Squared algorithm
/*! Works for non-symmetric linear systems. Convergence not guaranteed.
		Supposedly faster convergence than BiCG \a when \a converging.

		<b>Matrix-vector mults/iter: </b> 2
		<b>Preconditionner calls/iter: </b> 2
		<b>Storage requirements: </b> 7n
	*/
template < typename Matrix, typename Preconditioner, typename Traits>
struct CGS : public KrylovSolverBase< CGS, Matrix, Preconditioner, Traits>
{

	BOGUS_MAKE_KRYLOV_SOLVER_HEADER( CGS )

	template < typename RhsT, typename ResT >
	Scalar solve( const RhsT &b, ResT &x ) const
	{
		typedef typename Traits::template MutableClone< RhsT >::Type Vector ;
		Vector r ;

		Scalar res = init( m_A, b, x, r ) * m_scale ;
		if( res < m_tol ) return res ;

		const Vector r0h = r ;

		Vector u = r, p = r, q ;
		Scalar rho1, rho0, alpha, beta ;

		Vector y( m_A.cols() ), nu ( m_A.rows() ) ;

		for( unsigned k = 0 ; k < m_maxIters ; ++k )
		{
			rho1 = r0h.dot( r ) ;

			if( bogus::NumTraits< Scalar >::isZero( rho1 ) )
				break ;

			if( k > 0 )
			{
				beta = rho1 / rho0 ;
				u = r + beta *  q ;
				p = u + beta * (q + beta*p) ;
			}


			m_P.template apply< false >( p, y ) ;
			mv_set< false >( m_A, y, nu ) ;
			alpha = rho1 / r0h.dot( nu ) ;
			q = u - alpha*nu ;

			m_P.template apply< false >( u+q, y ) ;
			mv_set< false >( m_A, y, nu ) ;

			x += alpha * y ;
			r -= alpha * nu ;

			res = r.squaredNorm() * m_scale;
			this->m_callback.trigger( k, res ) ;
			if( res < m_tol ) break ;

			rho0 = rho1 ;
		}

		return res ;
	}
} ;

// GMRES
//! Solves ( m_A * \p x = \p b ) using the (restarted) Generalized Minimum Residual
/*!
		\param restart If non-zero, use the GMRES(m) restarted algorithm. Lower the storage cost,
		but slows-down or even forbid the convergence of the algorithm.

		Works for non-symmetric linear systems.
		Probably the more robust method for non symmetric systems, but with the highest storage
		cost.

		<b>Matrix-vector mults/iter: </b> 1
		<b>Preconditionner calls/iter: </b> 1
		<b>Other ops/iter: </b> 1 k*k triangular solve, 2 k*n m-v mult, 1 k*k m-v mult
		<b>Storage requirements: </b> 2*m*( n + m )

	*/
template < typename Matrix, typename Preconditioner, typename Traits>
struct GMRES : public KrylovSolverBase< GMRES, Matrix, Preconditioner, Traits>
{
protected:
	BOGUS_MAKE_KRYLOV_SOLVER_TYPEDEFS( GMRES )

	unsigned m_restart ;
public:
	GMRES( const Matrix &A,
				const Preconditioner &P,
				const Signal< unsigned, Scalar > &callback,
				Scalar tol, unsigned maxIters,
				unsigned restart = 0 )
		: Base( A, P, callback, tol, maxIters ), m_restart( restart )
	{}

	GMRES &setRestart( unsigned restart )
	{
		m_restart = restart ;
		return *this ;
	}

	template < typename RhsT, typename ResT >
	Scalar solve( const RhsT &b, ResT &x ) const
	{
		typedef typename Traits::template MutableClone< RhsT >::Type Vector ;
		typedef typename Traits::DynMatrix WorkMatrix ;
		typedef typename Traits::DynVector WorkVector ;
		typedef typename LocalProblemTraits< 2, Scalar >::Matrix Matrix22 ;

		Vector r ;

		Scalar res = init( m_A, b, x, r ) * m_scale ;
		if( res < m_tol ) return res ;

		const unsigned n = b.rows() ;

		const unsigned restart = ( m_restart == 0 ) ? n : m_restart ;
		const unsigned m = std::min( restart, m_maxIters ) ;

		// Allocate working memory

		WorkMatrix H ( std::max( n, m+1 ), m + 1 ) ;		// Hessenberg
		WorkMatrix V ( n, m + 1 ) ; // Krylov subspace basis

		WorkMatrix U ( m+1, m ) ;   // Upper triangular matrix
		WorkMatrix O ( m+1, m+1 ) ; // Orthogonal matrix such that U = O*H
		Matrix22 G ;			// Givens rotation

		WorkVector g(m+1),     // O * beta * e1
				y(m+1) ;    // Res of least-square

		// Restart loop
		unsigned globalIter = 0 ;
		do
		{
			typename WorkMatrix::ColXpr v0 ( V.col(0) ) ;
			m_P.template apply< false >( r, v0 ) ;
			Scalar beta = v0.norm() ;
			v0 /= beta ; // ~ r.normalized()

			O(0,0) = 1 ; // Initialize O to identity
			g(0,0) = beta ;

			unsigned k ;
			for( k = 0 ; k < m && res >= m_tol ; ++k )
			{

				// 1 - Arnoldi iteration
				typename WorkMatrix::ColXpr v ( V.col(k+1) ) ;
				mv_set< false >( m_A, V.col(k), r ) ;
				m_P.template apply< false >( r, v ) ;

				H.col(k  ).head( k+1 ) = V.leftCols(k+1).transpose() * v ;
				H.row(k+1).head( k   ).setZero() ;

				v -= V.leftCols( k+1 ) * H.col( k ).head( k+1 );

				const Scalar vhn = v.norm() ;
				H(k+1, k) = vhn ;

				if( vhn > bogus::NumTraits< Scalar >::epsilon( ) ) //If vhn is zero, the algorithm shall stop at this step
					V.col(k+1) /= vhn ;

				// 2 - Least squares
				// a. Grow orthogonal matrix O and vector g = O * res0 * (1,0,0, ... )'

				O.row( k+1 ).head( k+1 ).setZero() ;     // Set last row to zero
				O.col( k+1 ).head( k+1 ).setZero() ;     // Set last col to zero
				O ( k+1, k+1 ) = 1 ;
				g ( k+1 )      = 0 ;

				// a' Store temporary, before-rotation, not yet upper-triangular new U
				U.col(k).head( k+1 ) = O.topLeftCorner( k+1, k+1 ) * H.col(k).head( k+1 ) ;
				U( k+1, k ) = vhn ;

				// b. Apply givens rotation
				G.row(0) = U.col(k).template segment< 2 >( k ).transpose() ;

				const Scalar l = G.row(0).norm() ;
				G.row(0) /= l ;

				G(1, 0) = -G( 0, 1 ) ;
				G(1, 1) =  G( 0, 0 ) ;

				O.block( k, 0, 2, k+2 ).applyOnTheLeft( G );
				g.template segment< 2 >( k ).applyOnTheLeft( G ) ;

				// c. Update upper-triagular matrix U = O * H
				U.row( k+1 ).head( k+1 ).setZero() ; // Last U line to zero
				U(k,k) = l ;

				// d. Solve triangular system
				y = U.topLeftCorner( k+1, k+1 ).template triangularView< Eigen::Upper >().solve( g.head( k+1 ) ) ;

				// 3 - Update residual
				res = g( k+1 ) * g( k+1 ) * m_scale ;

				//			std::cout << " ==== Iteration " << globalIter << " + " << k <<std::endl
				//			          << "H" << std::endl
				//			          << H.topLeftCorner( k+2, k+1 )<< std::endl
				//			          << "O" << std::endl
				//			          << O.topLeftCorner( k+2, k+2 )<< std::endl
				//			          << "U" << std::endl
				//			          << U.topLeftCorner( k+2, k+1 )<< std::endl
				//			          << "V" << std::endl
				//			          << V.leftCols( k+2 ) << std::endl
				//			          << "Orthogonality" << std::endl
				//			          << ( O.topLeftCorner(k+2,k+2) * O.topLeftCorner(k+2,k+2).transpose() - Matrix::Identity( k+2, k+2 ) ).squaredNorm() << std::endl
				//			          << "Equality" << std::endl
				//			          << ( O.topLeftCorner(k+2,k+2) * H.topLeftCorner(k+2,k+1) - U.topLeftCorner(k+2,k+1) ).squaredNorm()<< std::endl
				//			          << "Solve" << std::endl
				//			          << ( U.topLeftCorner( k+1, k+1 )*y - g.segment(0,k+1) ).transpose()<< std::endl
				//			          << "res" << std::endl
				//			          << g(k+1) << std::endl
				//			          << ( m_A*x - b ).norm()
				//			          << std::endl ;


				this->m_callback.trigger( k + globalIter, res ) ;
			}

			x += V.leftCols( k ) * y.head( k ) ;
			globalIter += m ;

			if( res < m_tol || globalIter >= m_maxIters  )
				break ;

			// Restart

			r = b - m_A*x ;
			res = r.squaredNorm() * m_scale ;

		} while( res >= m_tol ) ;

		return res ;
	}
} ;

//! Solves ( m_A * \p x = \p b ) using the transpose-free Quasi Minimal Reisual method
/*! Works for non-symmetric linear systems. Convergence not guaranteed.
		Supposedly less erratic convergence than BiCG method, but faster than BiCGSTAB.

		<b>Matrix-vector mults/iter: </b> 2
		<b>Preconditionner calls/iter: </b> 2
		<b>Storage requirements: </b> 7n

		\warning This function returns an approximation of the residual instead of the real one
	*/
template < typename Matrix, typename Preconditioner, typename Traits>
struct TFQMR : public KrylovSolverBase< TFQMR, Matrix, Preconditioner, Traits>
{

	BOGUS_MAKE_KRYLOV_SOLVER_HEADER( TFQMR )

	template < typename RhsT, typename ResT >
	Scalar solve( const RhsT &b, ResT &x ) const
	{
		typedef typename Traits::template MutableClone< RhsT >::Type Vector ;
		Vector r ;

		Scalar res = init( m_A, b, x, r ) * m_scale ;
		if( res < m_tol ) return res ;

		Vector u = r, w = r ;
		const Vector& r0h = r ;

		Vector d ( m_A.rows() ), v ( m_A.rows() ),
				nu ( m_A.rows() ), y ( m_A.rows() ) ;

		d.setZero() ;

		m_P.template apply< false >( u, y ) ;
		mv_set< false >( m_A, y, v ) ;

		Scalar tau2 = res/m_scale ;
		Scalar rho = tau2, rho1 ; // rho = r.dot( r0h )
		Scalar theta2 = 0, eta = 0, alpha = 0,  beta, c2 ;

		for( unsigned k = 0 ; k < m_maxIters ; ++k )
		{


			const bool odd = k&1u;
			if( !odd )
			{
				alpha = rho / r0h.dot( v ) ;
			}

			m_P.template apply< false >( u, y ) ;
                        mv_set< false >( m_A, y, nu ) ;
			w -= alpha * nu ;
			d = y + ( theta2 / alpha ) * eta * d ;

			theta2 = w.squaredNorm() ;
			theta2 /= tau2 ;

			c2 = 1. / ( 1 + theta2 ) ;
			tau2 *= theta2 * c2 ;
			eta = c2 * alpha ;

			x += eta * d ;
			// m_A.template multiply< false >( d, r, -eta, 1 ) ;
			res = tau2 * m_scale ; // Approx

			this->m_callback.trigger( k, res ) ;
			if( res < m_tol ) break ;

			if( odd )
			{
				rho1 = w.dot( r0h ) ;
				beta = rho1 / rho ;

				u = w + beta * u ;

				v = nu + beta * v ;
				m_P.template apply< false >( u, y ) ;
                                v *= beta ;
                                mv_add< false >( m_A, y, v, 1 ) ;

				rho = rho1 ;
			} else {
				u -= alpha * v ;
			}

			//		rho0 = rho1 ;
		}

		return res ;
	}
};


} //namespace krylov

} //namespace bogus

#endif
