#ifndef BOGUS_KRYLOV_METHODS_HPP
#define BOGUS_KRYLOV_METHODS_HPP

#include "../Utils/LinearSolverBase.hpp"
#include "../Utils/Signal.hpp"
#include "../BlockSolvers.fwd.hpp"

#define BOGUS_KRYLOV_METHODS \
    BOGUS_PROCESS_KRYLOV_METHOD(CG      )\
    BOGUS_PROCESS_KRYLOV_METHOD(BiCG    )\
    BOGUS_PROCESS_KRYLOV_METHOD(BiCGSTAB)\
    BOGUS_PROCESS_KRYLOV_METHOD(CGS     )\
    BOGUS_PROCESS_KRYLOV_METHOD(GMRES   )\
    BOGUS_PROCESS_KRYLOV_METHOD(TFQMR   )


namespace bogus {

template< template < typename, typename, typename, typename > class Method,
          typename Matrix, typename Preconditioner, typename Vector, typename Scalar >
struct KrilovSolverBase
        : public LinearSolverBase< Method< Matrix, Preconditioner, Vector, Scalar > >
{
    typedef Method< Matrix, Preconditioner, Vector, Scalar > Derived ;
    typedef LinearSolverBase< Derived > Base ;

    const Matrix &m_A ;
    const Preconditioner &m_P ;
    const Signal< unsigned, Scalar >  &m_callback ;

    unsigned m_maxIters;
    Scalar m_tol ;

    KrilovSolverBase( const Matrix &A,
                      const Preconditioner &P,
                      const Signal< unsigned, Scalar > callback,
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

private:
    Scalar m_scale ;
} ;

// Forward decl + create traits

#define BOGUS_PROCESS_KRYLOV_METHOD( MethodName ) \
    namespace krylov {\
        template< typename Matrix, typename Preconditioner, typename Vector, typename Scalar > \
        struct MethodName ;\
    } \
    template< typename Matrix, typename Preconditioner, typename Vector, typename Scalar > \
    struct LinearSolverTraits< krylov::MethodName< Matrix, Preconditioner, Vector, Scalar > > \
    { \
      typedef Matrix MatrixType ; \
      template < typename RhsT > struct Result { \
          typedef Vector Type ; \
      } ; \
    } ;

BOGUS_KRYLOV_METHODS
#undef BOGUS_PROCESS_KRYLOV_METHOD


#define BOGUS_MAKE_KRYLOV_SOLVER_HEADER( MethodName ) \
    typedef KrilovSolverBase< MethodName, Matrix, Preconditioner, Vector, Scalar > Base ; \
    \
    MethodName( const Matrix &A, \
    const Preconditioner &P, \
    const Signal< unsigned, Scalar > callback, \
    Scalar tol, unsigned maxIters ) \
    : Base( A, P, callback, tol, maxIters ) \
{}  \



namespace krylov {


// Conjugate Gradient

//! Solves ( m_matrix * \p x = \p b ) using the Conjugate Gradient algorithm
/*! For symmetric matrices only. Converges for positive definite linear systems.

        <b>Matrix-vector mults/iter: </b> 1
        <b>Preconditionner calls/iter: </b> 1
        <b>Storage requirements: </b> 4n
    */
template < typename Matrix, typename Preconditioner, typename Vector, typename Scalar >
struct ConjugateGradient : public KrilovSolverBase< ConjugateGradient, Matrix, Preconditioner, Vector, Scalar >
{

    BOGUS_MAKE_KRYLOV_SOLVER_HEADER( ConjugateGradient )

    template < typename RhsT, typename ResT >
    Scalar solve( const RhsT &b, ResT &x ) const
    {
        Vector r ;

        Scalar res = init( b, x, r ) ;
        if( res < m_tol ) return res ;

        Vector z( r.rows() ) ;
        m_preconditioner.template apply< false >( r, z ) ;
        Vector p = z;

        Scalar zr0 = r.dot( z ) ;
        Scalar zr1 ;

        Vector Mp( m_matrix->rows() ) ;

        for( unsigned k = 0 ; k < m_maxIters ; ++k )
        {
            m_matrix->template multiply< false >( p, Mp ) ;
            const Scalar alpha = zr0 / ( p.dot( Mp ) ) ;
            x += alpha * p ;
            r -= alpha * Mp ;

            res = r.squaredNorm() * m_scale ;
            this->m_callback.trigger( k, res ) ;
            if( res < m_tol ) break ;

            m_preconditioner.template apply< false >( r, z ) ;
            zr1 = z.dot( r ) ;

            p = z + ( zr1 / zr0 ) * p ;

            zr0 = zr1 ;
        }

        return res ;
    }
} ;

// BiConjugate Gradient
//! Solves ( m_matrix * \p x = \p b ) using the BiConjugate Gradient algorithm
/*! Works for non-symmetric linear systems. Convergence not guaranteed.
        Requires ability to perform transpose multiplication and preconditioning

        <b>Matrix-vector mults/iter: </b> 2 ( inc. 1 transpose )
        <b>Preconditionner calls/iter: </b> 2 ( inc. 1 transpose )
        <b>Storage requirements: </b> 8n
    */
template < typename Matrix, typename Preconditioner, typename Vector, typename Scalar >
struct BiCG : public KrilovSolverBase< BiCG, Matrix, Preconditioner, Vector, Scalar >
{

    BOGUS_MAKE_KRYLOV_SOLVER_HEADER( BiCG )

    template < typename BlockMatrixType, template< typename BlockMatrixT > class PreconditionerType >
    template < typename RhsT, typename ResT >
    typename Krylov< BlockMatrixType, PreconditionerType >::Scalar
    Krylov< BlockMatrixType, PreconditionerType >::solve_BiCG( const RhsT &b, ResT &x ) const
    {
        typedef typename GlobalProblemTraits::DynVector Vector ;
        Vector r ;

        Scalar res = init( b, x, r ) ;
        if( res < m_tol ) return res ;

        Vector x_ = x ;

        Vector r_ = b ;
        m_matrix->template multiply< true >( x_, r_, 1, 1 );

        Vector p( r.rows() ), p_ ( r_.rows() ) ;
        m_preconditioner.template apply< false >( r,  p ) ;
        m_preconditioner.template apply< true >( r_, p_ ) ;

        Vector z = p, z_ = p_ ;
        Scalar zr0 = r_.dot( z ) ;
        Scalar zr1 ;

        Vector Mp( m_matrix->rows() ) ;

        for( unsigned k = 0 ; k < m_maxIters ; ++k )
        {
            m_matrix->template multiply< false >( p, Mp, 1, 0 ) ;
            const Scalar alpha = zr0 / ( p_.dot( Mp ) ) ;
            x  += alpha * p  ;
            x_ += alpha * p_ ;

            r  -= alpha * Mp ;
            m_matrix->template multiply< true >( p_, r_, -alpha, 1 ) ;

            res = r.squaredNorm() * m_scale ;
            this->m_callback.trigger( k, res ) ;
            if( res < m_tol ) break ;

            m_preconditioner.template apply< false >( r, z ) ;
            zr1 = z.dot( r_ ) ;

            const Scalar beta = ( zr1 / zr0 ) ;

            p = z + beta * p ;
            m_preconditioner.template apply< true >( r_, z_ ) ;
            p_ = z_ + beta * p_ ;

            zr0 = zr1 ;
        }


        return res ;
    }
} ;

// BiCG STAB
//! Solves ( m_matrix * \p x = \p b ) using the BiConjugate Gradient stabilized algorithm
/*! Works for non-symmetric linear systems. Convergence not guaranteed.
        Supposedly less erratic convergence than BiCG method.

        <b>Matrix-vector mults/iter: </b> 2
        <b>Preconditionner calls/iter: </b> 2
        <b>Storage requirements: </b> 8n
    */
template < typename Matrix, typename Preconditioner, typename Vector, typename Scalar >
struct BiCGSTAB : public KrilovSolverBase< BiCGSTAB, Matrix, Preconditioner, Vector, Scalar >
{

    BOGUS_MAKE_KRYLOV_SOLVER_HEADER( BiCGSTAB )

    template < typename BlockMatrixType, template< typename BlockMatrixT > class PreconditionerType >
    template < typename RhsT, typename ResT >
    typename Krylov< BlockMatrixType, PreconditionerType >::Scalar
    Krylov< BlockMatrixType, PreconditionerType >::solve_BiCGSTAB( const RhsT &b, ResT &x ) const
    {
        typedef typename GlobalProblemTraits::DynVector Vector ;
        Vector r ;

        Scalar res = init( b, x, r ) ;
        if( res < m_tol ) return res ;

        const Vector r0h = r ;

        Scalar rho0 = 1, rho1 ;
        Scalar alpha = 1, w = 1 ;

        Vector nu = Vector::Zero( r.rows() );
        Vector p = nu ;
        Vector s, t ( m_matrix->rows() ) ;
        Vector y( r.rows() ), z( t.rows() ) ;

        for( unsigned k = 0 ; k < m_maxIters ; ++k )
        {
            rho1 = r0h.dot( r ) ;

            const Scalar beta = ( rho1 / rho0 ) * ( alpha / w ) ;
            p = r + beta * ( p - w * nu ) ;
            m_preconditioner.template apply< false >( p, y ) ;
            m_matrix->template multiply< false >( y, nu ) ;

            alpha = rho1 / r0h.dot( nu ) ;
            s = r - alpha * nu ;
            m_preconditioner.template apply< false >( s, z ) ;
            m_matrix->template multiply< false >( z, t ) ;

            const Scalar nt2 = t.squaredNorm() ;
            if ( nt2 < NumTraits< Scalar >::epsilon( ) )
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
//! Solves ( m_matrix * \p x = \p b ) using the Conjugate Gradient Squared algorithm
/*! Works for non-symmetric linear systems. Convergence not guaranteed.
        Supposedly faster convergence than BiCG \a when \a converging.

        <b>Matrix-vector mults/iter: </b> 2
        <b>Preconditionner calls/iter: </b> 2
        <b>Storage requirements: </b> 7n
    */
template < typename Matrix, typename Preconditioner, typename Vector, typename Scalar >
struct CGS : public KrilovSolverBase< CGS, Matrix, Preconditioner, Vector, Scalar >
{

    BOGUS_MAKE_KRYLOV_SOLVER_HEADER( CGS )


    template < typename BlockMatrixType, template< typename BlockMatrixT > class PreconditionerType >
    template < typename RhsT, typename ResT >
    typename Krylov< BlockMatrixType, PreconditionerType >::Scalar
    Krylov< BlockMatrixType, PreconditionerType >::solve_CGS( const RhsT &b, ResT &x ) const
    {
        typedef typename GlobalProblemTraits::DynVector Vector ;
        Vector r ;

        Scalar res = init( b, x, r ) ;
        if( res < m_tol ) return res ;

        const Vector r0h = r ;

        Vector u = r, p = r, q ;
        Scalar rho1, rho0, alpha, beta ;

        Vector y( m_matrix->cols() ), nu ( m_matrix->rows() ) ;

        for( unsigned k = 0 ; k < m_maxIters ; ++k )
        {
            rho1 = r0h.dot( r ) ;

            if( NumTraits< Scalar >::isZero( rho1 ) )
                break ;

            if( k > 0 )
            {
                beta = rho1 / rho0 ;
                u = r + beta *  q ;
                p = u + beta * (q + beta*p) ;
            }


            m_preconditioner.template apply< false >( p, y ) ;
            m_matrix->template multiply< false >( y, nu ) ;
            alpha = rho1 / r0h.dot( nu ) ;
            q = u - alpha*nu ;

            m_preconditioner.template apply< false >( u+q, y ) ;
            m_matrix->template multiply< false >( y, nu ) ;

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
//! Solves ( m_matrix * \p x = \p b ) using the (restarted) Generalized Minimum Residual
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
template < typename Matrix, typename Preconditioner, typename Vector, typename Scalar >
struct CMRES : public KrilovSolverBase< GMRES, Matrix, Preconditioner, Vector, Scalar >
{

    BOGUS_MAKE_KRYLOV_SOLVER_HEADER( GMRES )


    template < typename BlockMatrixType, template< typename BlockMatrixT > class PreconditionerType >
    template < typename RhsT, typename ResT >
    typename Krylov< BlockMatrixType, PreconditionerType >::Scalar
    Krylov< BlockMatrixType, PreconditionerType >::solve_GMRES( const RhsT &b, ResT &x, unsigned restart ) const
    {
        typedef typename GlobalProblemTraits::DynVector Vector ;
        typedef typename GlobalProblemTraits::DynMatrix Matrix ;
        typedef typename LocalProblemTraits< 2, Scalar >::Matrix Matrix22 ;


        Vector r ;

        Scalar res = init( b, x, r ) ;
        if( res < m_tol ) return res ;

        const unsigned n = b.rows() ;

        if( restart == 0 ) restart = n ;
        const unsigned m = std::min( restart, m_maxIters ) ;

        // Allocate working memory

        Matrix H ( std::max( n, m+1 ), m + 1 ) ;		// Hessenberg
        Matrix V ( n, m + 1 ) ; // Krylov subspace basis

        Matrix U ( m+1, m ) ;   // Upper triangular matrix
        Matrix O ( m+1, m+1 ) ; // Orthogonal matrix such that U = O*H
        Matrix22 G ;			// Givens rotation

        Vector g(m+1),     // O * beta * e1
                y(m+1) ;    // Res of least-square

        // Restart loop
        unsigned globalIter = 0 ;
        do
        {
            typename Matrix::ColXpr v0 ( V.col(0) ) ;
            m_preconditioner.template apply< false >( r, v0 ) ;
            Scalar beta = v0.norm() ;
            v0 /= beta ; // ~ r.normalized()

            O(0,0) = 1 ; // Initialize O to identity
            g(0,0) = beta ;

            unsigned k ;
            for( k = 0 ; k < m && res >= m_tol ; ++k )
            {

                // 1 - Arnoldi iteration
                typename Matrix::ColXpr v ( V.col(k+1) ) ;
                m_matrix->template multiply< false >( V.col(k), r ) ;
                m_preconditioner.template apply< false >( r, v ) ;

                H.col(k  ).head( k+1 ) = V.leftCols(k+1).transpose() * v ;
                H.row(k+1).head( k   ).setZero() ;

                v -= V.leftCols( k+1 ) * H.col( k ).head( k+1 );

                const Scalar vhn = v.norm() ;
                H(k+1, k) = vhn ;

                if( vhn > NumTraits< Scalar >::epsilon( ) ) //If vhn is zero, the algorithm shall stop at this step
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
                //			          << ( (*m_matrix)*x - b ).norm()
                //			          << std::endl ;


                this->m_callback.trigger( k + globalIter, res ) ;
            }

            x += V.leftCols( k ) * y.head( k ) ;
            globalIter += restart ;

            if( res < m_tol || globalIter >= m_maxIters  )
                break ;

            // Restart

            r = b - (*m_matrix)*x ;
            res = r.squaredNorm() * m_scale ;

        } while( res >= m_tol ) ;

        return res ;
    }
} ;

//! Solves ( m_matrix * \p x = \p b ) using the transpose-free Quasi Minimal Reisual method
/*! Works for non-symmetric linear systems. Convergence not guaranteed.
        Supposedly less erratic convergence than BiCG method, but faster than BiCGSTAB.

        <b>Matrix-vector mults/iter: </b> 2
        <b>Preconditionner calls/iter: </b> 2
        <b>Storage requirements: </b> 7n

        \warning This function returns an approximation of the residual instead of the real one
    */
template < typename Matrix, typename Preconditioner, typename Vector, typename Scalar >
struct TFQMR : public KrilovSolverBase< TFQMR, Matrix, Preconditioner, Vector, Scalar >
{

    BOGUS_MAKE_KRYLOV_SOLVER_HEADER( TFQMR )

    template < typename BlockMatrixType, template< typename BlockMatrixT > class PreconditionerType >
    template < typename RhsT, typename ResT >
    typename Krylov< BlockMatrixType, PreconditionerType >::Scalar
    Krylov< BlockMatrixType, PreconditionerType >::solve_TFQMR( const RhsT &b, ResT &x ) const
    {
        typedef typename GlobalProblemTraits::DynVector Vector ;
        Vector r ;

        Scalar res = init( b, x, r ) ;
        if( res < m_tol ) return res ;

        Vector u = r, w = r ;
        const Vector& r0h = r ;

        Vector d ( m_matrix->rows() ), v ( m_matrix->rows() ),
                nu ( m_matrix->rows() ), y ( m_matrix->rows() ) ;

        d.setZero() ;

        m_preconditioner.template apply< false >( u, y ) ;
        m_matrix->template multiply< false >( y, v ) ;

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

            m_preconditioner.template apply< false >( u, y ) ;
            m_matrix->template multiply< false >( y, nu ) ;
            w -= alpha * nu ;
            d = y + ( theta2 / alpha ) * eta * d ;

            theta2 = w.squaredNorm() ;
            theta2 /= tau2 ;

            c2 = 1. / ( 1 + theta2 ) ;
            tau2 *= theta2 * c2 ;
            eta = c2 * alpha ;

            x += eta * d ;
            // m_matrix->template multiply< false >( d, r, -eta, 1 ) ;
            res = tau2 * m_scale ; // Approx

            this->m_callback.trigger( k, res ) ;
            if( res < m_tol ) break ;

            if( odd )
            {
                rho1 = w.dot( r0h ) ;
                beta = rho1 / rho ;

                u = w + beta * u ;

                v = nu + beta * v ;
                m_preconditioner.template apply< false >( u, y ) ;
                m_matrix->template multiply< false >( y, v, 1, beta ) ;

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
