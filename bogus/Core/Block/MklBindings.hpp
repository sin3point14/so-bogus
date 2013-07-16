#ifndef BOGUS_BLOCK_MKL_BINDINGS_HPP
#define BOGUS_BLOCK_MKL_BINDINGS_HPP

#include "SparseBlockMatrix.hpp"
#include "CompressedSparseBlockIndex.hpp"

#include <mkl.h>
#include <mkl_spblas.h>

// Creates a compile error with boost
#ifdef P4
#undef P4
#endif

#include <iostream>

namespace bogus
{

namespace mkl
{

template< typename Scalar >
struct bindings {} ;

template< >
struct bindings< double >
{
    typedef double Scalar ;
    static void bsrmv (char *transa, MKL_INT *m, MKL_INT *k, MKL_INT *lb, Scalar *alpha, char *matdescra, Scalar  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, Scalar *x, Scalar *beta, Scalar *y)
    {
        mkl_dbsrmv( transa, m, k, lb, alpha,
                    matdescra, val, indx, pntrb, pntre,
                    x, beta, y ) ;
    }
} ;

/* FIXME sbsrmv no supported on older MKL versions ; do version detection
template< >
struct bindings< float >
{
    typedef float Scalar ;
    static void bsrmv (char *transa, MKL_INT *m, MKL_INT *k, MKL_INT *lb, Scalar *alpha, char *matdescra, Scalar  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, Scalar *x, Scalar *beta, Scalar *y)
    {
        mkl_sbsrmv( transa, m, k, lb, alpha,
                    matdescra, val, indx, pntrb, pntre,
                    x, beta, y ) ;
    }
} ; */


template< typename Scalar, typename BlockPtr >
void bsrmv(
        bool Symmetric, bool Transpose, MKL_INT Dimension,
        const SparseBlockIndex< true, MKL_INT, BlockPtr >& index, const Scalar *data,
        const Scalar *rhs, int rhsCols, Scalar *res, Scalar alpha, Scalar beta )
{

    char matdescra[4] = { Symmetric ? 'S' : 'G', 'L', 'N', 'C'} ;

    int m = index.outerSize() ;
    int k = index.innerSize() ;

    char transa = Transpose ? 'T' : 'N' ;
    int lb = Dimension ;
    Scalar *a = const_cast< Scalar* > ( data ) ;
    Scalar *x = const_cast< Scalar* > ( rhs ) ;
    Scalar *y = res ;

    MKL_INT* rowIndex = const_cast< MKL_INT* >( index.rowIndex() ) ;
    MKL_INT* columns  = const_cast< MKL_INT* >( index.columns () ) ;

    MKL_INT* pntrb = rowIndex ;
    MKL_INT* pntre = rowIndex+1 ;
    MKL_INT* indx = columns  ;

    for( int i = 0 ; i < rhsCols ; ++i )
    {
        bindings< Scalar >::bsrmv( &transa, &m, &k, &lb, &alpha,
                                   matdescra, a, indx, pntrb, pntre,
                                   x + i*Dimension*k, &beta, y + i*Dimension*k ) ;
    }
}

}

template <>
struct SparseBlockMatrixOpProxy< true, true, double, MKL_INT >
{
    typedef double Scalar ;

    template < bool Transpose, typename Derived, typename RhsT, typename ResT >
    static void multiply( const SparseBlockMatrixBase< Derived >& matrix,  const RhsT& rhs, ResT& res,
                          Scalar alpha, Scalar beta )
    {
        typedef BlockMatrixTraits< Derived > Traits ;

        mkl::bsrmv< Scalar >
                ( Traits::is_symmetric, Transpose, Derived::RowsPerBlock,
                  matrix.majorIndex(), data_pointer( matrix.data()[0] ),
                  rhs.data(), rhs.cols(), res.data(), alpha, beta ) ;
    }

    template < typename Derived, typename RhsT, typename ResT >
    static void splitRowMultiply( const SparseBlockMatrixBase< Derived >& matrix, typename Derived::Index row, const RhsT& rhs, ResT& res  )
    {
        typedef BlockMatrixTraits< Derived > Traits ;
        SparseBlockSplitRowMultiplier< Traits::is_symmetric, !Traits::is_col_major >
                ::splitRowMultiply( matrix, row, rhs, res )  ;
    }
} ;


} //namespace bogus

#endif
