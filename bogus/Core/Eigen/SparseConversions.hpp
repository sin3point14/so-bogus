#ifndef BOGUS_EIGEN_SPARSE_CONVERSIONS_HPP
#define BOGUS_EIGEN_SPARSE_CONVERSIONS_HPP

#include "../Block/SparseBlockMatrix.hpp"

#include <map>

// Conversions to Sparse Matrix
namespace bogus
{

template< typename EigenDerived, typename BogusDerived >
void convert( const Eigen::SparseMatrixBase< EigenDerived >& source,
			  SparseBlockMatrixBase< BogusDerived >& dest )
{
	typedef BlockMatrixTraits< BogusDerived > Traits ;
    typedef typename Traits::Index Index ;
    typedef typename Traits::BlockPtr BlockPtr ;
    const Index RowsPerBlock = Traits::BlockType::RowsAtCompileTime ;
	const Index ColsPerBlock = Traits::BlockType::ColsAtCompileTime ;

	assert( RowsPerBlock != (Index) -1 ) ;
	assert( ColsPerBlock != (Index) -1 ) ;

	assert( ( (bool) Eigen::SparseMatrixBase< EigenDerived >::IsRowMajor ) ^
			( (bool) Traits::is_col_major ) ) ;

	assert( 0 == ( source.rows() % RowsPerBlock ) ) ;
	assert( 0 == ( source.cols() % ColsPerBlock ) ) ;

	dest.clear() ;
	dest.setRows( source.rows() / RowsPerBlock, RowsPerBlock ) ;
	dest.setCols( source.cols() / ColsPerBlock, ColsPerBlock ) ;

	const Index blockSize = Traits::is_col_major ? ColsPerBlock : RowsPerBlock ;
	for( Index outer = 0 ; outer < dest.majorIndex().outerSize() ; ++outer )
	{
		// I - compute non-zero blocks
        std::map < Index, BlockPtr > nzBlocks ;

		for( unsigned i = 0 ; i < blockSize ; ++i )
		{
            for( typename EigenDerived::InnerIterator innerIt( source.derived(), outer*blockSize + i ) ;
				 innerIt ; ++innerIt )
			{
				const Index blockId = (Index) ( innerIt.index() ) / blockSize  ;
                if( Traits::is_symmetric && blockId > outer ) break ;
				nzBlocks[ blockId ] = 0 ;
			}
		}

        // II - Insert them in block mat
        for( typename std::map< Index, BlockPtr >::iterator bIt = nzBlocks.begin() ; bIt != nzBlocks.end() ; ++bIt )
		{
            bIt->second = (BlockPtr) dest.nBlocks() ;
            dest.insertBackOuterInner( outer, bIt->first ).setZero() ;
        }

        // III - copy values
		for( unsigned i = 0 ; i < blockSize ; ++i )
		{
            for( typename EigenDerived::InnerIterator innerIt( source.derived(), outer*blockSize + i ) ;
				 innerIt ; ++innerIt )
			{
				const Index blockId = (Index) ( innerIt.index() ) / blockSize  ;
                if( Traits::is_symmetric && blockId > outer ) break ;
				const Index bin = innerIt.index() - blockId * blockSize  ;
				const Index row = Traits::is_col_major ? bin : i ;
				const Index col = Traits::is_col_major ? i : bin ;

                dest.block( nzBlocks[ blockId ] ) ( row, col ) = innerIt.value() ;
			}
		}

	}

	dest.finalize() ;

}

} // naemspace bogus


#endif
