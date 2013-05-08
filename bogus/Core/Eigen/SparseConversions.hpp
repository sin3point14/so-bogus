#ifndef BOGUS_EIGEN_SPARSE_CONVERSIONS_HPP
#define BOGUS_EIGEN_SPARSE_CONVERSIONS_HPP

#include "../Block/SparseBlockMatrix.hpp"

#include <map>
#include <iostream>

// Conversions to Sparse Matrix
namespace bogus
{

template< typename EigenDerived, typename BogusDerived >
void convert( const Eigen::SparseMatrixBase< EigenDerived >& source,
			  SparseBlockMatrixBase< BogusDerived >& dest )
{
	typedef BlockMatrixTraits< BogusDerived > Traits ;
	typedef typename Traits::Index Index ;
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
		std::map < Index, Index > nzBlocks ;

		for( unsigned i = 0 ; i < blockSize ; ++i )
		{
			for( typename EigenDerived::InnerIterator innerIt( source, outer*blockSize + i ) ;
				 innerIt ; ++innerIt )
			{
				const Index blockId = (Index) ( innerIt.index() ) / blockSize  ;
				//if( Traits::is_symmetric && blockId > outer ) break ;
				nzBlocks[ blockId ] = 0 ;
			}
		}

		typename BlockContainerTraits< typename Traits::BlockType >::Type values ;
		values.resize( nzBlocks.size() ) ;

		Index blockId = 0 ;
		for( typename std::map< Index, Index >::iterator bIt = nzBlocks.begin() ; bIt != nzBlocks.end() ; ++bIt )
		{
			values[ blockId ].setZero() ;
			bIt->second = blockId ++ ;
		}

		// II - copy values
		for( unsigned i = 0 ; i < blockSize ; ++i )
		{
			for( typename EigenDerived::InnerIterator innerIt( source, outer*blockSize + i ) ;
				 innerIt ; ++innerIt )
			{
				const Index blockId = (Index) ( innerIt.index() ) / blockSize  ;
				//if( Traits::is_symmetric && blockId > outer ) break ;
				const Index bin = innerIt.index() - blockId * blockSize  ;
				const Index row = Traits::is_col_major ? bin : i ;
				const Index col = Traits::is_col_major ? i : bin ;

				std::cout << blockId << " " << outer << " " << i << " " << innerIt.value() << std::endl ;
				values[ nzBlocks[ blockId ] ] ( row, col ) = innerIt.value() ;
			}
		}


		// III - Insert them in block mat
		for( typename std::map< Index, Index >::const_iterator bIt = nzBlocks.begin() ; bIt != nzBlocks.end() ; ++bIt )
		{
			dest.insertBackOuterInner( outer, bIt->first ) = values[ bIt->second ] ;
		}

	}

	dest.finalize() ;

}

} // naemspace bogus


#endif
