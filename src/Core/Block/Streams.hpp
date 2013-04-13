#ifndef BOGUS_BLOCK_STREAMS_HPP
#define BOGUS_BLOCK_STREAMS_HPP

#include <iostream>
#include "SparseBlockMatrix.hpp"

template < typename Derived >
std::ostream& operator<<( std::ostream &out, const bogus::SparseBlockMatrixBase< Derived > &sbm )
{
	out << " Total rows: " << sbm.rows() << " / cols: " << sbm.cols() << std::endl ;
	for ( unsigned i = 0 ; i < sbm.majorIndex().outerSize() ; ++ i )
	{
		out << i << ": " ;
		for( typename bogus::SparseBlockMatrixBase< Derived >::SparseIndexType::InnerIterator it( sbm.majorIndex(), i ) ;
			 it ; ++ it )
		{
			out << "(" << it.inner() << ";" << it.ptr() << ")" ;
		}
		out << std::endl ;
	}
	out << " Blocks (" << sbm.nBlocks() << ")" << std::endl ;
	for ( unsigned i = 0 ; i < sbm.nBlocks() ; ++ i )
	{
		out << sbm.block(i) << std::endl ;
		out << " --- " << std::endl ;
	}
	return out ;
}

#endif
