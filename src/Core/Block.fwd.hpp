#ifndef BOGUS_BLOCK_FWD_HPP
#define BOGUS_BLOCK_FWD_HPP

namespace bogus {

namespace flags
{
	enum {
		NONE = 0,
		COMPRESSED = 0x1,
		COL_MAJOR = 0x2,
		SYMMETRIC = 0x4
	} ;
}

template< typename Derived >
struct BlockMatrixTraits ;

template < typename Derived >
struct BlockObjectBase ;

template < typename Derived >
class BlockMatrixBase ;

template < typename Derived >
class SparseBlockMatrixBase ;

template < typename BlockT, int Flags = flags::NONE >
class SparseBlockMatrix  ;


}

#endif
