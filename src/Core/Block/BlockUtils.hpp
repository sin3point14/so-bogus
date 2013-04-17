#ifndef BOGUS_BLOCK_UTILS
#define BOGUS_BLOCK_UTILS

namespace bogus {

template < bool Symmetric, bool DoTranspose >
struct BlockTranspose{
	template < typename BlockT >
	static const BlockT& get( const BlockT& src, bool )
	{ return src ; }
} ;
template <  >
struct BlockTranspose< false, true > {
	//SFINAE

	template < typename BlockT >
	static typename BlockT::ConstTransposeReturnType get( const BlockT& src, bool )
	{ return src.transpose() ; }

	template < typename BlockT >
	static typename BlockTransposeTraits< typename BlockT::Base >::ReturnType get( const BlockT& src, bool )
	{ return transpose_block( src ) ; }

	template < typename BlockT >
	static typename BlockTransposeTraits< BlockT >::ReturnType get( const BlockT& src, bool )
	{ return transpose_block( src ) ; }

	template < typename BlockT >
	static typename SelfTransposeTraits< BlockT >::ReturnType get( const BlockT& src, bool )
	{ return src ; }
} ;

template < bool DoTranspose >
struct BlockTranspose< true, DoTranspose > {
	template < typename BlockT >
	static BlockT get( const BlockT& src, bool afterDiag )
	{ return afterDiag ? BlockT( transpose_block( src ) ) : src ; }
} ;


}

#endif
