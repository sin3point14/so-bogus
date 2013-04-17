#ifndef BOGUS_BLOCK_UTILS
#define BOGUS_BLOCK_UTILS

namespace bogus {

template < bool DoTranspose >
struct BlockGetter {
	template < typename BlockT >
	static const BlockT& get( const BlockT& src, bool = false )
	{ return src ; }
} ;
template < >
struct BlockGetter< true > {
	//SFINAE

	template < typename BlockT >
	static typename BlockT::ConstTransposeReturnType get( const BlockT& src, bool = false )
	{ return src.transpose() ; }

	template < typename BlockT >
	static typename BlockTransposeTraits< typename BlockT::Base >::ReturnType get( const BlockT& src, bool = false )
	{ return transpose_block( src ) ; }

	template < typename BlockT >
	static typename BlockTransposeTraits< BlockT >::ReturnType get( const BlockT& src, bool = false )
	{ return transpose_block( src ) ; }

	template < typename BlockT >
	static typename SelfTransposeTraits< BlockT >::ReturnType get( const BlockT& src, bool = false )
	{ return src ; }
} ;

template < bool RuntimeCheck, bool DoTranspose >
struct BlockTransposeOption : public BlockGetter< DoTranspose >{
} ;

template < bool IgnoredDoTranspose >
struct BlockTransposeOption< true, IgnoredDoTranspose > {
	template < typename BlockT >
	static BlockT get( const BlockT& src, bool doTranspose = false )
	{ return doTranspose ? BlockT( transpose_block( src ) ) : src ; }
} ;


}

#endif
