/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver          
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>                       
 *
 * This Source Code Form is subject to the terms of the Mozilla Public 
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


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
