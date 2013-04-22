/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver          
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>                       
 *
 * This Source Code Form is subject to the terms of the Mozilla Public 
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef BOGUS_BLOCK_SERIALIZATION_HPP
#define BOGUS_BLOCK_SERIALIZATION_HPP

#include "SparseBlockMatrix.hpp"

namespace boost {
namespace serialization {

template<typename Archive, typename Index, typename BlockPtr  >
inline void serialize(
	   Archive & ar,
	   bogus::SparseBlockIndex<false, Index, BlockPtr > & index,
	   const unsigned int file_version
   )
{
	(void) file_version ;
	ar & index.valid ;
	ar & index.innerOffsets ;
	ar & index.outer ;
}

template<typename Archive, typename Index, typename BlockPtr  >
inline void serialize(
	   Archive & ar,
	   bogus::SparseBlockIndex<true, Index, BlockPtr > & index,
	   const unsigned int file_version
   )
{
	(void) file_version ;
	ar & index.valid ;
	ar & index.innerOffsets ;
	ar & index.base ;
	ar & index.inner ;
	ar & index.outer ;
}

} //serialization
} //boost

namespace bogus {

template < typename Derived >
template < typename Archive >
void SparseBlockMatrixBase< Derived >::serialize(
	   Archive & ar,
	   const unsigned int file_version
   )
{
	(void) file_version ;

	ar & m_rows ;
	ar & m_cols ;
	ar & m_blocks ;
	ar & m_nBlocks ;
	ar & m_majorIndex ;
	ar & m_minorIndex ;
	ar & m_transposeIndex ;

}

} //bogus

#endif
