/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver          
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>                       
 *
 * This Source Code Form is subject to the terms of the Mozilla Public 
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Standard shared_ptr -- BOGUS_SHARED_PTR_NS should be set to std, std::tr1 or boost
#ifdef BOGUS_SHARED_PTR_NS
#ifndef BOGUS_SHARED_PTR
#include <memory>
#define BOGUS_SHARED_PTR( Type, Name ) BOGUS_SHARED_PTR_NS::shared_ptr< Type > Name
#endif
#endif

#ifndef BOGUS_NAIVE_SHARED_PTR_HPP
#define BOGUS_NAIVE_SHARED_PTR_HPP

#ifndef BOGUS_SHARED_PTR
#define BOGUS_SHARED_PTR( Type, Name ) NaiveSharedPtr< Type > Name
#endif

namespace bogus {

//! Naive reference-counting Shared Pointer
/*! There's no reason to use this implementation over any other one from boost, tr1, or c++11,
  except for c++98 compability. Fortunately it is only used for an edge case ( see EigenSparseLinearSolvers.hpp )
*/
template < typename T >
class NaiveSharedPtr
{
private:

	T* m_instance ;
	unsigned* m_refCount ;

public:

	NaiveSharedPtr( T * instance = 0 )
	{
		acquire( instance ) ;
	}

	NaiveSharedPtr( const NaiveSharedPtr< T >& rhs )
		: m_refCount( NULL )
	{
		add_ref( rhs ) ;
	}

	const NaiveSharedPtr& operator=( const NaiveSharedPtr< T >& rhs )
	{
		if( this != &rhs )
		{
			release() ;
			add_ref( rhs ) ;
		}
		return *this ;
	}

	~NaiveSharedPtr()
	{
		release() ;
	}

	void reset( T * instance = 0 )
	{
		release() ;
		acquire( instance ) ;
	}

	void release()
	{
		if( m_instance )
		{
			if( 0 == --*m_refCount )
			{
				delete m_instance ;
				delete m_refCount ;
			}

			m_instance = NULL ;
			m_refCount = NULL ;
		}
	}

	T& operator * ( ) {
		return *m_instance;
	}
	T const& operator * ( ) const {
		return *m_instance;
	}
	T* operator -> ( ) {
		return m_instance;
	}
	T const* operator -> ( ) const {
		return m_instance;
	}
	T* get() {
		return m_instance ;
	}
	const T* get() const {
		return m_instance ;
	}

	operator bool() const {
		return m_instance ;
	}

private:

	void add_ref( const NaiveSharedPtr< T >& rhs )
	{
		m_instance = rhs.m_instance ;
		if( m_instance )
		{
			m_refCount = rhs.m_refCount ;
			++*m_refCount ;
		}
	}

	void acquire( T* instance )
	{
		m_instance = instance ;
		if( m_instance ) {
			m_refCount = new unsigned ;
			*m_refCount = 1 ;
		} else {
			m_refCount = NULL ;
		}
	}
} ;

} //naemspace bogus

#endif
