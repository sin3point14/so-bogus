/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef BOGUS_CPP_TOOLS_HPP
#define BOGUS_CPP_TOOLS_HPP

#include <iostream>

namespace bogus
{

// Swap template parameters if DoSwap is true

template < bool DoSwap, typename First_, typename Second_ >
struct TypeSwapIf
{
    typedef First_  First  ;
    typedef Second_ Second ;
} ;

template < typename First_, typename Second_ >
struct TypeSwapIf< true, First_, Second_ >
{
    typedef First_  Second  ;
    typedef Second_ First ;
} ;

template < bool DoSwap, int First_, int Second_ >
struct SwapIf
{
    enum { First = First_, Second = Second_  } ;
} ;

template < int First_, int Second_ >
struct SwapIf< true, First_, Second_ >
{
    enum { First = Second_, Second = First_  } ;
} ;


// Enable if (for SFINAE use )

template < bool Condition, typename ReturnType_ = void >
struct EnableIf
{
} ;

template < typename ReturnType_ >
struct EnableIf< true, ReturnType_ >
{
    typedef ReturnType_ ReturnType ;
} ;

template < bool Condition, typename ReturnType_ = void >
struct DisableIf
{
} ;

template < typename ReturnType_ >
struct DisableIf< false, ReturnType_ >
{
    typedef ReturnType_ ReturnType ;
} ;

template < typename BaseType >
struct HasReturnType
{
private:
    enum { True = 1, False = 2 } ;
    typedef char  TrueType[  True ] ;
    typedef char FalseType[ False ] ;

    template< typename T >
    static const  TrueType& check( const typename T::ReturnType* ) ;
    template< typename >
    static const FalseType& check( ... ) ;
public:
    enum { Value = ( True == sizeof( check< BaseType >( 0 ) ) ) } ;
} ;


// Static assertions

template < bool Assertion >
struct StaticAssert
{
    enum {
        BLOCKS_MUST_BE_SQUARE_OR_HAVE_DYNAMIC_DIMENSIONS,
        BLOCKS_MUST_HAVE_FIXED_DIMENSIONS,
        MATRICES_ORDERING_IS_INCONSISTENT,
        TRANSPOSE_OF_FACTORIZATION_MAKES_NO_SENSE_IN_THIS_CONTEXT,
        TRANSPOSE_IS_NOT_DEFINED_FOR_THIS_BLOCK_TYPE
    } ;
} ;

template < >
struct StaticAssert< false >
{
} ;

#define BOGUS_STATIC_ASSERT( test, message ) (void) StaticAssert< test >::message

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
    void reset( T * instance = 0 )
    {
        release() ;
        acquire( instance ) ;
    }

    NaiveSharedPtr( const NaiveSharedPtr< T >& rhs )
    {
        *this = rhs ;
    }

    const NaiveSharedPtr& operator=( const NaiveSharedPtr< T >& rhs )
    {
        if( this != &rhs )
        {
            release() ;
            m_instance = rhs.m_instance ;
            if( m_instance )
            {
                m_refCount = rhs.m_refCount ;
                ++*m_refCount ;
            }
        }
        return *this ;
    }

    ~NaiveSharedPtr()
    {
        release() ;
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

} ;

#ifndef BOGUS_SHARED_PTR
#define BOGUS_SHARED_PTR( Type, Name ) NaiveSharedPtr< Type > Name
#endif

} //namespace bogus


#endif



