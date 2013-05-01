/* This file is part of so-bogus, a block-sparse Gauss-Seidel solver
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef BOGUS_SIGNAL_HPP
#define BOGUS_SIGNAL_HPP

#include <list>

namespace bogus {

template < typename Derived >
struct SignalTraits
{ } ;

template< typename Arg1, typename Arg2 = void >
struct Signal ;

template< typename Derived >
class SignalBase
{
	typedef SignalTraits< Derived > Traits ;

public:
	virtual ~SignalBase()
	{
		disconnectAll();
	}

	void disconnectAll() ;

	void connect( typename Traits::Function::Type func ) ;

	template <typename T >
	void connect( T& object, typename Traits::template Method< T >::Type member_func ) ;

protected:
	typedef std::list< typename Traits::Callable* > Callables ;
	Callables  m_callees ;

} ;

template < typename Arg1, typename Arg2 >
struct SignalTraits< Signal< Arg1, Arg2 > >
{
	struct Callable
	{
		virtual ~Callable() {}
		virtual void call( Arg1, Arg2 ) = 0 ;
	} ;

	struct Function : public Callable
	{
		typedef  void (*Type)( Arg1, Arg2 ) ;
		Type func ;
		Function ( Type _func ) : func( _func ) {}
		virtual void call( Arg1 arg1, Arg2 arg2 ) { func( arg1, arg2 ) ; }
	} ;
	template< typename T >
	struct Method : public Callable
	{
		typedef  void (T::*Type)( Arg1, Arg2 ) ;
		T& obj ;
		Type func ;
		Method ( T& _obj, Type _func ) : obj( _obj ), func( _func ) {}
		virtual void call( Arg1 arg1, Arg2 arg2 ) { (obj.*func)( arg1, arg2 ) ; }
	} ;
} ;

template< typename Arg >
struct SignalTraits< Signal< Arg, void > >
{
	struct Callable
	{
		virtual ~Callable() {}
		virtual void call( Arg ) = 0 ;
	} ;

	struct Function : public Callable
	{
		typedef  void (*Type)( Arg ) ;
		Type func ;
		Function ( Type _func ) : func( _func ) {}
		virtual void call( Arg arg ) { func( arg ) ; }
	} ;
	template< typename T >
	struct Method : public Callable
	{
		typedef  void (T::*Type)( Arg ) ;
		T& obj ;
		Type func ;
		Method ( T& _obj, Type _func ) : obj( _obj ), func( _func ) {}
		virtual void call( Arg arg ) { (obj.*func)( arg ) ; }
	} ;

} ;

template< typename Arg1, typename Arg2 >
struct Signal : public SignalBase< Signal< Arg1, Arg2 > >
{
	void trigger( Arg1 arg1, Arg2 arg2 ) const
	{
		typedef SignalBase< Signal< Arg1, Arg2 > > Base ;

		for( typename Base::Callables::const_iterator it = this->m_callees.begin() ; it != this->m_callees.end() ; ++it )
		{ (*it)->call( arg1, arg2 ) ; }
	}

} ;

template< typename Arg >
struct Signal< Arg, void > : public SignalBase< Signal< Arg, void > >
{
	void trigger( Arg arg ) const
	{
		typedef SignalBase< Signal< Arg, void > > Base ;

		for( typename Base::Callables::const_iterator it = this->m_callees.begin() ; it != this->m_callees.end() ; ++it )
		{ (*it)->call( arg ) ; }
	}
} ;

template< typename Derived >
void SignalBase< Derived >::disconnectAll() {
	for( typename Callables::iterator it = m_callees.begin() ; it != m_callees.end() ; ++it )
	{
		delete *it ;
	}
	m_callees.clear() ;
}

template< typename Derived >
void SignalBase< Derived >::connect( typename Traits::Function::Type func )
{
	m_callees.push_back( new typename Traits::Function( func ) );
}

template< typename Derived >
template <typename T >
void SignalBase< Derived >::connect( T& object, typename Traits::template Method< T >::Type member_func )
{
	m_callees.push_back( new typename Traits::template Method< T >( object, member_func ) );
}


} // namespace bogus

#endif
