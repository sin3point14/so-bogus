/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_LOCK_HPP
#define BOGUS_LOCK_HPP

#ifdef BOGUS_DONT_PARALLELIZE
class Lock {
public:
	int  *for_abi_compat ;

	struct Guard {
		explicit Guard( Lock& ) {}
	} ;
};
#else

#include <omp.h>

class Lock {
public:
	struct Guard {
		explicit Guard( const Lock& lock )
			: lockPtr( lock.ptr() )
		{
			omp_set_lock( lockPtr ) ;
		}

		~Guard()
		{
			omp_unset_lock( lockPtr ) ;
		}
	private:
		Guard(const Guard &guard) ;
		Guard& operator=(const Guard &guard) ;

		omp_lock_t * lockPtr ;
	} ;

	Lock()
		: m_lock( new omp_lock_t )
	{
		omp_init_lock( m_lock ) ;
	}

	Lock( const Lock& )
		: m_lock( new omp_lock_t )
	{
		omp_init_lock( m_lock ) ;
	}

	Lock& operator=( const Lock& )
	{
		return *this ;
	}

	~Lock()
	{
		delete m_lock ;
	}

	omp_lock_t* ptr() const { return m_lock ; }

private:
	omp_lock_t* m_lock ;
};

#endif

#endif
