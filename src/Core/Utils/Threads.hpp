#ifndef BOGUS_UTILS_THREADS_HPP
#define BOGUS_UTILS_THREADS_HPP

#include "../Block/Constants.hpp"

#ifndef BOGUS_DONT_PARALLELIZE
#include <omp.h>
#endif

namespace bogus {

#ifdef BOGUS_DONT_PARALLELIZE
	struct WithMaxThreads {
		int nThreads() const { return 1 ; }
	} ;

#else
	struct WithMaxThreads {

		WithMaxThreads( int maxThreads )
			: m_prevMaxThreads( omp_get_max_threads() )
			,  m_newMaxThreads( maxThreads == 0 ? m_prevMaxThreads : maxThreads )
		{
			omp_set_num_threads( m_newMaxThreads ) ;
		}

		~WithMaxThreads() {
			omp_set_num_threads( m_prevMaxThreads ) ;
		}


		int nThreads() const { return m_newMaxThreads ; }

	private:
		const int m_prevMaxThreads ;
		const int m_newMaxThreads ;

	};
#endif

} // bogus

#endif
