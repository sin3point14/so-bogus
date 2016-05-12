/*
 * Any copyright is dedicated to the Public Domain.
 * http://creativecommons.org/publicdomain/zero/1.0/
*/

#ifndef D6_TESTS_RESIDUAL_INFO_HH
#define D6_TESTS_RESIDUAL_INFO_HH

#include "Core/Utils/Signal.hpp"

#include <string>

class ResidualInfo {

public:
	explicit ResidualInfo( bool verbose = false )
		: m_verbose( verbose )
	{ }

	void setVerbose( bool verbose ) {
		m_verbose = verbose ;
	}

	void ack( unsigned iter, double err ) {
		if( m_verbose ) {
			std::cout << m_meth << ": \t" << iter << "\t ==> " << err << std::endl ;
		}
	}

	void setMethodName( const std::string& meth ) {
		m_meth = meth ;
	}

	void bindTo( bogus::Signal<unsigned, double> &signal ) {
		signal.connect( *this, &ResidualInfo::ack );
	}

private:

	bool m_verbose ;
	std::string m_meth ;

};


#endif
