/*
 * Any copyright is dedicated to the Public Domain.
 * http://creativecommons.org/publicdomain/zero/1.0/
*/

#ifndef BOGUS_FCLIB_SOLVER_HPP
#define BOGUS_FCLIB_SOLVER_HPP

#include <Core/Block.impl.hpp>
#include <Core/Block.io.hpp>
#include <Core/BlockSolvers/ProjectedGradient.impl.hpp>
#include <Core/BlockSolvers/GaussSeidel.impl.hpp>

#include <Core/Utils/Timer.hpp>

#include <Interfaces/FrictionProblem.hpp>

namespace  bogus {
namespace fclib {

//! Solver configuration
struct Options {
	enum Algorithm {
		GaussSeidel,
		ProjectedGradient
	};

	int maxThreads ;        //!< Maximum number of threads that the GS will use.

	unsigned maxIters ;     //!< Max number of iterations. 0 means GS's default
	unsigned cadouxIters ;  //!< If staticProblem is false and cadouxIters is greater than zero, use the Cadoux algorithm to solve the friction problem.

	double tolerance ;      //!< Solver tolerance
	bool useInfinityNorm ;  //!< Whether to use the infinity norm to evaluate the residual of the friction problem,

	Algorithm algorithm ;

	// Solver-specific options
	double gsRegularization ; //!< GS proximal regularization coefficient
	bool   gsColoring ;       //!< Use coloring for parallel GS; slower but deterministic

	projected_gradient::Variant pgVariant ;


	Options()
	: maxThreads(0), maxIters(100000), cadouxIters(0),
	  tolerance(1.e-16), useInfinityNorm( false ),
	  algorithm( GaussSeidel ),
	  gsRegularization( 0 ), gsColoring( false ),
	  pgVariant( projected_gradient::SPG )
	{
	}
};

//! Solver result
struct Stats {
	bool verbose ;

	// Last iteration info
	unsigned nIters ;
	double   error  ;
	double   time   ;

	// Connect there for info about each iteration
	Signal<unsigned, double, double> callback ;
	Timer    timer ;

	explicit Stats( bool v = false )
	    : verbose(v), nIters(0), error(-1), time(0)
	{}

	void reset() {
		timer.reset() ;
	}

	void update( unsigned it, double err )
	{
		nIters = it ;
		error  = err ;
		time   = timer.elapsed() ;
		if( verbose )
			std::cout << it << ":\t" << err << "\t" << time << std::endl ;
		callback.trigger( it, err, time );
	}
};

template< unsigned Dimension >
double solveDual( const DualFrictionProblem<Dimension> &dual,
                         const Options& options,
                         Eigen::VectorXd &r, Eigen::VectorXd &u,
                         Stats& stats )
{

	Signal< unsigned, double > callback ;
	callback.connect( stats, &Stats::update );
	stats.reset();

	double res = -1 ;

	if( options.algorithm == Options::GaussSeidel ) {
		typename bogus::DualFrictionProblem< Dimension >::GaussSeidelType gs ;

		gs.setTol( options.tolerance ) ;
		gs.useInfinityNorm( options.useInfinityNorm );
		gs.setMaxIters( options.maxIters ) ;
		gs.setAutoRegularization( options.gsRegularization ) ;
		gs.setMaxThreads( options.maxThreads );
		gs.coloring().update( options.gsColoring, dual.W );

		if( options.cadouxIters > 0 ) {
			res = dual.solveCadoux( gs, r.data(), options.cadouxIters, &callback ) ;
		} else {
			gs.callback().connect( callback );
			res = dual.solveWith( gs, r.data(), false ) ;
		}
	} else {
		typename bogus::DualFrictionProblem< Dimension >::ProjectedGradientType pg ;
		pg.setTol( options.tolerance ) ;
		pg.useInfinityNorm( options.useInfinityNorm );
		pg.setMaxIters( options.maxIters ) ;
		pg.setDefaultVariant( options.pgVariant );

		if( options.cadouxIters > 0 ) {
			res = dual.solveCadoux( pg, r.data(), options.cadouxIters, &callback ) ;
		} else {
			pg.callback().connect( callback );
			res = dual.solveWith( pg, r.data(), false ) ;
		}
	}


	u = dual.W * r + dual.b ;

	return res ;
}

} // ns fclib

} // ns bogus

#endif
