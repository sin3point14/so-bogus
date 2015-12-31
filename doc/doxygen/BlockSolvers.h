/*!
  \file BlockSolvers.h
  \brief High level documentation for the Core/BlockSolvers module
*/

namespace bogus {

/*!

\page block_solvers BlockSolvers
\tableofcontents

\note This module is released under the terms of the <a href="http://mozilla.org/MPL/2.0/">Mozilla Public License version 2.0</a>

\section block_basics Basics

To use this library,
\code
#include <bogus/Core/BlockSolvers.impl.hpp>
// or if you only plan to use a specific part
#include <bogus/Core/BlockSolvers/Krylov.impl.hpp>
#include <bogus/Core/BlockSolvers/GaussSeidel.impl.hpp>
\endcode


The \ref block_solvers module is a collection of solvers operating on \ref block matrices.

At the moment, those solvers are:
 - \ref block_solvers_is
 - \ref block_solvers_gs
 - \ref block_solvers_pg

\section block_solvers_is Iterative Linear Solvers

A few Krylov methods are available through the Krylov class, as well as some naive preconditioners.
For the full list of available solvers, see the krylov::Method enum or the krylov::solvers namespace.

Here is some code solving a very simple system without preconditioning:
\code
  // Building a block matrix
  typedef bogus::SparseBlockMatrix< Eigen::Matrix3d > Mat ;
  Mat sbm ;
  sbm.setRows( 1, 3 ) ;
  sbm.setCols( 1, 3 ) ;
  sbm.insertBack(0,0) = Eigen::Vector3d::Constant( 2 ).asDiagonal() ;
  sbm.finalize() ;

  // Right-hand-side and result vectors
  Eigen::Vector3d rhs, res ;
  rhs.setOnes( ) ;

  // Solving
  bogus::Krylov< Mat > cg( sbm ) ;
  cg.solve_CG( rhs, res ) ;
  //or
  cg.solve_GMRES( rhs, res ) ;
  //or
  cg.solve( rhs, res, bogus::krylov::CGS ) ;
\endcode

If we wanted to use a preconditioner
\code
  bogus::Krylov< Mat, bogus::DiagonalPreconditioner > pcg( sbm ) ;
  //or
  bogus::Krylov< Mat, bogus::DiagonalLDLTPreconditioner > ldltcg( sbm ) ;

  pcg.solve_CG( rhs, res ) ;
\endcode

... or sparse matrices

\code
  typedef Eigen::SparseMatrix< double > SparseBlock ;
  typedef bogus::SparseBlockMatrix< SparseBlock > SparseMat ;
  SparseMat ssbm ;
  ssbm.cloneStructure( sbm ) ;
  ssbm.block(0) =  sbm.block(0).sparseView() ;

  // DiagonalLDLTPreconditioner on SparseMatrix blocks requires Eigen 3.1+
  bogus::Krylov< SparseMat, bogus::DiagonalLDLTPreconditioner > sldltcg( ssbm ) ;
  err = sldltcg.solve_CG( rhs, res ) ;

\endcode

Iterative linear solvers may also be used in a matrix-free version -- that is, 
without explictely computing the system matrix. Suppose that we want to solve
\f$ J J^T x = b \f$, this can be done as

\code
  typedef bogus::SparseBlockMatrix< BlockType > JType ;
  typedef bogus::Product< JType, bogus::Transpose< JType > > Prod ; // Or use c++11 to infer the correct type
  Prod W = J * J.transpose() ; 
  
  bogus::Krylov< Prod >( W ).solve( b, x ) ;
\endcode

Alternatively, a Krylov object may be converted to a particular method-object which
inherits from LinearSolverBase.
This means it can be used as a BlockType of a SparseBlockMatrix, and can offer
more configuration options, such as setting the 'restart' option for the
krylov::GMRES method.

\code
  typedef bogus::Krylov< Mat, bogus::DiagonalPreconditioner > KrylovType ;
  KrylovType krylov( sbm ) ;

  krylov.asGMRES().setRestart( 10 ).solve( rhs, res ) ;

  // Creating a SparseBlockMatrix of GMRES objects

  typedef typename KrylovType::GMRESType GMRES ;

  bogus::SparseBlockMatrix< GMRES > gmresMat ;
  // [..] Set rows, etc
  gmResmat.block( 0 ) = krylov.asGMRES() ;
  // [..]
  //
  // This performs no convergence check, and is probably a very bad idea
  res = gmresMat * rhs ;

\endcode

\section block_solvers_gs Projected Gauss Seidel

Constrained linear systems can be solved using the GaussSeidel class,
as long as a corresponding NSLaw is available. SOCLaw, from the \ref soc module, is an example of such NSlaw.

Example code for 3D Coulomb friction (requires the \ref soc module):
\code

//Construct the dual matrix
typedef bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::SYMMETRIC > WType ;
WType W = primal.H * ( primal.MInv * primal.H.transpose() ) ;

//Friction coefficients
const double * mu = { .... } ;

bogus::GaussSeidel< WType > gs( W ) ;
// Optional: connect our function to the callback, so we can monitor the GaussSeidel's convergence
gs.callback().connect( &ackCurrentResidual );

Eigen::VectorXd x( W.rows() ) ;
double res = gs.solve( bogus::Coulomb3D( n, mu ), b, x ) ;


\endcode

\section block_solvers_pg Projected Gradient

A projected gradient algorithm is also implemented in the ProjectedGradient class.
Its interface is very similar to that of the above GaussSeidel. A few variants of the
algorithm, such as the Nesterov \cite Nesterov1983 acceleration, are implemented; they can be selected using the ProjectedGradient::setDefaultVariant()
method or using a template parameter. See \ref projected_gradient::Variant for more information.

\code
bogus::ProjectedGradient< WType > pg( W ) ;
res = pg.solve( bogus::SOC3D( n, mu ), b, x ) ;

// ... or explicitely chose a variant

res = pg.solve< bogus::projected_gradient::APGD >( bogus::SOC3D( n, mu ), b, x ) ;
\endcode

\note This algorithm can only be used to solve constrained quadratic optimization problems.
Coulomb friction does not belong to this class, but Linear and Cone Complementarity problems do.

Once again, this algrithms can be used in a matrix-free fashion
\code
  typedef bogus::SparseBlockMatrix< BlockType > JType ;
  typedef bogus::Product< JType, bogus::Transpose< JType > > Prod ; // Or use c++11 to infer the correct type
  Prod W = J * J.transpose() ; 
  
  bogus::ProjectedGradient< Prod > pg( W ) ;
  res = pg.solve( bogus::SOC3D( n, mu ), b, x ) ;
\endcode
*/

}

