/*! 
  \file BlockSolvers.h 
  \brief High level documentation for the Core/BlockSolvers module
*/

namespace bogus {

/*! 

\page block_solvers BlockSolvers
\tableofcontents

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

\section block_solvers_is Iterative Linear Solvers

A few Krylov methods are available through the Krylov class, as well as some naive preconditioners.
	
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

*/

}

