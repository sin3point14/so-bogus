/*! 
  \file BlockSolvers.h 
  \brief High level documentation for the Core/BlockSolvers module
*/

namespace bogus {

/*! 

\page block_solvers BlockSolvers
\tableofcontents

\section block_basics Basics

To use the library, 
\code
#include <bogus/Core/BlockSolvers.impl.hpp>
\endcode

The \ref block_solvers module is a collection of solvers operating on \ref block matrices.

At the moment, those solvers are:
 - \ref block_solvers_cg
 - \ref block_solvers_gs

\section block_solvers_cg Conjugate Gradient

A few conjugate gradient variants are available through the ConjugateGradient class, as well as a few naive preconditioners.
	
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
  bogus::ConjugateGradient< Mat > cg( sbm ) ;
  cg.solve( rhs, res ) ;
  //or 
  cg.solve_BiCG( rhs, res ) ;
  //or
  cg.solve_BiCGSTAB( rhs, res ) ;
\endcode 

If we wanted to use a preconditioner
\code
  bogus::ConjugateGradient< Mat, bogus::DiagonalPreconditioner > pcg( sbm ) ;
  //or
  bogus::ConjugateGradient< Mat, bogus::DiagonalLDLTPreconditioner > ldltcg( sbm ) ;

  pcg.solve( rhs, res ) ;
\endcode 

... or sparse matrices

\code
  typedef Eigen::SparseMatrix< double > SparseBlock ;
  typedef bogus::SparseBlockMatrix< SparseBlock > SparseMat ;
  SparseMat ssbm ;
  ssbm.cloneStructure( sbm ) ;
  ssbm.block(0) =  sbm.block(0).sparseView() ;

  // DiagonalLDLTPreconditioner on SparseMatrix blocks requires Eigen 3.1+
  bogus::ConjugateGradient< SparseMat, bogus::DiagonalLDLTPreconditioner > sldltcg( ssbm ) ;
  err = sldltcg.solve( rhs, res ) ;
  
\endcode 

\section block_solvers_gs Projected Gauss Seidel

Constrained linear systems can be solved using the GaussSeidel class,
as long as a corresponding NSLaw is available. SOCLaw, from the \ref soc module, is an example of such NSlaw.

Example code for 3D Coulomb friction (requires the \ref soc module):
\code

//Construct the dual matrix
typedef bogus::SparseBlockMatrix< Eigen::Matrix3d, bougs::flags::SYMMETRIC | bogus::flags::COMPRESSED > WType ;
WType W = primal.H * primal.MInvHt ;

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

