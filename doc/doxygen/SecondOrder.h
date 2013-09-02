/*! 
  \file SecondOrder.h 
  \brief High level documentation for the Core/SecondOrder module
*/

namespace bogus {

/*! 

\page soc Second Order
\tableofcontents

\note This module is released under the terms of the <a href="http://www.gnu.org/licenses/gpl-2.0.html">GNU General Public License version 2</a>

\section soc_basics Basics

To use the library, 
\code
#include <bogus/Extra/SecondOrder.impl.hpp>
\endcode

The \ref soc module provide tools for dealing with Second Order Cone complementarity problems. 
Combining it with the GaussSeidel implementation from the \ref block_solvers module yields efficient
solvers for Coulomb friction and \b SOCQP problems.

The only user-facing class of this module is SOCLaw, which can be used as the first argument of GaussSeidel::solve().
It basically implements the algorithms described in the appendices of \cite DDB11 . 

As such, the user interested in those implementations may want to look at 
 - FischerBurmeister and FBBaseFunction for the computation of the \b MFB function 
 - LocalSOCSolver specializations for the hybrid algorithm and the degree 4 polynomial computation.

Experimenting with these different algorithm for solving the 1-contact problem may be done by specifying the local_soc_solver::Strategy template argument of the SOCLaw class.


*/

}

