#ifndef BOGUS_MECHE_INTERFACE_HPP
#define BOGUS_MECHE_INTERFACE_HPP

namespace bogus
{

class MecheFrictionProblem
{
public:
	MecheFrictionProblem()
		: m_data( 0 )
	{}
	~MecheFrictionProblem() ;

	void fromPrimal (
			unsigned int NObj, //!< number of subsystems
			const unsigned int * ndof, //!< array of size \a NObj, the number of degree of freedom of each subsystem
			const double *const * MassMat, //!< array of pointers to the mass matrix of each subsystem
			const double * f_in, //!< the constant term in \f$ M v + f= {}^t \! H r \f$
			unsigned int n_in, //!< number of contact points
			const double * mu_in, //!< array of size \a n giving the friction coeffs
			const double * E_in, //!< array of size \f$ n \times d \times d \f$ giving the \a n normals followed by the \a n tangent vectors (and by again \a n tangent vectors if \a d is 3). Said otherwise, \a E is a \f$ (nd) \times d \f$ matrix, stored column-major, formed by \a n blocks of size \f$ d \times d \f$ with each block being an orthogonal matrix (the transition matrix from the world space coordinates \f$ (x_1, x_2, x_3) \f$ to the local coordinates \f$ (x_N, x_{T1}, x_{T2}) \f$
			const double * w_in, //!< array of size \a nd, the constant term in \f$ u = H v + w \f$
			const int * const ObjA, //!< array of size \a n, the first object involved in the \a i-th contact (must be an internal object) (counted from 0)
			const int * const ObjB, //!< array of size \a n, the second object involved in the \a i-th contact (-1 for an external object) (counted from 0)
			const double *const HA[], //!< array of size \a n, containing pointers to a dense, colum-major matrix of size <c> d*ndof[ObjA[i]] </c> corresponding to the H-matrix of <c> ObjA[i] </c>
			const double *const HB[] //!< array of size \a n, containing pointers to a dense, colum-major matrix of size <c> d*ndof[ObjA[i]] </c> corresponding to the H-matrix of <c> ObjB[i] </c> (\c NULL for an external object)
			);

	double solve(
			double * r, //!< length \a nd : initialization for \a r (in world space coordinates) + used to return computed r
			double * v, //!< length \a m: to return computed v ( or NULL if not needed )
			bool deterministic = false,       //!< Whether the Gauss-Seidel should be eterministic
			double tol = 0.,                  //!< Gauss-Seidel tolerance. 0. means GS's default
			unsigned maxIters = 0             //!< Max number of iterations. 0 means GS's default
			);

protected:
	struct Data  ;

	Data* m_data ;

} ;

}

#endif
