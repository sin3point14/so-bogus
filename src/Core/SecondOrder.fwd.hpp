#ifndef BOGUS_SECOND_ORDER_FWD_HPP
#define BOGUS_SECOND_ORDER_FWD_HPP

namespace bogus
{

namespace local_soc_solver
{
enum Strategy
{
	PureNewton,
	PureEnumerative,
	Hybrid,
	RevHybrid
} ;
}

template < typename LocalMatrixType, bool DeSaxceCOV,
		   local_soc_solver::Strategy Strat = local_soc_solver::Hybrid  >
class SOCLaw ;

}

#endif
