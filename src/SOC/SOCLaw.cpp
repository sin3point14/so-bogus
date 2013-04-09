#include "SOCLaw.impl.hpp"

namespace bogus
{

template class SOCLaw< Eigen::Matrix2d, false >  ;
template class SOCLaw< Eigen::Matrix3d, false >  ;
template class SOCLaw< Eigen::Matrix2d,  true >  ;
template class SOCLaw< Eigen::Matrix3d,  true, local_soc_solver::PureNewton >  ;
template class SOCLaw< Eigen::Matrix3d,  true, local_soc_solver::PureEnumerative >  ;
template class SOCLaw< Eigen::Matrix3d,  true, local_soc_solver::Hybrid >  ;
template class SOCLaw< Eigen::Matrix3d,  true, local_soc_solver::RevHybrid >  ;

}
