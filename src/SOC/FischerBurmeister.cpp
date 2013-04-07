#include "FischerBurmeister.impl.hpp"

namespace bogus
{

  template struct FBBaseFunction< 2u, double > ;
  template struct FBBaseFunction< 3u, double > ;
  template class FischerBurmeister< 2u, double, false > ;
  template class FischerBurmeister< 3u, double, false > ;
  template class FischerBurmeister< 2u, double, true  > ;
  template class FischerBurmeister< 3u, double, true  > ;

}

