#include "FischerBurmeister.impl.hpp"

namespace bogus
{

  template struct FBBaseFunction< 2u, double > ;
  template struct FBBaseFunction< 3u, double > ;
  template struct FischerBurmeister< 2u, double, false > ;
  template struct FischerBurmeister< 3u, double, false > ;
  template struct FischerBurmeister< 2u, double, true  > ;
  template struct FischerBurmeister< 3u, double, true  > ;

}

