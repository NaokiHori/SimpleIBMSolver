#include "domain.h"
#include "fluid.h"
#include "arrays/fluid.h"



/**
 * @brief update boundary values of y velocity
 * @param[in   ] domain : information about domain decomposition and size
 * @param[inout] uy     : y velocity
 * @return              : error code
 */
int fluid_update_boundaries_uy(const domain_t * restrict domain, double * restrict uy){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  /* ! update halo values ! 1 ! */
  domain_communicate_halo_uy_like(domain, uy);
  /* ! set boundary values ! 4 ! */
  for(int j = 1; j <= jsize; j++){
    UY(      0, j) = 0.; // no-slip
    UY(isize+1, j) = 0.; // no-slip
  }
  return 0;
}

