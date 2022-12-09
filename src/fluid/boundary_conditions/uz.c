
#include "domain.h"
#include "fluid.h"
#include "arrays/fluid.h"


/**
 * @brief update boundary values of z velocity
 * @param[in   ] domain : information about domain decomposition and size
 * @param[inout] uz     : z velocity
 * @return              : error code
 */
int fluid_update_boundaries_uz(const domain_t * restrict domain, double * restrict uz){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  /* ! update halo values ! 1 ! */
  domain_communicate_halo_uz_like(domain, uz);
  /* ! set boundary values ! 6 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      UZ(      0, j, k) = 0.; // no-slip
      UZ(isize+1, j, k) = 0.; // no-slip
    }
  }
  return 0;
}

