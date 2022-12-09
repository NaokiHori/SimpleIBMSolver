#include "domain.h"
#include "temperature.h"
#include "arrays/temperature.h"



/**
 * @brief update boundary values of temperature
 * @param[in   ] domain : information about domain decomposition and size
 * @param[inout] temp   : temperature
 * @return              : error code
 */
int temperature_update_boundaries_temp(const domain_t * restrict domain, double * restrict temp){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  /* ! update halo values ! 1 ! */
  domain_communicate_halo_p_like(domain, temp);
  /* ! set boundary values of temp ! 8 ! */
  // temp_xm - temp_xp should be 1,
  //   since governing equations are normalised based on the assumption
  const double temp_xm = +0.5;
  const double temp_xp = -0.5;
  for(int j = 1; j <= jsize; j++){
    TEMP(      0, j) = temp_xm;
    TEMP(isize+1, j) = temp_xp;
  }
  return 0;
}

