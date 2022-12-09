
#include "common.h"
#include "domain.h"
#include "fluid.h"
#include "arrays/fluid.h"
#include "internal.h"


/**
 * @brief correct uz using scalar potential \psi
 * @param[in   ] domain    : information about domain decomposition and size
 * @param[in   ] prefactor : pre-factor in front of grad \psi
 * @param[inout] fluid     : scalar potential \psi (in), uz (out)
 * @return                 : error code
 */
int fluid_correct_velocity_uz(const domain_t * restrict domain, const double prefactor, fluid_t * restrict fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double dz  = domain->dz;
  const double * restrict psi = fluid->psi;
  double * restrict uz = fluid->uz;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        /* ! correct ! 3 ! */
        double psi_zm = PSI(i  , j  , k-1);
        double psi_zp = PSI(i  , j  , k  );
        UZ(i, j, k) -= prefactor*(psi_zp-psi_zm)/dz;
      }
    }
  }
  /* ! update boundary and halo values ! 1 ! */
  fluid_update_boundaries_uz(domain, uz);
  return 0;
}

