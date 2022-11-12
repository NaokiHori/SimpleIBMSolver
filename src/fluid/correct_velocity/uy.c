#include "common.h"
#include "domain.h"
#include "fluid.h"
#include "arrays/fluid.h"
#include "internal.h"


#if NDIMS == 2

/**
 * @brief correct uy using scalar potential \psi
 * @param[in   ] domain    : information about domain decomposition and size
 * @param[in   ] prefactor : pre-factor in front of grad \psi
 * @param[inout] fluid     : scalar potential \psi (in), uy (out)
 * @return                 : error code
 */
int fluid_correct_velocity_uy(const domain_t * restrict domain, const double prefactor, fluid_t * restrict fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double dy  = domain->dy;
  const double * restrict psi = fluid->psi;
  double * restrict uy = fluid->uy;
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      /* ! correct ! 3 ! */
      double psi_ym = PSI(i  , j-1);
      double psi_yp = PSI(i  , j  );
      UY(i, j) -= prefactor*(psi_yp-psi_ym)/dy;
    }
  }
  /* ! update boundary and halo values ! 1 ! */
  fluid_update_boundaries_uy(domain, uy);
  return 0;
}

#else // NDIMS == 3

/**
 * @brief correct uy using scalar potential \psi
 * @param[in   ] domain    : information about domain decomposition and size
 * @param[in   ] prefactor : pre-factor in front of grad \psi
 * @param[inout] fluid     : scalar potential \psi (in), uy (out)
 * @return                 : error code
 */
int fluid_correct_velocity_uy(const domain_t * restrict domain, const double prefactor, fluid_t * restrict fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double dy  = domain->dy;
  const double * restrict psi = fluid->psi;
  double * restrict uy = fluid->uy;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        /* ! correct ! 3 ! */
        double psi_ym = PSI(i  , j-1, k  );
        double psi_yp = PSI(i  , j  , k  );
        UY(i, j, k) -= prefactor*(psi_yp-psi_ym)/dy;
      }
    }
  }
  /* ! update boundary and halo values ! 1 ! */
  fluid_update_boundaries_uy(domain, uy);
  return 0;
}

#endif // NDIMS
