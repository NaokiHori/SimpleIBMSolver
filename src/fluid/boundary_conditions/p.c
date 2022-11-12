#include "domain.h"
#include "fluid.h"
#include "arrays/fluid.h"


#if NDIMS == 2

/**
 * @brief update boundary values of pressure and scalar potential
 * @param[in   ] domain : information about domain decomposition and size
 * @param[inout] p      : pressure or scalar potential
 * @return              : error code
 */
int fluid_update_boundaries_p(const domain_t * restrict domain, double * restrict p){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  /* ! update halo values ! 1 ! */
  domain_communicate_halo_p_like(domain, p);
  /* ! set boundary values ! 4 ! */
  for(int j = 1; j <= jsize; j++){
    P(      0, j) = P(    1, j); // Neumann
    P(isize+1, j) = P(isize, j); // Neumann
  }
  return 0;
}

#else // NDIMS == 3

/**
 * @brief update boundary values of pressure and scalar potential
 * @param[in   ] domain : information about domain decomposition and size
 * @param[inout] p      : pressure or scalar potential
 * @return              : error code
 */
int fluid_update_boundaries_p(const domain_t * restrict domain, double * restrict p){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  /* ! update halo values ! 1 ! */
  domain_communicate_halo_p_like(domain, p);
  /* ! set boundary values ! 6 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      P(      0, j, k) = P(    1, j, k); // Neumann
      P(isize+1, j, k) = P(isize, j, k); // Neumann
    }
  }
  return 0;
}

#endif // NDIMS
