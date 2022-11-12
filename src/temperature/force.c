#include "config.h"
#include "domain.h"
#include "temperature.h"
#include "arrays/temperature.h"


#if NDIMS == 2

/**
 * @brief compute buoyancy force
 * @param[in   ] domain      : information related to MPI domain decomposition
 * @param[inout] temperature : temperature (in), buoyancy force (out)
 * @return                   : error code
 */
int temperature_compute_force(const domain_t * restrict domain, temperature_t * restrict temperature){
  const bool add_buoyancy = config.get_bool("add_buoyancy");
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict temp = temperature->temp;
  double * restrict tempforcex = temperature->tempforcex;
  if(add_buoyancy){
    /* ! buoyancy force acting only to wall-normal (x) direction ! 7 ! */
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        double temp_xm = TEMP(i-1, j  );
        double temp_xp = TEMP(i  , j  );
        TEMPFORCEX(i, j) = 0.5*temp_xm+0.5*temp_xp;
      }
    }
  }else{
    /* ! no forcing ! 5 ! */
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        TEMPFORCEX(i, j) = 0.;
      }
    }
  }
  return 0;
}

#else // NDIMS == 3

/**
 * @brief compute buoyancy force
 * @param[in   ] domain      : information related to MPI domain decomposition
 * @param[inout] temperature : temperature (in), buoyancy force (out)
 * @return                   : error code
 */
int temperature_compute_force(const domain_t * restrict domain, temperature_t * restrict temperature){
  const bool add_buoyancy = config.get_bool("add_buoyancy");
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict temp = temperature->temp;
  double * restrict tempforcex = temperature->tempforcex;
  if(add_buoyancy){
    /* ! buoyancy force acting only to wall-normal (x) direction ! 9 ! */
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 2; i <= isize; i++){
          double temp_xm = TEMP(i-1, j  , k  );
          double temp_xp = TEMP(i  , j  , k  );
          TEMPFORCEX(i, j, k) = 0.5*temp_xm+0.5*temp_xp;
        }
      }
    }
  }else{
    /* ! no forcing ! 7 ! */
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 2; i <= isize; i++){
          TEMPFORCEX(i, j, k) = 0.;
        }
      }
    }
  }
  return 0;
}

#endif // NDIMS
