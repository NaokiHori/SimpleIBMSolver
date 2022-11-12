#include <math.h>
#include "config.h"
#include "common.h"
#include "domain.h"
#include "fluid.h"
#include "arrays/fluid.h"


#if NDIMS == 2

/**
 * @brief update pressure using scalar potential \psi
 * @param[in   ] domain : information related to MPI domain decomposition
 * @param[in   ] rkstep : Runge-Kutta step
 * @param[in   ] dt     : time step size
 * @param[inout] fluid  : scalar potential \psi (in), pressure (out)
 * @return              : error code
 */
int fluid_update_pressure(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict psi = fluid->psi;
  double * restrict p = fluid->p;
  /* ! add correction ! 5 ! */
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      P(i, j) += PSI(i, j);
    }
  }
  /* ! correction for implicit treatment in x ! 18 ! */
  if(config.get_bool("implicitx")){
    const double Ra = config.get_double("Ra");
    const double Pr = config.get_double("Pr");
    const double Re = sqrt(Ra)/sqrt(Pr);
    const double gamma = RKCOEFS[rkstep].gamma;
    const double * restrict dxf = domain->dxf;
    const double * restrict dxc = domain->dxc;
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        double dpsidx_xm = (-PSI(i-1, j  )+PSI(i  , j  ))/DXC(i  );
        double dpsidx_xp = (-PSI(i  , j  )+PSI(i+1, j  ))/DXC(i+1);
        P(i, j) +=
          -(gamma*dt)/(2.*Re)*(
              +(dpsidx_xp-dpsidx_xm)/DXF(i)
        );
      }
    }
  }
  /* ! correction for implicit treatment in y ! 17 ! */
  if(config.get_bool("implicity")){
    const double Ra = config.get_double("Ra");
    const double Pr = config.get_double("Pr");
    const double Re = sqrt(Ra)/sqrt(Pr);
    const double gamma = RKCOEFS[rkstep].gamma;
    const double dy = domain->dy;
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        double dpsidy_ym = (-PSI(i  , j-1)+PSI(i  , j  ))/dy;
        double dpsidy_yp = (-PSI(i  , j  )+PSI(i  , j+1))/dy;
        P(i, j) +=
          -(gamma*dt)/(2.*Re)*(
              +(dpsidy_yp-dpsidy_ym)/dy
        );
      }
    }
  }
  /* ! boundary and halo values are updated ! 1 ! */
  fluid_update_boundaries_p(domain, p);
  return 0;
}

#else // NDIMS == 3

/**
 * @brief update pressure using scalar potential \psi
 * @param[in   ] domain : information related to MPI domain decomposition
 * @param[in   ] rkstep : Runge-Kutta step
 * @param[in   ] dt     : time step size
 * @param[inout] fluid  : scalar potential \psi (in), pressure (out)
 * @return              : error code
 */
int fluid_update_pressure(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict psi = fluid->psi;
  double * restrict p = fluid->p;
  /* ! add correction ! 7 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        P(i, j, k) += PSI(i, j, k);
      }
    }
  }
  /* ! correction for implicit treatment in x ! 20 ! */
  if(config.get_bool("implicitx")){
    const double Ra = config.get_double("Ra");
    const double Pr = config.get_double("Pr");
    const double Re = sqrt(Ra)/sqrt(Pr);
    const double gamma = RKCOEFS[rkstep].gamma;
    const double * restrict dxf = domain->dxf;
    const double * restrict dxc = domain->dxc;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          double dpsidx_xm = (-PSI(i-1, j  , k  )+PSI(i  , j  , k  ))/DXC(i  );
          double dpsidx_xp = (-PSI(i  , j  , k  )+PSI(i+1, j  , k  ))/DXC(i+1);
          P(i, j, k) +=
            -(gamma*dt)/(2.*Re)*(
                +(dpsidx_xp-dpsidx_xm)/DXF(i)
          );
        }
      }
    }
  }
  /* ! correction for implicit treatment in y ! 19 ! */
  if(config.get_bool("implicity")){
    const double Ra = config.get_double("Ra");
    const double Pr = config.get_double("Pr");
    const double Re = sqrt(Ra)/sqrt(Pr);
    const double gamma = RKCOEFS[rkstep].gamma;
    const double dy = domain->dy;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          double dpsidy_ym = (-PSI(i  , j-1, k  )+PSI(i  , j  , k  ))/dy;
          double dpsidy_yp = (-PSI(i  , j  , k  )+PSI(i  , j+1, k  ))/dy;
          P(i, j, k) +=
            -(gamma*dt)/(2.*Re)*(
                +(dpsidy_yp-dpsidy_ym)/dy
          );
        }
      }
    }
  }
  /* ! correction for implicit treatment in z ! 19 ! */
  if(config.get_bool("implicitz")){
    const double Ra = config.get_double("Ra");
    const double Pr = config.get_double("Pr");
    const double Re = sqrt(Ra)/sqrt(Pr);
    const double gamma = RKCOEFS[rkstep].gamma;
    const double dz = domain->dz;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          double dpsidz_zm = (-PSI(i  , j  , k-1)+PSI(i  , j  , k  ))/dz;
          double dpsidz_zp = (-PSI(i  , j  , k  )+PSI(i  , j  , k+1))/dz;
          P(i, j, k) +=
            -(gamma*dt)/(2.*Re)*(
                +(dpsidz_zp-dpsidz_zm)/dz
          );
        }
      }
    }
  }
  /* ! boundary and halo values are updated ! 1 ! */
  fluid_update_boundaries_p(domain, p);
  return 0;
}

#endif // NDIMS
