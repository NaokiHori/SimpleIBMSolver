#include <string.h>
#include <math.h>
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "temperature.h"
#include "arrays/fluid.h"
#include "arrays/temperature.h"



/**
 * @brief comute right-hand-side of Runge-Kutta scheme
 * @param[in   ] domain      : information related to domain decomposition and size
 * @param[in   ] rkstep      : Runge-Kutta step
 * @param[in   ] fluid       : velocity
 * @param[inout] temperature : temperature field (in), RK source terms (inout)
 * @return                   : error code
 */
int temperature_compute_rhs(const domain_t * restrict domain, const int rkstep, const fluid_t * restrict fluid, temperature_t * restrict temperature){
  const int implicitx = config.get_bool("implicitx");
  const int implicity = config.get_bool("implicity");
  const double Ra = config.get_double("Ra");
  const double Pr = config.get_double("Pr");
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict dxf = domain->dxf;
  const double * restrict dxc = domain->dxc;
  const double dy = domain->dy;
  const double * restrict ux = fluid->ux;
  const double * restrict uy = fluid->uy;
  const double * restrict temp = temperature->temp;
  double * restrict srctempa = temperature->srctempa;
  double * restrict srctempb = temperature->srctempb;
  double * restrict srctempg = temperature->srctempg;
  /* ! previous k-step source term of temp is copied ! 3 ! */
  if(rkstep != 0){
    memcpy(srctempb, srctempa, SRCTEMPA_SIZE_0 * SRCTEMPA_SIZE_1 * sizeof(double));
  }
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      /* ! gradient of T in x ! 2 ! */
      double dtdx_xm = (-TEMP(i-1, j  )+TEMP(i  , j  ))/DXC(i  );
      double dtdx_xp = (-TEMP(i  , j  )+TEMP(i+1, j  ))/DXC(i+1);
      /* ! gradient of T in y ! 2 ! */
      double dtdy_ym = (-TEMP(i  , j-1)+TEMP(i  , j  ))/dy;
      double dtdy_yp = (-TEMP(i  , j  )+TEMP(i  , j+1))/dy;
      /* advection */
      /* ! T is transported by ux ! 10 ! */
      double advx;
      {
        double c_xm = DXC(i  )/(2.*DXF(i));
        double c_xp = DXC(i+1)/(2.*DXF(i));
        double ux_xm = UX(i  , j  );
        double ux_xp = UX(i+1, j  );
        advx =
          -c_xm*ux_xm*dtdx_xm
          -c_xp*ux_xp*dtdx_xp;
      }
      /* ! T is transported by uy ! 8 ! */
      double advy;
      {
        double uy_ym = UY(i  , j  );
        double uy_yp = UY(i  , j+1);
        advy =
          -0.5*uy_ym*dtdy_ym
          -0.5*uy_yp*dtdy_yp;
      }
      /* diffusion */
      /* ! T is diffused in x ! 5 ! */
      double difx;
      {
        const double prefactor = 1./sqrt(Pr)/sqrt(Ra);
        difx = prefactor*(dtdx_xp-dtdx_xm)/DXF(i);
      }
      /* ! T is diffused in y ! 5 ! */
      double dify;
      {
        const double prefactor = 1./sqrt(Pr)/sqrt(Ra);
        dify = prefactor*(dtdy_yp-dtdy_ym)/dy;
      }
      /* summation */
      /* ! summation of explicit terms ! 5 ! */
      SRCTEMPA(i, j) =
        +advx
        +advy
        +(1.-1.*implicitx)*difx
        +(1.-1.*implicity)*dify;
      /* ! summation of implicit terms ! 3 ! */
      SRCTEMPG(i, j) =
        +(  +1.*implicitx)*difx
        +(  +1.*implicity)*dify;
    }
  }
  return 0;
}

