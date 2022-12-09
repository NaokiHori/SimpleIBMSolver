#include <string.h>
#include <math.h>
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "temperature.h"
#include "arrays/fluid.h"
#include "arrays/temperature.h"
#include "internal.h"



/**
 * @brief comute right-hand-side of Runge-Kutta scheme of ux
 * @param[in   ] domain      : information related to MPI domain decomposition
 * @param[in   ] rkstep      : Runge-Kutta step
 * @param[inout] fluid       : velocity and pressure (in), RK source terms (inout)
 * @param[in   ] temperature : buoyancy force (in)
 * @return                   : error code
 */
int compute_rhs_ux(const domain_t * restrict domain, const int rkstep, fluid_t * restrict fluid, const temperature_t * restrict temperature){
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
  const double * restrict p = fluid->p;
  const double * restrict tempforcex = temperature->tempforcex;
  double * restrict srcuxa = fluid->srcuxa;
  double * restrict srcuxb = fluid->srcuxb;
  double * restrict srcuxg = fluid->srcuxg;
  /* ! previous k-step source term is copied ! 3 ! */
  if(rkstep != 0){
    memcpy(srcuxb, srcuxa, SRCUXA_SIZE_0 * SRCUXA_SIZE_1 * sizeof(double));
  }
  for(int j = 1; j <= jsize; j++){
    for(int i = 2; i <= isize; i++){
      /* ! velocity-gradient tensor L_xx ! 2 ! */
      double duxdx_xm = (-UX(i-1, j  )+UX(i  , j  ))/DXF(i-1);
      double duxdx_xp = (-UX(i  , j  )+UX(i+1, j  ))/DXF(i  );
      /* ! velocity-gradient tensor L_xy ! 2 ! */
      double duxdy_ym = (-UX(i  , j-1)+UX(i  , j  ))/dy;
      double duxdy_yp = (-UX(i  , j  )+UX(i  , j+1))/dy;
      /* advection */
      /* ! transported by ux ! 10 ! */
      double advx;
      {
        double c_xm = DXF(i-1)/(2.*DXC(i));
        double c_xp = DXF(i  )/(2.*DXC(i));
        double ux_xm = 0.5*UX(i-1, j  )+0.5*UX(i  , j  );
        double ux_xp = 0.5*UX(i  , j  )+0.5*UX(i+1, j  );
        advx =
          -c_xm*ux_xm*duxdx_xm
          -c_xp*ux_xp*duxdx_xp;
      }
      /* ! transported by uy ! 10 ! */
      double advy;
      {
        double c_xm = DXF(i-1)/(2.*DXC(i));
        double c_xp = DXF(i  )/(2.*DXC(i));
        double uy_ym = c_xm*UY(i-1, j  )+c_xp*UY(i  , j  );
        double uy_yp = c_xm*UY(i-1, j+1)+c_xp*UY(i  , j+1);
        advy =
          -0.5*uy_ym*duxdy_ym
          -0.5*uy_yp*duxdy_yp;
      }
      /* diffusion */
      /* ! diffused in x ! 5 ! */
      double difx;
      {
        const double prefactor = sqrt(Pr)/sqrt(Ra);
        difx = prefactor*(duxdx_xp-duxdx_xm)/DXC(i);
      }
      /* ! diffused in y ! 5 ! */
      double dify;
      {
        const double prefactor = sqrt(Pr)/sqrt(Ra);
        dify = prefactor*(duxdy_yp-duxdy_ym)/dy;
      }
      /* ! pressure gradient ! 6 ! */
      double pre;
      {
        double p_xm = P(i-1, j  );
        double p_xp = P(i  , j  );
        pre = -(p_xp-p_xm)/DXC(i);
      }
      /* ! buoyancy force ! 1 ! */
      double buo = TEMPFORCEX(i, j);
      /* ! summation of explicit terms ! 6 ! */
      SRCUXA(i, j) =
        +advx
        +advy
        +(1.-1.*implicitx)*difx
        +(1.-1.*implicity)*dify
        +buo;
      /* ! summation of implicit terms ! 4 ! */
      SRCUXG(i, j) =
        +(  +1.*implicitx)*difx
        +(  +1.*implicity)*dify
        +pre;
    }
  }
  return 0;
}

