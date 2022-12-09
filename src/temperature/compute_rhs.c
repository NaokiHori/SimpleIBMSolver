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
  const int implicitz = config.get_bool("implicitz");
  const double Ra = config.get_double("Ra");
  const double Pr = config.get_double("Pr");
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict dxf = domain->dxf;
  const double * restrict dxc = domain->dxc;
  const double dy = domain->dy;
  const double dz = domain->dz;
  const double * restrict ux = fluid->ux;
  const double * restrict uy = fluid->uy;
  const double * restrict uz = fluid->uz;
  const double * restrict temp = temperature->temp;
  double * restrict srctempa = temperature->srctempa;
  double * restrict srctempb = temperature->srctempb;
  double * restrict srctempg = temperature->srctempg;
  /* ! previous k-step source term of temp is copied ! 3 ! */
  if(rkstep != 0){
    memcpy(srctempb, srctempa, SRCTEMPA_SIZE_0 * SRCTEMPA_SIZE_1 * SRCTEMPA_SIZE_2 * sizeof(double));
  }
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        /* ! gradient of T in x ! 2 ! */
        double dtdx_xm = (-TEMP(i-1, j  , k  )+TEMP(i  , j  , k  ))/DXC(i  );
        double dtdx_xp = (-TEMP(i  , j  , k  )+TEMP(i+1, j  , k  ))/DXC(i+1);
        /* ! gradient of T in y ! 2 ! */
        double dtdy_ym = (-TEMP(i  , j-1, k  )+TEMP(i  , j  , k  ))/dy;
        double dtdy_yp = (-TEMP(i  , j  , k  )+TEMP(i  , j+1, k  ))/dy;
        /* ! gradient of T in z ! 2 ! */
        double dtdz_zm = (-TEMP(i  , j  , k-1)+TEMP(i  , j  , k  ))/dz;
        double dtdz_zp = (-TEMP(i  , j  , k  )+TEMP(i  , j  , k+1))/dz;
        /* advection */
        /* ! T is transported by ux ! 10 ! */
        double advx;
        {
          double c_xm = DXC(i  )/(2.*DXF(i));
          double c_xp = DXC(i+1)/(2.*DXF(i));
          double ux_xm = UX(i  , j  , k  );
          double ux_xp = UX(i+1, j  , k  );
          advx =
            -c_xm*ux_xm*dtdx_xm
            -c_xp*ux_xp*dtdx_xp;
        }
        /* ! T is transported by uy ! 8 ! */
        double advy;
        {
          double uy_ym = UY(i  , j  , k  );
          double uy_yp = UY(i  , j+1, k  );
          advy =
            -0.5*uy_ym*dtdy_ym
            -0.5*uy_yp*dtdy_yp;
        }
        /* ! T is transported by uz ! 8 ! */
        double advz;
        {
          double uz_zm = UZ(i  , j  , k  );
          double uz_zp = UZ(i  , j  , k+1);
          advz =
            -0.5*uz_zm*dtdz_zm
            -0.5*uz_zp*dtdz_zp;
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
        /* ! T is diffused in z ! 5 ! */
        double difz;
        {
          const double prefactor = 1./sqrt(Pr)/sqrt(Ra);
          difz = prefactor*(dtdz_zp-dtdz_zm)/dz;
        }
        /* summation */
        /* ! summation of explicit terms ! 7 ! */
        SRCTEMPA(i, j, k) =
          +advx
          +advy
          +advz
          +(1.-1.*implicitx)*difx
          +(1.-1.*implicity)*dify
          +(1.-1.*implicitz)*difz;
        /* ! summation of implicit terms ! 4 ! */
        SRCTEMPG(i, j, k) =
          +(  +1.*implicitx)*difx
          +(  +1.*implicity)*dify
          +(  +1.*implicitz)*difz;
      }
    }
  }
  return 0;
}

