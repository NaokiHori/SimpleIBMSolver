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
 * @brief comute right-hand-side of Runge-Kutta scheme of uy
 * @param[in   ] domain : information related to MPI domain decomposition
 * @param[in   ] rkstep : Runge-Kutta step
 * @param[inout] fluid  : velocity and pressure (in), RK source terms (inout)
 * @return              : error code
 */
int compute_rhs_uy(const domain_t * restrict domain, const int rkstep, fluid_t * restrict fluid){
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
  const double * restrict p = fluid->p;
  double * restrict srcuya = fluid->srcuya;
  double * restrict srcuyb = fluid->srcuyb;
  double * restrict srcuyg = fluid->srcuyg;
  /* ! previous k-step source term is copied ! 3 ! */
  if(rkstep != 0){
    memcpy(srcuyb, srcuya, SRCUYA_SIZE_0 * SRCUYA_SIZE_1 * SRCUYA_SIZE_2 * sizeof(double));
  }
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        /* ! velocity-gradient tensor L_yx ! 2 ! */
        double duydx_xm = (-UY(i-1, j  , k  )+UY(i  , j  , k  ))/DXC(i  );
        double duydx_xp = (-UY(i  , j  , k  )+UY(i+1, j  , k  ))/DXC(i+1);
        /* ! velocity-gradient tensor L_yy ! 2 ! */
        double duydy_ym = (-UY(i  , j-1, k  )+UY(i  , j  , k  ))/dy;
        double duydy_yp = (-UY(i  , j  , k  )+UY(i  , j+1, k  ))/dy;
        /* ! velocity-gradient tensor L_yz ! 2 ! */
        double duydz_zm = (-UY(i  , j  , k-1)+UY(i  , j  , k  ))/dz;
        double duydz_zp = (-UY(i  , j  , k  )+UY(i  , j  , k+1))/dz;
        /* advection */
        /* ! transported by ux ! 10 ! */
        double advx;
        {
          double c_xm = DXC(i  )/(2.*DXF(i));
          double c_xp = DXC(i+1)/(2.*DXF(i));
          double ux_xm = 0.5*UX(i  , j-1, k  )+0.5*UX(i  , j  , k  );
          double ux_xp = 0.5*UX(i+1, j-1, k  )+0.5*UX(i+1, j  , k  );
          advx =
            -c_xm*ux_xm*duydx_xm
            -c_xp*ux_xp*duydx_xp;
        }
        /* ! transported by uy ! 8 ! */
        double advy;
        {
          double uy_ym = 0.5*UY(i  , j-1, k  )+0.5*UY(i  , j  , k  );
          double uy_yp = 0.5*UY(i  , j  , k  )+0.5*UY(i  , j+1, k  );
          advy =
            -0.5*uy_ym*duydy_ym
            -0.5*uy_yp*duydy_yp;
        }
        /* ! transported by uz ! 8 ! */
        double advz;
        {
          double uz_zm = 0.5*UZ(i  , j-1, k  )+0.5*UZ(i  , j  , k  );
          double uz_zp = 0.5*UZ(i  , j-1, k+1)+0.5*UZ(i  , j  , k+1);
          advz =
            -0.5*uz_zm*duydz_zm
            -0.5*uz_zp*duydz_zp;
        }
        /* diffusion */
        /* ! diffused in x ! 5 ! */
        double difx;
        {
          const double prefactor = sqrt(Pr)/sqrt(Ra);
          difx = prefactor*(duydx_xp-duydx_xm)/DXF(i);
        }
        /* ! diffused in y ! 5 ! */
        double dify;
        {
          const double prefactor = sqrt(Pr)/sqrt(Ra);
          dify = prefactor*(duydy_yp-duydy_ym)/dy;
        }
        /* ! diffused in z ! 5 ! */
        double difz;
        {
          const double prefactor = sqrt(Pr)/sqrt(Ra);
          difz = prefactor*(duydz_zp-duydz_zm)/dz;
        }
        /* ! pressure gradient ! 6 ! */
        double pre;
        {
          double p_ym = P(i  , j-1, k  );
          double p_yp = P(i  , j  , k  );
          pre = -(p_yp-p_ym)/dy;
        }
        /* ! summation of explicit terms ! 7 ! */
        SRCUYA(i, j, k) =
          +advx
          +advy
          +advz
          +(1.-1.*implicitx)*difx
          +(1.-1.*implicity)*dify
          +(1.-1.*implicitz)*difz;
        /* ! summation of implicit terms ! 5 ! */
        SRCUYG(i, j, k) =
          +(  +1.*implicitx)*difx
          +(  +1.*implicity)*dify
          +(  +1.*implicitz)*difz
          +pre;
      }
    }
  }
  return 0;
}

