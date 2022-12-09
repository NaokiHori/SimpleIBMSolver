#include <stdbool.h>
#include <math.h>
#include "sdecomp.h"
#include "common.h"
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "linear_system.h"
#include "arrays/fluid.h"
#include "internal.h"


static linear_system_t *linear_system = NULL;


#define DUY(I, J) (duy[ ((J)-1) * (isize) + ((I)-1) ])

/**
 * @brief update uy
 * @param[in   ] domain : information about domain decomposition and size
 * @param[in   ] rkstep : Runge-Kutta step
 * @param[in   ] dt     : time step size
 * @param[inout] fluid  : Runge-Kutta source terms (in), velocity (out)
 * @return              : error code
 */
int fluid_update_velocity_uy(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid){
  /* ! initalise linear solver ! 7 ! */
  if(linear_system == NULL){
    // size of linear system to be solved for "uy": glisize x gljsize
    const int glisize = domain->glsizes[0];
    const int gljsize = domain->glsizes[1];
    const int glsizes_uy[NDIMS] = {glisize, gljsize};
    linear_system = init_linear_system(domain->sdecomp, glsizes_uy);
  }
  /* ! compute increments ! 19 ! */
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const double alpha = RKCOEFS[rkstep].alpha;
    const double beta  = RKCOEFS[rkstep].beta;
    const double gamma = RKCOEFS[rkstep].gamma;
    const double * restrict srcuya = fluid->srcuya;
    const double * restrict srcuyb = fluid->srcuyb;
    const double * restrict srcuyg = fluid->srcuyg;
    double * restrict duy = linear_system->x1pncl;
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        DUY(i, j) =
          +alpha*dt*SRCUYA(i, j)
          +beta *dt*SRCUYB(i, j)
          +gamma*dt*SRCUYG(i, j);
      }
    }
  }
  /* ! pre-factor in front of coefficients ! 7 ! */
  double prefactor;
  {
    const double gamma = RKCOEFS[rkstep].gamma;
    const double Ra = config.get_double("Ra");
    const double Pr = config.get_double("Pr");
    prefactor = (gamma*dt*sqrt(Pr))/(2.*sqrt(Ra));
  }
  /* ! solve linear system in x direction ! 16 ! */
  if(config.get_bool("implicitx")){
    /* set diagonal components of linear system */
    const int isize = domain->glsizes[0];
    for(int i = 1; i <= isize; i++){
      const double * restrict dxf = domain->dxf;
      const double * restrict dxc = domain->dxc;
      double val_l = 1./DXF(i  )/DXC(i  );
      double val_u = 1./DXF(i  )/DXC(i+1);
      double val_c = - val_l - val_u;
      // N.B. i = 1 should be mapped to index 0 for tdm_[luc]
      linear_system->tdm_l[i-1] =    - prefactor * val_l;
      linear_system->tdm_c[i-1] = 1. - prefactor * val_c;
      linear_system->tdm_u[i-1] =    - prefactor * val_u;
    }
    linear_system_solve_in_x(linear_system);
  }
  /* ! transpose x1pencil to y1pencil ! 8 ! */
  const bool needs_x_y_transpose = config.get_bool("implicity");
  if(needs_x_y_transpose){
    sdecomp_transpose_execute(
        linear_system->transposer_x1_to_y1,
        linear_system->x1pncl,
        linear_system->y1pncl
    );
  }
  /* ! solve linear system in y direction ! 15 ! */
  if(config.get_bool("implicity")){
    /* set diagonal components of the linear system */
    const int jsize = domain->glsizes[1];
    for(int j = 1; j <= jsize; j++){
      const double dy = domain->dy;
      double val_l = 1./dy/dy;
      double val_u = 1./dy/dy;
      double val_c = - val_l - val_u;
      // N.B. j = 1 should be mapped to index 0 for tdm_[luc]
      linear_system->tdm_l[j-1] =    - prefactor * val_l;
      linear_system->tdm_c[j-1] = 1. - prefactor * val_c;
      linear_system->tdm_u[j-1] =    - prefactor * val_u;
    }
    linear_system_solve_in_y(linear_system);
  }
  /* ! transpose y1pencil to x1pencil ! 7 ! */
  if(needs_x_y_transpose){
    sdecomp_transpose_execute(
        linear_system->transposer_y1_to_x1,
        linear_system->y1pncl,
        linear_system->x1pncl
    );
  }
  /* ! update ! 12 ! */
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const double * restrict duy = linear_system->x1pncl;
    double * restrict uy = fluid->uy;
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        UY(i, j) += DUY(i, j);
      }
    }
    fluid_update_boundaries_uy(domain, uy);
  }
  return 0;
}

#undef DUY


int fluid_update_velocity_finalise_uy(void){
  linear_system_finalise(linear_system);
  return 0;
}

