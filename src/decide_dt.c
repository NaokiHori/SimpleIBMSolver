
#include <math.h>
#include <float.h>
#include <mpi.h>
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "temperature.h"
#include "decide_dt.h"
#include "arrays/fluid.h"


/**
 * @brief decide time step size which can integrate the equations stably
 * @param[in] domain : information about domain decomposition and size
 * @param[in] fluid  : velocity
 * @return           : time step size
 */
double decide_dt(const domain_t *domain, const fluid_t *fluid){
  /* ! safety factors for two constraints ! 2 ! */
  const double coef_dt_adv = config.get_double("coef_dt_adv");
  const double coef_dt_dif = config.get_double("coef_dt_dif");
  const double Ra = config.get_double("Ra");
  const double Pr = config.get_double("Pr");
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double *dxc = domain->dxc;
  const double dy   = domain->dy;
  /* ! fluid-induced, temperature-induced diffusivities ! 1 ! */
  const double diffusivity = fmin(sqrt(Ra)/sqrt(Pr), sqrt(Ra)*sqrt(Pr));
  const double *ux = fluid->ux;
  const double *uy = fluid->uy;
  /* ! advective constraint, which will be used as dt ! 25 ! */
  double dt_adv;
  {
    // sufficiently small number
    const double small = 1.e-8;
    // grid-size/velocity
    dt_adv = 1.; // max dt_adv
    // ux
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        // to avoid zero-division
        double vel = fmax(fabs(UX(i, j)), small);
        dt_adv = fmin(dt_adv, DXC(i)/vel);
      }
    }
    // uy
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // to avoid zero-division
        double vel = fmax(fabs(UY(i, j)), small);
        dt_adv = fmin(dt_adv, dy/vel);
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &dt_adv, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt_adv *= coef_dt_adv;
  }
  /* ! diffusive constraints in each direction ! 12 ! */
  double dt_dif_x, dt_dif_y;
  {
    // find minimum grid size in x direction
    double dx = DBL_MAX;
    for(int i = 2; i <= isize; i++){
      dx = fmin(dx, DXC(i));
    }
    dt_dif_x = 0.5/NDIMS*diffusivity*pow(dx, 2.);
    dt_dif_y = 0.5/NDIMS*diffusivity*pow(dy, 2.);
    dt_dif_x *= coef_dt_dif;
    dt_dif_y *= coef_dt_dif;
  }
  /* ! find dt ! 7 ! */
  double dt = dt_adv;
  if(!config.get_bool("implicitx")){
    dt = fmin(dt, dt_dif_x);
  }
  if(!config.get_bool("implicity")){
    dt = fmin(dt, dt_dif_y);
  }
  return dt;
}

