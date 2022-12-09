#include <stdio.h>
#include <math.h>
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "temperature.h"
#include "arrays/fluid.h"
#include "arrays/temperature.h"
#include "fileio.h"
#include "internal.h"



/**
 * @brief compute Nusselt number based on heat flux on the walls
 * @param[in] domain      : information related to MPI domain decomposition
 * @param[in] temperature : temperature
 * @return                : Nusselt number
 */
static double compute_nu_wall(const domain_t *domain, const temperature_t *temperature){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double ly = domain->lengths[1];
  const double *dxc = domain->dxc;
  const double dy = domain->dy;
  const double *temp = temperature->temp;
  /* ! heat flux on the walls ! 10 ! */
  double retval = 0.;
  // integral on the walls
  for(int j = 1; j <= jsize; j++){
    double dtdx0 = (TEMP(      1, j)-TEMP(    0, j))/DXC(      1);
    double dtdx1 = (TEMP(isize+1, j)-TEMP(isize, j))/DXC(isize+1);
    retval -= 0.5*(dtdx0+dtdx1)*dy;
  }
  MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // average on the walls
  retval /= ly;
  return retval;
}

/**
 * @brief compute Nusselt number based on total energy injection
 * @param[in] domain      : information related to MPI domain decomposition
 * @param[in] fluid       : wall-normal velocity (ux)
 * @param[in] temperature : temperature
 * @return                : Nusselt number
 */
static double compute_nu_inje(const domain_t *domain, const fluid_t *fluid, const temperature_t *temperature){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double lx = domain->lengths[0];
  const double ly = domain->lengths[1];
  const double *dxc = domain->dxc;
  const double dy = domain->dy;
  const double Ra = config.get_double("Ra");
  const double Pr = config.get_double("Pr");
  const double *ux = fluid->ux;
  const double *temp = temperature->temp;
  /* ! energy injection ! 13 ! */
  double retval = 0.;
  // integral in the whole domain
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize+1; i++){
      double integrand = UX(i, j)*0.5*(TEMP(i-1, j)+TEMP(i, j));
      double cellsize = DXC(i)*dy;
      retval += sqrt(Pr)*sqrt(Ra)*integrand*cellsize;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // average in the whole space
  retval /= lx*ly;
  retval += 1.;
  return retval;
}

/**
 * @brief compute Nusselt number based on kinetic dissipation
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] fluid  : velocity
 * @return           : Nusselt number
 */
static double compute_nu_disk(const domain_t *domain, const fluid_t *fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double lx = domain->lengths[0];
  const double ly = domain->lengths[1];
  const double *dxf = domain->dxf;
  const double *dxc = domain->dxc;
  const double dy = domain->dy;
  const double Ra = config.get_double("Ra");
  const double Pr = config.get_double("Pr");
  const double *ux = fluid->ux;
  const double *uy = fluid->uy;
  double retval = 0.;
  /* ! x-momentum contribution ! 35 ! */
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize+1; i++){
      double cellsize = DXC(i)*dy;
      // local dissipation at UX(i, j)
      double localval = 0.;
      {
        // x-negative contribution
        if(i != 1){
          double dux   = UX(i  , j  )-UX(i-1, j  );
          double duxdx = dux/DXF(i-1);
          localval += 0.5/DXC(i)*dux*duxdx;
        }
        // x-positive contribution
        if(i != isize+1){
          double dux   = UX(i+1, j  )-UX(i  , j  );
          double duxdx = dux/DXF(i  );
          localval += 0.5/DXC(i)*dux*duxdx;
        }
        // y-negative contribution
        {
          double dux   = UX(i  , j  )-UX(i  , j-1);
          double duxdy = dux/dy;
          localval += 0.5/dy*dux*duxdy;
        }
        // y-positive contribution
        {
          double dux   = UX(i  , j+1)-UX(i  , j  );
          double duxdy = dux/dy;
          localval += 0.5/dy*dux*duxdy;
        }
      }
      localval *= sqrt(Pr)/sqrt(Ra);
      retval += sqrt(Pr)*sqrt(Ra)*localval*cellsize;
    }
  }
  /* ! y-momentum contribution ! 37 ! */
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      double cellsize = DXF(i)*dy;
      // local dissipation at UY(i, j)
      double localval = 0.;
      {
        // x-negative contribution
        {
          double c = i ==    1 ? 1. : 0.5;
          double duy   = UY(i  , j  )-UY(i-1, j  );
          double duydx = duy/DXC(i  );
          localval += c/DXF(i)*duy*duydx;
        }
        // x-positive contribution
        {
          double c = i == isize ? 1. : 0.5;
          double duy   = UY(i+1, j  )-UY(i  , j  );
          double duydx = duy/DXC(i+1);
          localval += c/DXF(i)*duy*duydx;
        }
        // y-negative contribution
        {
          double duy   = UY(i  , j  )-UY(i  , j-1);
          double duydy = duy/dy;
          localval += 0.5/dy*duy*duydy;
        }
        // y-positive contribution
        {
          double duy   = UY(i  , j+1)-UY(i  , j  );
          double duydy = duy/dy;
          localval += 0.5/dy*duy*duydy;
        }
      }
      localval *= sqrt(Pr)/sqrt(Ra);
      retval += sqrt(Pr)*sqrt(Ra)*localval*cellsize;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // average in the whole domain
  retval /= lx*ly;
  retval += 1.;
  return retval;
}

/**
 * @brief compute Nusselt number based on thermal dissipation
 * @param[in] domain      : information related to MPI domain decomposition
 * @param[in] temperature : temperature
 * @return                : Nusselt number
 */
static double compute_nu_dish(const domain_t *domain, const temperature_t *temperature){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double lx = domain->lengths[0];
  const double ly = domain->lengths[1];
  const double *dxf = domain->dxf;
  const double *dxc = domain->dxc;
  const double dy = domain->dy;
  const double Ra = config.get_double("Ra");
  const double Pr = config.get_double("Pr");
  const double *temp = temperature->temp;
  double retval = 0.;
  /* ! thermal energy dissipation rate ! 41 ! */
  // integral in the whole domain
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      double cellsize = DXF(i)*dy;
      // local dissipation at TEMP(i, j)
      double localval = 0.;
      {
        // x-negative contribution
        {
          double c = i ==    1 ? 1. : 0.5;
          double dtemp   = TEMP(i  , j  )-TEMP(i-1, j  );
          double dtempdx = dtemp/DXC(i  );
          localval += c/DXF(i)*dtemp*dtempdx;
        }
        // x-positive contribution
        {
          double c = i == isize ? 1. : 0.5;
          double dtemp   = TEMP(i+1, j  )-TEMP(i  , j  );
          double dtempdx = dtemp/DXC(i+1);
          localval += c/DXF(i)*dtemp*dtempdx;
        }
        // y-negative contribution
        {
          double dtemp   = TEMP(i  , j  )-TEMP(i  , j-1);
          double dtempdy = dtemp/dy;
          localval += 0.5/dy*dtemp*dtempdy;
        }
        // y-positive contribution
        {
          double dtemp   = TEMP(i  , j+1)-TEMP(i  , j  );
          double dtempdy = dtemp/dy;
          localval += 0.5/dy*dtemp*dtempdy;
        }
      }
      localval *= 1./sqrt(Pr)/sqrt(Ra);
      retval += sqrt(Pr)*sqrt(Ra)*localval*cellsize;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // average in the whole domain
  retval /= lx*ly;
  return retval;
}


/**
 * @brief compute Nusselt number based on various definitions
 * @param[in] fname       : file name to which the log is written
 * @param[in] domain      : information related to MPI domain decomposition
 * @param[in] time        : current simulation time
 * @param[in] fluid       : velocity
 * @param[in] temperature : temperature
 * @return                : error code
 */
int check_nusselt(const char fname[], const domain_t *domain, const double time, const fluid_t *fluid, const temperature_t *temperature){
  /* compute temperature Nusselt number based on several definintions */
  // compute Nu from heat flux on the walls
  double nu_wall = compute_nu_wall(domain, temperature);
  // compute Nu from energy injection
  double nu_inje = compute_nu_inje(domain, fluid, temperature);
  // compute Nu from kinetic energy dissipation rate
  double nu_eps_k = compute_nu_disk(domain, fluid);
  // comoute Nu from thermal energy dissipation rate
  double nu_eps_h = compute_nu_dish(domain, temperature);
  int myrank;
  MPI_Comm_rank(domain->sdecomp->comm_cart, &myrank);
  if(myrank == 0){
    FILE *fp = fileio_fopen(fname, "a");
    if(fp != NULL){
      fprintf(fp, "%8.2f % 18.15e % 18.15e % 18.15e % 18.15e\n", time, nu_wall, nu_inje, nu_eps_k, nu_eps_h);
      fileio_fclose(fp);
    }
  }
  return 0;
}

