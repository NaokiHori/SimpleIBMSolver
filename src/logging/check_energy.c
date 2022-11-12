#include <stdio.h>
#include <math.h>
#include "domain.h"
#include "fluid.h"
#include "temperature.h"
#include "arrays/fluid.h"
#include "arrays/temperature.h"
#include "fileio.h"
#include "internal.h"


#if NDIMS == 2

/**
 * @brief compute total kinetic and thermal energies
 * @param[in] fname       : file name to which the log is written
 * @param[in] domain      : information related to MPI domain decomposition
 * @param[in] time        : current simulation time
 * @param[in] fluid       : velocity
 * @param[in] temperature : temperature
 * @return                : error code
 */
int check_energy(const char fname[], const domain_t *domain, const double time, const fluid_t *fluid, const temperature_t *temperature){
  /* compute kinetic and thermal energies */
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double *dxf = domain->dxf;
  const double *dxc = domain->dxc;
  const double dy = domain->dy;
  const double *ux = fluid->ux;
  const double *uy = fluid->uy;
  const double *temp = temperature->temp;
  // velocity in each dimension and thermal
  double quantities[NDIMS+1] = {0.};
  /* ! compute quadratic quantity in x direction ! 6 ! */
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize+1; i++){
      double cellsize = DXC(i)*dy;
      quantities[0] += 0.5*pow(UX(i, j), 2.)*cellsize;
    }
  }
  /* ! compute quadratic quantity in y direction ! 6 ! */
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      double cellsize = DXF(i)*dy;
      quantities[1] += 0.5*pow(UY(i, j), 2.)*cellsize;
    }
  }
  /* ! compute thermal energy ! 6 ! */
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      double cellsize = DXF(i)*dy;
      quantities[2] += 0.5*pow(TEMP(i, j), 2.)*cellsize;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, quantities, NDIMS+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  int myrank;
  MPI_Comm_rank(domain->sdecomp->comm_cart, &myrank);
  if(myrank == 0){
    FILE *fp = fileio_fopen(fname, "a");
    if(fp != NULL){
      fprintf(fp, "%8.2f % 18.15e % 18.15e % 18.15e\n", time, quantities[0], quantities[1], quantities[2]);
      fileio_fclose(fp);
    }
  }
  return 0;
}

#else // NDIMS == 3

/**
 * @brief compute total kinetic and thermal energies
 * @param[in] fname       : file name to which the log is written
 * @param[in] domain      : information related to MPI domain decomposition
 * @param[in] time        : current simulation time
 * @param[in] fluid       : velocity
 * @param[in] temperature : temperature
 * @return                : error code
 */
int check_energy(const char fname[], const domain_t *domain, const double time, const fluid_t *fluid, const temperature_t *temperature){
  /* compute kinetic and thermal energies */
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double *dxf = domain->dxf;
  const double *dxc = domain->dxc;
  const double dy = domain->dy;
  const double dz = domain->dz;
  const double *ux = fluid->ux;
  const double *uy = fluid->uy;
  const double *uz = fluid->uz;
  const double *temp = temperature->temp;
  // velocity in each dimension and thermal
  double quantities[NDIMS+1] = {0.};
  /* ! compute quadratic quantity in x direction ! 8 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize+1; i++){
        double cellsize = DXC(i)*dy*dz;
        quantities[0] += 0.5*pow(UX(i, j, k), 2.)*cellsize;
      }
    }
  }
  /* ! compute quadratic quantity in y direction ! 8 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        double cellsize = DXF(i)*dy*dz;
        quantities[1] += 0.5*pow(UY(i, j, k), 2.)*cellsize;
      }
    }
  }
  /* ! compute quadratic quantity in z direction ! 8 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        double cellsize = DXF(i)*dy*dz;
        quantities[2] += 0.5*pow(UZ(i, j, k), 2.)*cellsize;
      }
    }
  }
  /* ! compute thermal energy ! 8 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        double cellsize = DXF(i)*dy*dz;
        quantities[3] += 0.5*pow(TEMP(i, j, k), 2.)*cellsize;
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, quantities, NDIMS+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  int myrank;
  MPI_Comm_rank(domain->sdecomp->comm_cart, &myrank);
  if(myrank == 0){
    FILE *fp = fileio_fopen(fname, "a");
    if(fp != NULL){
      fprintf(fp, "%8.2f % 18.15e % 18.15e % 18.15e % 18.15e\n", time, quantities[0], quantities[1], quantities[2], quantities[3]);
      fileio_fclose(fp);
    }
  }
  return 0;
}

#endif // NDIMS
