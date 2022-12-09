#include "domain.h"
#include "fluid.h"
#include "arrays/fluid.h"
#include "fileio.h"
#include "internal.h"



/**
 * @brief compute total momenta
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : information about domain decomposition and size
 * @param[in] time   : current simulation time
 * @param[in] fluid  : velocity
 * @return           : error code
 */
int check_momentum(const char fname[], const domain_t *domain, const double time, const fluid_t *fluid){
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
  double moms[NDIMS] = {0.};
  /* ! compute total x-momentum ! 8 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize+1; i++){
        double cellsize = DXC(i)*dy*dz;
        moms[0] += UX(i, j, k)*cellsize;
      }
    }
  }
  /* ! compute total y-momentum ! 8 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        double cellsize = DXF(i)*dy*dz;
        moms[1] += UY(i, j, k)*cellsize;
      }
    }
  }
  /* ! compute total z-momentum ! 8 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        double cellsize = DXF(i)*dy*dz;
        moms[2] += UZ(i, j, k)*cellsize;
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, moms, NDIMS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  int myrank;
  MPI_Comm_rank(domain->sdecomp->comm_cart, &myrank);
  if(myrank == 0){
    FILE *fp = fileio_fopen(fname, "a");
    if(fp != NULL){
      fprintf(fp, "%8.2f % 18.15e % 18.15e % 18.15e\n", time, moms[0], moms[1], moms[2]);
      fileio_fclose(fp);
    }
  }
  return 0;
}

