#include <stdio.h>
#include <math.h>
#include "domain.h"
#include "fluid.h"
#include "fileio.h"
#include "arrays/fluid.h"
#include "internal.h"


#if NDIMS == 2

/**
 * @brief compute local divergence and write the maximum value
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : domain information
 * @param[in] time   : current simulation time
 * @param[in] fluid  : velocity
 * @return           : error code
 */
int check_divergence(const char fname[], const domain_t *domain, const double time, const fluid_t *fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double *dxf = domain->dxf;
  const double dy = domain->dy;
  const double *ux = fluid->ux;
  const double *uy = fluid->uy;
  double divmax = 0.;
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      /* ! compute local divergence ! 7 ! */
      double ux_xm = UX(i  , j  );
      double ux_xp = UX(i+1, j  );
      double uy_ym = UY(i  , j  );
      double uy_yp = UY(i  , j+1);
      double div =
        (ux_xp-ux_xm)/DXF(i)
       +(uy_yp-uy_ym)/dy;
      /* ! check local maximum divergence ! 1 ! */
      divmax = fmax(divmax, fabs(div));
    }
  }
  /* ! collect information among all processes ! 1 ! */
  MPI_Allreduce(MPI_IN_PLACE, &divmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  /* ! information is dumped, similar for other functions ! 9 ! */
  int myrank;
  MPI_Comm_rank(domain->sdecomp->comm_cart, &myrank);
  if(myrank == 0){
    FILE *fp = fileio_fopen(fname, "a");
    if(fp != NULL){
      fprintf(fp, "%8.2f % .1e\n", time, divmax);
      fileio_fclose(fp);
    }
  }
  return 0;
}

#else // NDIMS == 3

/**
 * @brief compute local divergence and write the maximum value
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : domain information
 * @param[in] time   : current simulation time
 * @param[in] fluid  : velocity
 * @return           : error code
 */
int check_divergence(const char fname[], const domain_t *domain, const double time, const fluid_t *fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double *dxf = domain->dxf;
  const double dy = domain->dy;
  const double dz = domain->dz;
  const double *ux = fluid->ux;
  const double *uy = fluid->uy;
  const double *uz = fluid->uz;
  double divmax = 0.;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        /* ! compute local divergence ! 10 ! */
        double ux_xm = UX(i  , j  , k  );
        double ux_xp = UX(i+1, j  , k  );
        double uy_ym = UY(i  , j  , k  );
        double uy_yp = UY(i  , j+1, k  );
        double uz_zm = UZ(i  , j  , k  );
        double uz_zp = UZ(i  , j  , k+1);
        double div =
          (ux_xp-ux_xm)/DXF(i)
         +(uy_yp-uy_ym)/dy
         +(uz_zp-uz_zm)/dz;
        /* ! check local maximum divergence ! 1 ! */
        divmax = fmax(divmax, fabs(div));
      }
    }
  }
  /* ! collect information among all processes ! 1 ! */
  MPI_Allreduce(MPI_IN_PLACE, &divmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  /* ! information are output, similar for other functions ! 9 ! */
  int myrank;
  MPI_Comm_rank(domain->sdecomp->comm_cart, &myrank);
  if(myrank == 0){
    FILE *fp = fileio_fopen(fname, "a");
    if(fp != NULL){
      fprintf(fp, "%8.2f % .1e\n", time, divmax);
      fileio_fclose(fp);
    }
  }
  return 0;
}

#endif // NDIMS
