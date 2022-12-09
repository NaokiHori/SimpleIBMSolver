#include <mpi.h>
#include "common.h"
#include "domain.h"
#include "config.h"
#include "internal.h"



/**
 * @brief constructor of the structure
 * @return : structure being allocated and initalised
 */
domain_t *domain_init(void){
  /* ! initialise sdecomp to distribute MPI processes to spatial domain ! 9 ! */
  sdecomp_t *sdecomp = NULL;
  {
    // number of processes in each dimension
    // all zero: let MPI_Cart_create do
    const int dims[NDIMS] = {0, 0};
    // periodicity in each dimension
    const int periods[NDIMS] = {0, 1};
    sdecomp = sdecomp_init(MPI_COMM_WORLD, NDIMS, dims, periods);
  }
  /* ! get resolutions ! 2 ! */
  const int glisize = config.get_int("glisize");
  const int gljsize = config.get_int("gljsize");
  /* ! get domain sizes ! 2 ! */
  const double lx = config.get_double("lx");
  const double ly = config.get_double("ly");
  /* ! allocate and initialise x coordinates ! 4 ! */
  double *xf  = allocate_and_init_xf (glisize, lx);
  double *xc  = allocate_and_init_xc (glisize, xf);
  double *dxf = allocate_and_init_dxf(glisize, xf);
  double *dxc = allocate_and_init_dxc(glisize, xc);
  /* ! grid sizes in homogeneous directions ! 2 ! */
  const double dx = lx / glisize;
  const double dy = ly / gljsize;
  /* ! allocate and assign members ! 21 ! */
  // allocate
  domain_t *domain = common_calloc(1, sizeof(domain_t));
  domain->glsizes = common_calloc(NDIMS, sizeof(   int));
  domain->mysizes = common_calloc(NDIMS, sizeof(   int));
  domain->offsets = common_calloc(NDIMS, sizeof(   int));
  domain->lengths = common_calloc(NDIMS, sizeof(double));
  // assign
  domain->sdecomp = sdecomp;
  domain->glsizes[0] = glisize;
  domain->glsizes[1] = gljsize;
  domain->mysizes[0] = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X1PENCIL, SDECOMP_XDIR, glisize);
  domain->mysizes[1] = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X1PENCIL, SDECOMP_YDIR, gljsize);
  domain->offsets[0] = sdecomp_get_pencil_offset(sdecomp, SDECOMP_X1PENCIL, SDECOMP_XDIR, glisize);
  domain->offsets[1] = sdecomp_get_pencil_offset(sdecomp, SDECOMP_X1PENCIL, SDECOMP_YDIR, gljsize);
  domain->lengths[0] = lx;
  domain->lengths[1] = ly;
  domain->xf = xf;
  domain->xc = xc;
  domain->dxf = dxf;
  domain->dxc = dxc;
  domain->dx = dx;
  domain->dy = dy;
  return domain;
}

