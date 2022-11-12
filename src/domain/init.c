#include <mpi.h>
#include "common.h"
#include "domain.h"
#include "config.h"
#include "internal.h"


#if NDIMS == 2

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

#else // NDIMS == 3

#include <stdio.h>
#include <float.h>
#include <fftw3.h>

static sdecomp_t *optimise_sdecomp_init(void){
  // periodicity in each dimension
  const int periods[NDIMS] = {0, 1, 1};
  // number of processes in each dimension,
  //   which is to be optimised
  int dims_optimum[NDIMS] = {0, 0, 0};
  double wtime_optimum = DBL_MAX;
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  // factorise and decide dims
  for(int n = 1; n <= nprocs; n++){
    if(nprocs % n == 0){
      const int dims[NDIMS] = {1, n, nprocs/n};
      // sanitise, for all dimension, dims should not exceed the number of grid points
      bool legal = true;
      {
        int glsizes[NDIMS];
        glsizes[0] = config.get_int("glisize");
        glsizes[1] = config.get_int("gljsize")/2+1; // after FFT-ed
        glsizes[2] = config.get_int("glksize");
        for(int dim1 = 0; dim1 < NDIMS; dim1++){
          const int glsize = glsizes[dim1];
          for(int dim0 = 0; dim0 < NDIMS; dim0++){
            const int np = dims[dim0];
            if(np > glsize){
              legal = false;
            }
          }
        }
      }
      if(legal){
        // execute transpose and check how long they take
        // for the time being only transposes when solving Poisson equation are considered
        double wtimes[2] = {0., 0.};
        sdecomp_t *sdecomp = sdecomp_init(MPI_COMM_WORLD, NDIMS, dims, periods);
        int glsizes[NDIMS];
        glsizes[0] = config.get_int("glisize");
        glsizes[1] = config.get_int("gljsize");
        glsizes[2] = config.get_int("glksize");
        int sizes[NDIMS];
        // x1 to y1 (real)
        sizes[0] = glsizes[0];
        sizes[1] = glsizes[1];
        sizes[2] = glsizes[2];
        sdecomp_transpose_t *x1_to_y1_r = sdecomp_transpose_fwrd_init(sdecomp, SDECOMP_X1PENCIL, sizes, sizeof(      double),           MPI_DOUBLE);
        sdecomp_transpose_t *y1_to_x1_r = sdecomp_transpose_bwrd_init(sdecomp, SDECOMP_Y1PENCIL, sizes, sizeof(      double),           MPI_DOUBLE);
        // y1 to z1 (complex)
        sizes[0] = glsizes[0];
        sizes[1] = glsizes[1]/2+1;
        sizes[2] = glsizes[2];
        sdecomp_transpose_t *y1_to_z1_c = sdecomp_transpose_fwrd_init(sdecomp, SDECOMP_Y1PENCIL, sizes, sizeof(fftw_complex), MPI_C_DOUBLE_COMPLEX);
        sdecomp_transpose_t *z1_to_y1_c = sdecomp_transpose_bwrd_init(sdecomp, SDECOMP_Z1PENCIL, sizes, sizeof(fftw_complex), MPI_C_DOUBLE_COMPLEX);
        // z1 to x2 (complex)
        sizes[0] = glsizes[0];
        sizes[1] = glsizes[1]/2+1;
        sizes[2] = glsizes[2];
        sdecomp_transpose_t *z1_to_x2_c = sdecomp_transpose_fwrd_init(sdecomp, SDECOMP_Z1PENCIL, sizes, sizeof(fftw_complex), MPI_C_DOUBLE_COMPLEX);
        sdecomp_transpose_t *x2_to_z1_c = sdecomp_transpose_bwrd_init(sdecomp, SDECOMP_X2PENCIL, sizes, sizeof(fftw_complex), MPI_C_DOUBLE_COMPLEX);
        // execute transpose
        wtimes[0] = common_get_wtime();
        sdecomp_test_transpose_3d_x1_to_y1(sdecomp, x1_to_y1_r);
        sdecomp_test_transpose_3d_y1_to_z1(sdecomp, y1_to_z1_c);
        sdecomp_test_transpose_3d_z1_to_x2(sdecomp, z1_to_x2_c);
        sdecomp_test_transpose_3d_x2_to_z1(sdecomp, x2_to_z1_c);
        sdecomp_test_transpose_3d_z1_to_y1(sdecomp, z1_to_y1_c);
        sdecomp_test_transpose_3d_y1_to_x1(sdecomp, y1_to_x1_r);
        wtimes[1] = common_get_wtime();
        // clean-up tentative transpose plans
        sdecomp_transpose_finalise(x1_to_y1_r);
        sdecomp_transpose_finalise(y1_to_z1_c);
        sdecomp_transpose_finalise(z1_to_x2_c);
        sdecomp_transpose_finalise(x2_to_z1_c);
        sdecomp_transpose_finalise(z1_to_y1_c);
        sdecomp_transpose_finalise(y1_to_x1_r);
        // clean-up current sdecomp config
        sdecomp_finalise(sdecomp);
        // check time
        double wtime = wtimes[1] - wtimes[0];
        if(wtime < wtime_optimum){
          dims_optimum[0] = dims[0];
          dims_optimum[1] = dims[1];
          dims_optimum[2] = dims[2];
          wtime_optimum = wtime;
        }
        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        if(myrank == 0){
          const int nd = common_get_ndigits(nprocs);
          printf("dims: [%*d, %*d, %*d]: % .7e [sec]\n", nd, dims[0], nd, dims[1], nd, dims[2], wtime);
        }
      }
    }
  }
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if(myrank == 0){
    const int nd = common_get_ndigits(nprocs);
    printf("Conclusive domain decomposition: [%*d, %*d, %*d]\n", nd, dims_optimum[0], nd, dims_optimum[1], nd, dims_optimum[2]);
  }
  // create sdecomp which will be used in the main run
  sdecomp_t *sdecomp = sdecomp_init(MPI_COMM_WORLD, NDIMS, dims_optimum, periods);
  return sdecomp;
}

/**
 * @brief constructor of the structure
 * @return : structure being allocated and initalised
 */
domain_t *domain_init(void){
  /* ! initialise sdecomp to distribute MPI processes to spatial domain ! 1 ! */
  sdecomp_t *sdecomp = optimise_sdecomp_init();
  /* ! get resolutions ! 3 ! */
  const int glisize = config.get_int("glisize");
  const int gljsize = config.get_int("gljsize");
  const int glksize = config.get_int("glksize");
  /* ! get domain sizes ! 3 ! */
  const double lx = config.get_double("lx");
  const double ly = config.get_double("ly");
  const double lz = config.get_double("lz");
  /* ! allocate and initialise x coordinates ! 4 ! */
  double *xf  = allocate_and_init_xf (glisize, lx);
  double *xc  = allocate_and_init_xc (glisize, xf);
  double *dxf = allocate_and_init_dxf(glisize, xf);
  double *dxc = allocate_and_init_dxc(glisize, xc);
  /* ! grid sizes in homogeneous directions ! 3 ! */
  const double dx = lx / glisize;
  const double dy = ly / gljsize;
  const double dz = lz / glksize;
  /* ! allocate and assign members ! 26 ! */
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
  domain->glsizes[2] = glksize;
  domain->mysizes[0] = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X1PENCIL, SDECOMP_XDIR, glisize);
  domain->mysizes[1] = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X1PENCIL, SDECOMP_YDIR, gljsize);
  domain->mysizes[2] = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X1PENCIL, SDECOMP_ZDIR, glksize);
  domain->offsets[0] = sdecomp_get_pencil_offset(sdecomp, SDECOMP_X1PENCIL, SDECOMP_XDIR, glisize);
  domain->offsets[1] = sdecomp_get_pencil_offset(sdecomp, SDECOMP_X1PENCIL, SDECOMP_YDIR, gljsize);
  domain->offsets[2] = sdecomp_get_pencil_offset(sdecomp, SDECOMP_X1PENCIL, SDECOMP_ZDIR, glksize);
  domain->lengths[0] = lx;
  domain->lengths[1] = ly;
  domain->lengths[2] = lz;
  domain->xf = xf;
  domain->xc = xc;
  domain->dxf = dxf;
  domain->dxc = dxc;
  domain->dx = dx;
  domain->dy = dy;
  domain->dz = dz;
  return domain;
}

#endif // NDIMS
