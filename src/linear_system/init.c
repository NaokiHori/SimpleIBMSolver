#include <mpi.h>
#include "common.h"
#include "sdecomp.h"
#include "domain.h"
#include "linear_system.h"


static int max2(const int a, const int b){
  if(a > b){
    return a;
  }else{
    return b;
  }
}


static int max3(const int a, const int b, const int c){
  return max2(a, max2(b, c));
}

/**
 * @brief initialise linear solver to update field implicitly
 * @param[in] sdecomp : information about MPI decomposition
 * @param[in] glsizes : GLOBAL size of ARRAY (can differ for ux, uy, uz, p)
 * @return            : structure storing buffers and plans to solve linear systems in each direction
 */
linear_system_t *init_linear_system(sdecomp_t * restrict sdecomp, const int glsizes[NDIMS]){
  // allocate main structure
  linear_system_t *linear_system = common_calloc(1, sizeof(linear_system_t));
  linear_system->sdecomp = sdecomp;
  linear_system->glsizes = common_calloc(NDIMS, sizeof(int));
  for(int dim = 0; dim < NDIMS; dim++){
    linear_system->glsizes[dim] = glsizes[dim];
  }
  // structure storing information of MPI decomposition
  // LOCAL array sizes for all pencils
  const int x1pncl_isize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X1PENCIL, SDECOMP_XDIR, glsizes[0]);
  const int x1pncl_jsize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X1PENCIL, SDECOMP_YDIR, glsizes[1]);
  const int x1pncl_ksize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X1PENCIL, SDECOMP_ZDIR, glsizes[2]);
  const int y1pncl_isize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Y1PENCIL, SDECOMP_XDIR, glsizes[0]);
  const int y1pncl_jsize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Y1PENCIL, SDECOMP_YDIR, glsizes[1]);
  const int y1pncl_ksize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Y1PENCIL, SDECOMP_ZDIR, glsizes[2]);
  const int z1pncl_isize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Z1PENCIL, SDECOMP_XDIR, glsizes[0]);
  const int z1pncl_jsize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Z1PENCIL, SDECOMP_YDIR, glsizes[1]);
  const int z1pncl_ksize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Z1PENCIL, SDECOMP_ZDIR, glsizes[2]);
  /* allocate pencils */
  {
    linear_system->x1pncl = common_calloc(x1pncl_isize * x1pncl_jsize * x1pncl_ksize, sizeof(double));
    linear_system->y1pncl = common_calloc(y1pncl_isize * y1pncl_jsize * y1pncl_ksize, sizeof(double));
    linear_system->z1pncl = common_calloc(z1pncl_isize * z1pncl_jsize * z1pncl_ksize, sizeof(double));
  }
  /* parallel matrix transpose */
  {
    linear_system->transposer_x1_to_y1 = sdecomp_transpose_fwrd_init(sdecomp, SDECOMP_X1PENCIL, glsizes, sizeof(double), MPI_DOUBLE);
    linear_system->transposer_y1_to_x1 = sdecomp_transpose_bwrd_init(sdecomp, SDECOMP_Y1PENCIL, glsizes, sizeof(double), MPI_DOUBLE);
    linear_system->transposer_y1_to_z1 = sdecomp_transpose_fwrd_init(sdecomp, SDECOMP_Y1PENCIL, glsizes, sizeof(double), MPI_DOUBLE);
    linear_system->transposer_z1_to_y1 = sdecomp_transpose_bwrd_init(sdecomp, SDECOMP_Z1PENCIL, glsizes, sizeof(double), MPI_DOUBLE);
  }
  /* allocate tri-diagonal matrix solver */
  {
    // x, y, z tri-diagonal matrices are NOT treated simultaneously
    // so allocate buffer for the largest matrix and share buffers
    const int size = max3(x1pncl_isize, y1pncl_jsize, z1pncl_ksize);
    linear_system->tdm_l = common_calloc(size, sizeof(double));
    linear_system->tdm_c = common_calloc(size, sizeof(double));
    linear_system->tdm_u = common_calloc(size, sizeof(double));
  }
  return linear_system;
}

